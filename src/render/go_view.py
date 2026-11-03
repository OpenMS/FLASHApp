"""The GO enrichment cell of the Viewer, for a run and for a comparison.

The cell reads the proteins (run) or entries (comparison) picked in the other
cells of the grid from their shared StateManager (see src/render/go_terms.py)
and shows the Gene Ontology terms over-represented in them: a bar chart of the
most significant terms and a table of all of them.
"""

import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
import streamlit as st

from openms_insight import Table
from openms_insight.core.cache import CacheMissError

from src.common.common import show_fig
from src.render import go_terms

# Terms in the bar chart
TOP_TERMS = 15
# Results kept per session (one per selection)
_RESULTS_KEPT = 8

RUN_HINT = (
    'Select proteins to test their GO terms: tick them in the protein, proteoform or '
    'PrSM tables (the box in the header ticks every row shown, Shift-click ticks a range), '
    'or drag a box or lasso on the proteoform heatmap after choosing **Select** in its toolbar.'
)
COMPARISON_HINT = (
    'Select {nouns} to test their GO terms: drag a box or lasso on the volcano plot, or use '
    'its **All changed**, **Higher** and **Lower** buttons; tick rows in the comparison table; '
    'or select rows on the clustered heatmap. Shift adds to the selection.'
)

TERM_COLUMNS = [
    {'title': '#', 'field': 'rank', 'sorter': 'number', 'width': 50},
    {'title': 'GO term', 'field': 'term', 'sorter': 'string', 'minWidth': 220, 'hozAlign': 'left'},
    {'title': 'Selected', 'field': 'count', 'sorter': 'number',
     'headerTooltip': 'Selected proteins annotated with the term.'},
    {'title': 'Background', 'field': 'background', 'sorter': 'number',
     'headerTooltip': 'Proteins of the background annotated with the term.'},
    {'title': 'Ratio', 'field': 'ratio', 'sorter': 'string',
     'headerTooltip': 'Selected proteins with the term / selected proteins with any GO term.'},
    {'title': 'Bg ratio', 'field': 'bg_ratio', 'sorter': 'string',
     'headerTooltip': 'Background proteins with the term / background proteins with any GO term.'},
    {'title': 'Fold', 'field': 'fold', 'sorter': 'number',
     'formatter': 'fixed', 'formatterParams': {'precision': 2},
     'headerTooltip': 'Ratio / Bg ratio.'},
    {'title': 'p-value', 'field': 'pvalue', 'sorter': 'number',
     'formatter': 'scientific', 'formatterParams': {'precision': 2},
     'headerTooltip': "One-sided Fisher's exact test."},
    {'title': 'Adj. p-value', 'field': 'padj', 'sorter': 'number',
     'formatter': 'scientific', 'formatterParams': {'precision': 2},
     'headerTooltip': 'Benjamini-Hochberg adjusted, within the category.'},
    {'title': 'Proteins', 'field': 'proteins', 'sorter': 'string', 'minWidth': 200, 'hozAlign': 'left',
     'headerTooltip': 'The selected proteins annotated with the term.'},
]


def _digest(*parts):
    return hashlib.md5(json.dumps(parts, sort_keys=True, default=str).encode()).hexdigest()[:12]


def _result(file_manager, selected, background):
    """The enrichment of a selection, computed once per session."""
    store = st.session_state.setdefault('flash_go_results', {})
    digest = _digest(sorted(selected), sorted(background))
    if digest not in store:
        with st.spinner(f'Looking up the GO terms of {len(set(selected) | set(background))} proteins…'):
            store[digest] = go_terms.enrichment(file_manager, selected, background)
        while len(store) > _RESULTS_KEPT:
            store.pop(next(iter(store)))
    return digest, store[digest]


def _bar_chart(table, category):
    top = table.sort_values('p_adj').head(TOP_TERMS).iloc[::-1]
    score = -np.log10(top['p_adj'].clip(lower=1e-300))
    fig = go.Figure(go.Bar(
        x=score, y=top['GO_Term'], orientation='h',
        marker=dict(color=top['FoldEnrichment'], colorscale='Viridis', showscale=True,
                    colorbar=dict(title='Fold')),
        customdata=np.stack([top['GeneRatio'], top['p_adj']], axis=-1),
        hovertemplate='%{y}<br>Selected: %{customdata[0]}<br>adj. p: %{customdata[1]:.2e}<extra></extra>',
    ))
    fig.add_vline(x=-np.log10(0.05), line_dash='dash', line_color='#E74C3C',
                  annotation_text='adj. p = 0.05', annotation_position='top')
    fig.update_layout(
        title=f'{go_terms.CATEGORIES[category]}: top {len(top)} terms',
        xaxis_title='-log10(adjusted p-value)', yaxis_title=None,
        height=max(280, 60 + 24 * len(top)), margin=dict(l=10, r=10, t=50, b=40),
    )
    fig.update_yaxes(automargin=True)
    return fig


def _term_table(file_manager, scope, digest, category, frame):
    cache_path = Path(file_manager.cache_path, 'insight', 'go')
    cache_id = f'go-{scope}-{digest}-{category}'
    try:
        return Table(cache_id=cache_id, cache_path=str(cache_path))
    except CacheMissError:
        pass
    # Tables of earlier selections of this view are not shown again
    if cache_path.exists():
        for old in cache_path.glob(f'go-{scope}-*'):
            if not old.name.startswith(f'go-{scope}-{digest}-'):
                shutil.rmtree(old, ignore_errors=True)
    return Table(
        cache_id=cache_id, cache_path=str(cache_path), data=frame.lazy(),
        column_definitions=TERM_COLUMNS, index_field='rank', title=f'GO terms ({len(frame)})',
        initial_sort=[{'column': 'padj', 'dir': 'asc'}], layout='fitColumns',
        pagination=False, go_to_fields=['term'], float_precision='float64',
    )


def _clear(state_manager, identifier):
    state_manager.clear_selection(identifier)


def _render(file_manager, scope, state_manager, identifier, selected, background, key, picked_line):
    if state_manager.get_selection(identifier):
        head, clear = st.columns([5, 1], vertical_alignment='center')
        head.markdown(picked_line)
        clear.button(
            'Clear selection', key=f'{key}_clear', on_click=_clear, args=(state_manager, identifier),
            use_container_width=True,
        )
    else:
        st.info(picked_line)
    if not selected:
        return
    digest, (result, error, labels) = _result(file_manager, selected, background)
    if error:
        st.warning(error + ' Proteins already looked up are used; try again when the network is back.')
    if not result['fg_count']:
        st.info('None of the selected proteins has a GO annotation at MyGene.info.')
        return
    category = st.segmented_control(
        'GO category', list(go_terms.CATEGORIES), default='BP', key=f'{key}_category',
        format_func=go_terms.CATEGORIES.get, label_visibility='collapsed',
    ) or 'BP'
    table = result['categories'][category]
    st.caption(
        f"{result['fg_count']} of {result['fg_total']} selected and {result['bg_count']} of "
        f"{result['bg_total']} background proteins have GO terms (MyGene.info, by UniProt accession). "
        "One-sided Fisher's exact test, Benjamini-Hochberg adjusted per category."
    )
    if table is None or not len(table):
        st.info(f'The selected proteins have no {go_terms.CATEGORIES[category].lower()} terms.')
        return
    show_fig(_bar_chart(table, category), f'GO-{category}-{scope}')
    _term_table(file_manager, scope, digest, category, go_terms.term_frame(table, labels))(
        key=f'{key}_terms_{category}', state_manager=state_manager, height=360,
    )


def render_run(file_manager, dataset_id, state_manager, key):
    """GO enrichment of the proteins picked in a run, against all its
    identified proteins."""
    selected = [a for a in (state_manager.get_selection(go_terms.RUN_SET) or []) if a]
    background = go_terms.run_background(file_manager, dataset_id)
    if not background:
        st.info('This run has no identified proteins.')
        return
    if selected:
        line = f'**{len(selected)} of {len(background)} identified proteins selected.**'
    else:
        line = RUN_HINT
    _render(
        file_manager, _digest('run', dataset_id), state_manager, go_terms.RUN_SET,
        selected, background, key, line,
    )


def render_comparison(session, state_manager, key):
    """GO enrichment of the entries picked in a comparison, against every
    tested entry."""
    noun = {'feature': 'features', 'proteoform': 'proteoforms', 'protein': 'proteins'}[session.level]
    labels = state_manager.get_selection(go_terms.COMPARISON_SET) or []
    selected, background, picked, without = go_terms.comparison_selection(session, labels)
    if not background:
        st.info(f'None of the tested {noun} belongs to an identified protein, so GO terms cannot be tested.')
        return
    if picked:
        line = f'**{picked} {noun} selected: {len(selected)} of {len(background)} tested proteins.**'
        if without:
            line += f' {without} of them have no identified protein.'
    else:
        line = COMPARISON_HINT.format(nouns=noun)
    _render(
        session.file_manager, _digest('comparison', session.dataset_id, session.level), state_manager,
        go_terms.COMPARISON_SET, selected, background, key, line,
    )
