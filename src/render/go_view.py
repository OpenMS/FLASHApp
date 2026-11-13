"""The GO enrichment cell of the Viewer, for a run and for a comparison.

The cell reads the proteins (run) or entries (comparison) picked in the other
cells of the grid from their shared StateManager (see src/render/go_terms.py)
and shows the Gene Ontology terms over-represented in them: a bar chart of the
most significant terms and a table of all of them. While nothing is picked,
it shows the GO annotations of the whole run (or of every tested entry of a
comparison): how many proteins carry each term.
"""

import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import polars as pl
import streamlit as st

from openms_insight import Table
from openms_insight.core.cache import CacheMissError

from src.render import go_terms
from src.render.series import show_series

# Terms the bar chart can show
TERM_COUNTS = [10, 15, 20, 30, 50]
TOP_TERMS = 15
# Results kept per session (one per selection)
_RESULTS_KEPT = 8

RUN_HINT = (
    'Select proteins to test which GO terms are over-represented in them: use box or lasso select in the mode bar of the '
    'proteoform heatmap (Shift adds), or turn on the check boxes of the protein, proteoform or '
    'PrSM table with its check mark button and tick rows.'
)
COMPARISON_HINT = (
    'Select {nouns} to test which GO terms are over-represented in them: use **All changed**, **Higher** or **Lower** above '
    'the volcano plot, or box or lasso select in the mode bar of the volcano plot or the clustered '
    'heatmap (Shift adds); or turn on the check boxes of the comparison table with its check mark '
    'button and tick rows.'
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


COUNT_COLUMNS = [
    {'title': '#', 'field': 'rank', 'sorter': 'number', 'width': 50},
    {'title': 'GO term', 'field': 'term', 'sorter': 'string', 'minWidth': 220, 'hozAlign': 'left'},
    {'title': 'Proteins', 'field': 'count', 'sorter': 'number',
     'headerTooltip': 'Proteins annotated with the term.'},
    {'title': '% of annotated', 'field': 'share', 'sorter': 'number',
     'formatter': 'fixed', 'formatterParams': {'precision': 1},
     'headerTooltip': 'Proteins with the term, in percent of the proteins with any GO term.'},
    {'title': 'Annotated proteins', 'field': 'proteins', 'sorter': 'string', 'minWidth': 200,
     'hozAlign': 'left', 'headerTooltip': 'The proteins annotated with the term.'},
]


def _digest(*parts):
    return hashlib.md5(json.dumps(parts, sort_keys=True, default=str).encode()).hexdigest()[:12]


def _result(file_manager, selected, background, universe=None):
    """The enrichment of a selection, computed once per session. The GO
    terms of the proteins are looked up once, by the workflow or by the first
    selection of results from an earlier build; a selection only runs the
    tests."""
    store = st.session_state.setdefault('flash_go_results', {})
    digest = _digest(sorted(selected), sorted(background))
    if digest not in store:
        pool = list(universe if universe is not None else background) + list(selected)
        if go_terms.index_ready(file_manager, pool):
            message = 'Testing GO terms…'
        else:
            message = f'Looking up the GO terms of {len(set(pool))} proteins…'
        with st.spinner(message):
            store[digest] = go_terms.enrichment(file_manager, selected, background, universe=universe)
        while len(store) > _RESULTS_KEPT:
            store.pop(next(iter(store)))
    return digest, store[digest]


def _bar_chart(table, category, count):
    """The `count` most significant terms, as the arguments of a horizontal bar chart."""
    top = table.sort_values('p_adj').head(count)
    padj = top['p_adj'].astype(float).clip(lower=1e-300)
    data = pl.DataFrame({
        'term': top['GO_Term'].astype(str).tolist(),
        'score': (-np.log10(padj)).round(3).tolist(),
        'hover': [
            f'{t}: {c} selected ({r}), adj. p {p:.2e}, fold {f:.2f}'
            for t, c, r, p, f in zip(top['GO_Term'], top['Count'], top['GeneRatio'], padj, top['FoldEnrichment'])
        ],
    })
    return dict(
        data=data, x_column='term', y_column='score', mode='bars', orientation='h', show_legend=False,
        hover_column='hover', title=f'{go_terms.CATEGORIES[category]}: the most significant terms',
        x_label='GO term', y_label='-log10(adjusted p-value)',
    )


def _term_table(file_manager, scope, digest, category, frame, columns=TERM_COLUMNS,
                sort=({'column': 'padj', 'dir': 'asc'},)):
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
        column_definitions=list(columns), index_field='rank', title=f'GO terms ({len(frame)})',
        initial_sort=list(sort), layout='fitColumns',
        pagination=False, go_to_fields=['term'], float_precision='float64',
    )


def _clear(state_manager, identifier):
    # Set to None, not removed: the components keep a selection Python no
    # longer sends and would send it back
    state_manager.set_selection(identifier, None)


def _pick_terms(key):
    pick, size = st.columns([4, 1], vertical_alignment='bottom')
    with pick:
        category = st.segmented_control(
            'GO category', list(go_terms.CATEGORIES), default='BP', key=f'{key}_category',
            format_func=go_terms.CATEGORIES.get, label_visibility='collapsed',
        ) or 'BP'
    with size:
        count = st.selectbox(
            'Terms in the chart', TERM_COUNTS, index=TERM_COUNTS.index(TOP_TERMS), key=f'{key}_count',
            help='How many terms the bar chart shows. The table lists all of them.',
        )
    return category, count


def _count_chart(table, category, count, total):
    """The `count` most frequent terms, as the arguments of a horizontal bar chart."""
    top = table.head(count)
    counts = top['count'].to_list()
    data = pl.DataFrame({
        'term': top['term'].to_list(),
        'proteins': counts,
        'hover': [
            f'{t}: {c} proteins ({100.0 * c / total if total else 0.0:.1f}% of annotated)'
            for t, c in zip(top['term'].to_list(), counts)
        ],
    })
    return dict(
        data=data, x_column='term', y_column='proteins', mode='bars', orientation='h', show_legend=False,
        hover_column='hover', title=f'{go_terms.CATEGORIES[category]}: the most frequent terms',
        x_label='GO term', y_label='Proteins',
    )


def _chart(file_manager, spec, name, key):
    data = spec.pop('data')
    show_series(data, name, file_manager.cache_path, key, height=max(320, 100 + 30 * data.height), **spec)


def _render_all(file_manager, scope, state_manager, background, key, what, universe=None):
    """The GO annotations of every protein of `background`, while nothing is
    selected: how many proteins carry each term."""
    pool = list(universe if universe is not None else background)
    store = st.session_state.setdefault('flash_go_results', {})
    digest = _digest('all', sorted(background), sorted(pool))
    if digest not in store:
        if go_terms.index_ready(file_manager, pool):
            message = 'Counting GO terms…'
        else:
            message = f'Looking up the GO terms of {len(set(pool))} proteins…'
        with st.spinner(message):
            store[digest] = go_terms.run_terms(file_manager, background, universe=universe)
        while len(store) > _RESULTS_KEPT:
            store.pop(next(iter(store)))
    annotated, tables, error, labels = store[digest]
    if error:
        st.warning(error + ' Proteins already looked up are used; try again when the network is back.')
    if not annotated:
        st.info(f'No protein among {what} has a GO annotation at MyGene.info.')
        return
    category, count = _pick_terms(key)
    table = tables[category]
    st.caption(
        f'GO annotations of {what}: {annotated} of {len(set(background))} proteins have GO terms '
        '(MyGene.info, by UniProt accession).'
    )
    if table is None or not table.height:
        st.info(f'No protein has {go_terms.CATEGORIES[category].lower()} terms.')
        return
    _chart(file_manager, _count_chart(table, category, count, annotated), f'go-{scope}-all-{category}', f'{key}_all_chart')
    _term_table(
        file_manager, scope, f'all{digest}', category, go_terms.count_frame(table, labels, annotated),
        columns=COUNT_COLUMNS, sort=({'column': 'count', 'dir': 'desc'},),
    )(key=f'{key}_all_{category}', state_manager=state_manager, height=360)


def _render(file_manager, scope, state_manager, identifier, selected, background, key, picked_line,
            universe=None, what=None):
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
        if what:
            _render_all(file_manager, scope, state_manager, background, key, what, universe)
        return
    digest, (result, error, labels) = _result(file_manager, selected, background, universe)
    if error:
        st.warning(error + ' Proteins already looked up are used; try again when the network is back.')
    if not result['fg_count']:
        st.info('None of the selected proteins has a GO annotation at MyGene.info.')
        return
    category, count = _pick_terms(key)
    table = result['categories'][category]
    st.caption(
        f"{result['fg_count']} of {result['fg_total']} selected and {result['bg_count']} of "
        f"{result['bg_total']} background proteins have GO terms (MyGene.info, by UniProt accession). "
        "One-sided Fisher's exact test, Benjamini-Hochberg adjusted per category."
    )
    if table is None or not len(table):
        st.info(f'The selected proteins have no {go_terms.CATEGORIES[category].lower()} terms.')
        return
    _chart(file_manager, _bar_chart(table, category, count), f'go-{scope}-{category}', f'{key}_chart')
    _term_table(file_manager, scope, digest, category, go_terms.term_frame(table, labels))(
        key=f'{key}_terms_{category}', state_manager=state_manager, height=360,
    )


def render_run(file_manager, dataset_id, state_manager, key):
    """GO enrichment of the proteins picked in a run, against all its
    identified proteins."""
    selected = [a for a in (state_manager.get_selection(go_terms.RUN_SET) or []) if a]
    background = go_terms.stored_run_background(file_manager, dataset_id)
    if not background:
        st.info('This run has no identified proteins.')
        return
    if selected:
        line = f'**{len(selected)} of {len(background)} identified proteins selected.**'
    else:
        line = RUN_HINT
    _render(
        file_manager, _digest('run', dataset_id), state_manager, go_terms.RUN_SET,
        selected, background, key, line, what=f'all {len(background)} identified proteins of the run',
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
        universe=go_terms.comparison_universe(session),
        what=f'the proteins of all tested {noun}',
    )
