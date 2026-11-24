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
from pathlib import Path

import streamlit as st

from src.render import go_terms, teaching
from src.viz import go_panel

# The state identifier the chart writes its GO category to; the table follows it
CATEGORY_IDENTIFIER = go_panel.CATEGORY_IDENTIFIER
TERM_COUNTS = go_panel.TERM_COUNTS
TOP_TERMS = go_panel.TOP_TERMS
# Results kept per session (one per selection)
_RESULTS_KEPT = 8

RUN_HINT = (
    'Select proteins to test which GO terms are over-represented in them: use box or lasso select in the mode bar of the '
    'proteoform heatmap (Shift adds), or turn on the check boxes of the protein, proteoform or '
    'PrSM table with its check mark button and tick rows.'
)
COMPARISON_HINT = (
    'Select {nouns} to test which GO terms are over-represented in them: **Select** in the toolbar of the volcano plot '
    '(All changed, Higher, Lower), box or lasso select on the volcano plot or the clustered heatmap (Shift adds), '
    'or tick rows of the comparison table.'
)
_OFFLINE = ' Proteins already looked up are used; try again when the network is back.'


def _digest(*parts):
    return hashlib.md5(json.dumps(parts, sort_keys=True, default=str).encode()).hexdigest()[:12]


def _stored(digest, message, compute):
    store = st.session_state.setdefault('flash_go_results', {})
    if digest not in store:
        if message:
            with st.spinner(message):
                store[digest] = compute()
        else:
            store[digest] = compute()
        while len(store) > _RESULTS_KEPT:
            store.pop(next(iter(store)))
    return store[digest]


def _message(file_manager, pool):
    # The workflow looks the terms up after its run; only results of an
    # earlier build look them up here
    if go_terms.index_ready(file_manager, pool):
        return None
    return f'Looking up the GO terms of {len(set(pool))} proteins…'


def _render(file_manager, state_manager, identifier, selected, background, key, line, hint, universe=None):
    pool = list(universe if universe is not None else background)
    # Filled by whichever of the two the panel shows, before it draws its table
    labels = {}

    def enrichment():
        result, error, names = _stored(
            _digest(sorted(selected), sorted(background)), _message(file_manager, pool + list(selected)),
            lambda: go_terms.enrichment(file_manager, selected, background, universe=universe),
        )
        labels.update(names)
        return (result['categories'] if result['fg_count'] else None), (error + _OFFLINE if error else None)

    def counts():
        annotated, tables, error, names = _stored(
            _digest('all', sorted(background), sorted(pool)), _message(file_manager, pool),
            lambda: go_terms.run_terms(file_manager, background, universe=universe),
        )
        labels.update(names)
        return annotated, tables, (error + _OFFLINE if error else None)

    go_panel.show_go_panel(
        key, Path(file_manager.cache_path, 'insight'), state_manager, identifier, bool(selected), line,
        enrichment, counts, labels, hint=hint, help=teaching.help(),
    )


def render_run(file_manager, dataset_id, state_manager, key):
    """GO enrichment of the proteins picked in a run, against all its
    identified proteins."""
    selected = [a for a in (state_manager.get_selection(go_terms.RUN_SET) or []) if a]
    background = go_terms.stored_run_background(file_manager, dataset_id)
    if not background:
        st.info('This run has no identified proteins.')
        return
    line = f'**{len(selected)} of {len(background)} identified proteins selected.**'
    _render(file_manager, state_manager, go_terms.RUN_SET, selected, background, key, line, RUN_HINT)


def render_comparison(session, state_manager, key):
    """GO enrichment of the entries picked in a comparison, against every
    tested entry."""
    noun = {'feature': 'features', 'proteoform': 'proteoforms', 'protein': 'proteins'}[session.level]
    labels = state_manager.get_selection(go_terms.COMPARISON_SET) or []
    selected, background, picked, without = go_terms.comparison_selection(session, labels)
    if not background:
        st.info(f'None of the tested {noun} belongs to an identified protein, so GO terms cannot be tested.')
        return
    line = f'**{picked} {noun} selected: {len(selected)} of {len(background)} tested proteins.**'
    if without:
        line += f' {without} of them have no identified protein.'
    hint = COMPARISON_HINT.format(nouns=noun)
    if picked and not selected:
        # Picked entries without an identified protein: nothing to test
        hint = line
    _render(
        session.file_manager, state_manager, go_terms.COMPARISON_SET, selected, background, key, line, hint,
        universe=go_terms.comparison_universe(session),
    )
