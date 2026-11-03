"""Comparisons of runs in the Viewer.

A workflow run that quantified two or more files stores a comparison
(src/quant/session.py) next to its runs. The Viewer lists it in the experiment
selector as "Comparison: ..."; selecting it shows a grid of Quantification
components (src/render/quant_cells.py) with its own layout, edited in the same
window as the run layout. **Edit comparison** opens a window with the
conditions of the runs and the analysis settings.
"""

import streamlit as st

from src.quant import isobaric
from src.quant.session import LEVEL_NAMES, QuantSession, comparisons
from src.render.layout_editor import get_saved_layout, layout_editor
from src.render.layout_editor.catalog import QUANT_COMPONENTS, QUANT_DEFAULT_LAYOUT
from src.render.layout_editor.dialogs import _HELP, _own_file_manager

# One layout serves every comparison. Before the workflows were merged, a
# comparison of identified runs had its own (migrate_layout).
LAYOUT_DATASET = 'quant_layout'
LEGACY_LAYOUT_DATASET = 'flashtnt_quant_layout'


def migrate_layout(file_manager):
    """Keep the layout saved for comparisons of identified runs when no other
    comparison layout was saved; drop it either way."""
    if not file_manager.result_exists(LEGACY_LAYOUT_DATASET, 'layout'):
        return
    if not file_manager.result_exists(LAYOUT_DATASET, 'layout'):
        legacy = file_manager.get_results(LEGACY_LAYOUT_DATASET, 'layout')['layout']
        file_manager.store_data(LAYOUT_DATASET, 'layout', legacy)
    file_manager.remove_results(LEGACY_LAYOUT_DATASET)


def comparison_names(file_manager, tool=None):
    """Display name -> dataset id of the comparisons (of `tool`, default all), newest first."""
    names = {}
    for dataset_id in comparisons(file_manager, tool):
        manifest = file_manager.get_results(dataset_id, ['quant_manifest'])['quant_manifest']
        stamp = dataset_id.rsplit('_', 1)[-1]
        when = f'{stamp[:4]}-{stamp[4:6]}-{stamp[6:8]} {stamp[9:11]}:{stamp[11:13]}' if len(stamp) >= 13 else stamp
        if manifest.get('kind') == 'isobaric':
            what = f"{isobaric.label(manifest.get('plex', ''))}, {len(manifest['runs'])} samples"
        else:
            what = f"{len(manifest['runs'])} runs"
        names[f"Comparison: {what} ({when})"] = dataset_id
    return names


def comparison_tool(file_manager, dataset_id):
    """'flashtnt' for a comparison of identified runs, else 'flashdeconv'."""
    manifest = file_manager.get_results(dataset_id, ['quant_manifest'])['quant_manifest']
    return manifest.get('tool') or 'flashdeconv'


def layout(file_manager):
    saved = get_saved_layout(file_manager, LAYOUT_DATASET)
    return saved[0][0] if saved else QUANT_DEFAULT_LAYOUT


@st.dialog('Edit layout', width='large')
def _layout_dialog(file_manager):
    file_manager = _own_file_manager(file_manager)
    st.caption(
        _HELP.replace(', add experiments to compare several files', '')
        + ' This layout is used for every comparison of runs.'
    )
    if layout_editor(file_manager, LAYOUT_DATASET, QUANT_COMPONENTS, QUANT_DEFAULT_LAYOUT,
                     key='quant_layout_editor', max_experiments=1):
        st.rerun()


@st.dialog('Edit comparison', width='large')
def _comparison_dialog(file_manager, dataset_id):
    from src.render import quant_cells  # needs openms_insight
    session = QuantSession(_own_file_manager(file_manager), dataset_id)
    st.caption(
        'Which runs (or labelled samples) belong to which condition, which level to compare and how. '
        'The defaults suit most experiments; **Apply** recomputes the comparison.'
    )
    if quant_cells.settings_form(session):
        st.rerun()


def buttons(file_manager, dataset_id):
    """Edit comparison and Edit layout, right of the experiment selector."""
    if st.button('Edit comparison', icon=':material/tune:', key='edit_comparison', width='content',
                 help='Conditions of the runs, level and analysis settings'):
        _comparison_dialog(file_manager, dataset_id)
    if st.button('Edit layout', icon=':material/dashboard_customize:', key='edit_quant_layout',
                 width='content', help='Choose and arrange the components shown for a comparison'):
        _layout_dialog(file_manager)


def render(file_manager, dataset_id):
    """The comparison grid of the Viewer."""
    # openms_insight only where the comparison is drawn, so the layout helpers
    # above import without it (the unit tests do not install it)
    from openms_insight import StateManager, defer_reruns
    from src.render import quant_cells
    tool = comparison_tool(file_manager, dataset_id)
    session = QuantSession(file_manager, dataset_id)
    level = session.level
    if len(session.levels) == 1 and tool == 'flashtnt':
        st.caption('No PrSM identified a compared feature, so only features can be compared.')
    if session.isobaric:
        fractions = len(session.manifest.get('fractions', []))
        compared = (f"{len(session.runs)} {isobaric.label(session.manifest.get('plex', ''))} samples"
                    + (f' over {fractions} runs' if fractions > 1 else ''))
    else:
        compared = f'{len(session.runs)} runs'
    st.caption(f'Comparing **{LEVEL_NAMES[level].lower()}** of {compared}. '
               'Switch the level and set the conditions with **Edit comparison**.')
    # Selections are per comparison and level: an index means another entry on another level
    state_manager = StateManager(session_key=f'flash_quant_state_{tool}_{dataset_id}_{level}')
    with defer_reruns():
        for row_index, row in enumerate(layout(file_manager)):
            columns = st.columns(len(row))
            for col_index, (column, comp_name) in enumerate(zip(columns, row)):
                with column:
                    quant_cells.render_cell(
                        comp_name, session, state_manager, f'quant_{tool}_{level}_{row_index}_{col_index}',
                    )
