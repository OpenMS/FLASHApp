"""Comparisons of runs in the FLASHDeconv and FLASHTnT Viewers.

A workflow run that quantified two or more files stores a comparison
(src/quant/session.py) next to its runs. The Viewer lists it in the experiment
selector as "Comparison: ..."; selecting it shows a grid of Quantification
components (src/render/quant_cells.py) with its own layout, edited in the same
window as the run layout. **Edit comparison** opens a window with the
conditions of the runs and the analysis settings.
"""

import streamlit as st
from openms_insight import StateManager, defer_reruns

from src.quant.session import LEVEL_NAMES, QuantSession, comparisons
from src.render import quant_cells
from src.render.layout_editor import get_saved_layout, layout_editor
from src.render.layout_editor.catalog import QUANT_COMPONENTS, QUANT_DEFAULT_LAYOUT
from src.render.layout_editor.dialogs import _HELP, _own_file_manager

LAYOUT_DATASET = {'flashdeconv': 'quant_layout', 'flashtnt': 'flashtnt_quant_layout'}


def comparison_names(file_manager, tool):
    """Display name -> dataset id of the comparisons `tool` stored, newest first."""
    names = {}
    for dataset_id in comparisons(file_manager, tool):
        manifest = file_manager.get_results(dataset_id, ['quant_manifest'])['quant_manifest']
        stamp = dataset_id.rsplit('_', 1)[-1]
        when = f'{stamp[:4]}-{stamp[4:6]}-{stamp[6:8]} {stamp[9:11]}:{stamp[11:13]}' if len(stamp) >= 13 else stamp
        names[f"Comparison: {len(manifest['runs'])} runs ({when})"] = dataset_id
    return names


def layout(file_manager, tool):
    saved = get_saved_layout(file_manager, LAYOUT_DATASET[tool])
    return saved[0][0] if saved else QUANT_DEFAULT_LAYOUT


@st.dialog('Edit layout', width='large')
def _layout_dialog(file_manager, tool):
    file_manager = _own_file_manager(file_manager)
    st.caption(
        _HELP.replace(', add experiments to compare several files', '')
        + ' This layout is used for every comparison of runs.'
    )
    if layout_editor(file_manager, LAYOUT_DATASET[tool], QUANT_COMPONENTS, QUANT_DEFAULT_LAYOUT,
                     key=f'{tool}_quant_layout_editor', max_experiments=1):
        st.rerun()


@st.dialog('Edit comparison', width='large')
def _comparison_dialog(file_manager, dataset_id):
    session = QuantSession(_own_file_manager(file_manager), dataset_id)
    st.caption(
        'Which runs belong to which condition, which level to compare and how. '
        'The defaults suit most experiments; **Apply** recomputes the comparison.'
    )
    if quant_cells.settings_form(session):
        st.rerun()


def buttons(file_manager, dataset_id, tool):
    """Edit comparison and Edit layout, right of the experiment selector."""
    if st.button('Edit comparison', icon=':material/tune:', key=f'{tool}_edit_comparison', width='content',
                 help='Conditions of the runs, level and analysis settings'):
        _comparison_dialog(file_manager, dataset_id)
    if st.button('Edit layout', icon=':material/dashboard_customize:', key=f'{tool}_edit_quant_layout',
                 width='content', help='Choose and arrange the components shown for a comparison'):
        _layout_dialog(file_manager, tool)


def render(file_manager, dataset_id, tool):
    """The comparison grid of a Viewer."""
    session = QuantSession(file_manager, dataset_id)
    level = session.level
    if len(session.levels) == 1 and tool == 'flashtnt':
        st.caption('No PrSM identified a compared feature, so only features can be compared.')
    st.caption(f'Comparing **{LEVEL_NAMES[level].lower()}** of {len(session.runs)} runs. '
               'Switch the level and set the conditions with **Edit comparison**.')
    # Selections are per comparison and level: an index means another entry on another level
    state_manager = StateManager(session_key=f'flash_quant_state_{tool}_{dataset_id}_{level}')
    with defer_reruns():
        for row_index, row in enumerate(layout(file_manager, tool)):
            columns = st.columns(len(row))
            for col_index, (column, comp_name) in enumerate(zip(columns, row)):
                with column:
                    quant_cells.render_cell(
                        comp_name, session, state_manager, f'quant_{tool}_{level}_{row_index}_{col_index}',
                    )
