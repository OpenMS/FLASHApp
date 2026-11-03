"""The Layout Manager as a window on the Viewer.

``layout_button()`` draws an "Edit layout" button (the Viewer puts it right of
the experiment selector); pressing it opens the drag-and-drop editor in a
dialog over the Viewer. Saving or resetting reruns the page, which closes the
dialog and shows the new layout.
"""

import streamlit as st

from src.workflow.FileManager import FileManager
from src.render import topdown_layout
from src.render.layout_editor import layout_editor
from src.render.layout_editor.catalog import TNT_DEFAULT_LAYOUTS

EDITOR_KEY = 'topdown_layout_editor'
_START_KEY = 'topdown_layout_start'
DECONVOLUTION = 'Deconvolution'

_HELP = (
    "Drag components into the grid (up to 3 per row), add experiments to "
    "compare several files, then press **Save**."
)


def _own_file_manager(file_manager):
    """A FileManager on the same cache with its own SQLite connection.

    A dialog reruns as a fragment, possibly on another script thread, and an
    SQLite connection may only be used by the thread that opened it.
    """
    return FileManager(file_manager.workflow_dir, file_manager.cache_path)


def _start_from(name):
    """Put a default layout into the editor as an unsaved draft."""
    st.session_state[_START_KEY] = name
    revision = f'{EDITOR_KEY}_revision'
    st.session_state[revision] = st.session_state.get(revision, 0) + 1


def _defaults(file_manager, labelled=False):
    defaults = {DECONVOLUTION: topdown_layout.deconvolution_default(file_manager)}
    defaults.update(TNT_DEFAULT_LAYOUTS)
    if labelled:
        defaults = {
            name: topdown_layout.with_reporters(rows, topdown_layout.REPORTER_AFTER)
            for name, rows in defaults.items()
        }
    return defaults


@st.dialog("Edit layout", width="large")
def _dialog(file_manager):
    file_manager = _own_file_manager(file_manager)
    components = topdown_layout.components(file_manager)
    note = '' if topdown_layout.has_sequence(file_manager) else (
        " Submit a sequence on the **Sequence Input** page to show a sequence "
        "of your own with a run deconvolved only."
    )
    st.caption(
        _HELP + " Until a layout is saved, each experiment starts from the default "
        "of its run: **Deconvolution** for a run deconvolved only, **Native** or "
        "**Lysate** for an identified run, as it was searched, each with the reporter "
        "ions of a run with isobaric labels." + note + " Start from any of them:"
    )
    labelled = st.toggle(
        'With reporter ions (TMT, iTRAQ)', key='topdown_layout_start_labelled',
        value=bool(file_manager.get_results_list(['isobaric_dfs'])),
    )
    defaults = _defaults(file_manager, labelled)
    with st.container(horizontal=True):
        for name in defaults:
            st.button(f'{name} default', key=f'topdown_layout_start_{name}',
                      on_click=_start_from, args=(name,))

    # A save or reset ends the started draft
    submitted = st.session_state.get(EDITOR_KEY)
    if submitted and submitted.get('nonce') != st.session_state.get(f'{EDITOR_KEY}_handled'):
        st.session_state.pop(_START_KEY, None)
    start = st.session_state.get(_START_KEY)

    if layout_editor(file_manager, topdown_layout.LAYOUT_DATASET, components,
                     defaults[DECONVOLUTION], key=EDITOR_KEY,
                     start_layout=[defaults[start]] if start in defaults else None):
        st.rerun()


def layout_button(file_manager):
    """Draw the "Edit layout" button of the Viewer."""
    if st.button('Edit layout', icon=':material/dashboard_customize:',
                 key='topdown_edit_layout', width='content',
                 help='Choose and arrange the components shown on this page'):
        _dialog(file_manager)
