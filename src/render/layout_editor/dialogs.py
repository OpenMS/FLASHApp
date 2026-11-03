"""The Layout Manager as a window on the Viewer pages.

``layout_button()`` draws an "Edit layout" button (the Viewers put it right of
the experiment selector); pressing it opens the
drag-and-drop editor in a dialog over the Viewer. Saving or resetting reruns
the page, which closes the dialog and shows the new layout.
"""

import streamlit as st

from src.workflow.FileManager import FileManager
from src.render.layout_editor import layout_editor
from src.render.layout_editor.catalog import (
    DECONV_COMPONENTS, DECONV_SEQUENCE_COMPONENTS, DECONV_DEFAULT_LAYOUT,
    TNT_COMPONENTS, TNT_DEFAULT_LAYOUT, TNT_DEFAULT_LAYOUTS,
)

DECONV_EDITOR_KEY = 'deconv_layout_editor'
TNT_EDITOR_KEY = 'tnt_layout_editor'
_TNT_START_KEY = 'tnt_layout_start'

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


def deconv_components(file_manager):
    """Components and default layout offered for FLASHDeconv results."""
    components = DECONV_COMPONENTS
    default_layout = DECONV_DEFAULT_LAYOUT
    # "Sequence view" and "Internal fragment map" need a submitted sequence
    if file_manager.result_exists('sequence', 'sequence'):
        components = components + DECONV_SEQUENCE_COMPONENTS
        default_layout = default_layout + [['sequence_view']]
    return components, default_layout


@st.dialog("Edit layout", width="large")
def _deconv_dialog(file_manager):
    file_manager = _own_file_manager(file_manager)
    components, default_layout = deconv_components(file_manager)
    note = '' if len(components) > len(DECONV_COMPONENTS) else (
        " Submit a sequence on the **Sequence Input** page to unlock the "
        "sequence components."
    )
    st.caption(_HELP + note)
    if layout_editor(file_manager, 'layout', components, default_layout, key=DECONV_EDITOR_KEY):
        st.rerun()


def _tnt_start_from(sample_type):
    """Put a default layout into the editor as an unsaved draft."""
    st.session_state[_TNT_START_KEY] = sample_type
    revision = f'{TNT_EDITOR_KEY}_revision'
    st.session_state[revision] = st.session_state.get(revision, 0) + 1


@st.dialog("Edit layout", width="large")
def _tnt_dialog(file_manager):
    file_manager = _own_file_manager(file_manager)
    st.caption(
        _HELP + " Without a saved layout, the Viewer shows the **Native** default "
        "for native or targeted runs (few spectra) and the **Lysate** default for "
        "LC-MS/MS runs of a lysate, chosen per experiment. Start from either:"
    )
    with st.container(horizontal=True):
        for sample_type in TNT_DEFAULT_LAYOUTS:
            st.button(f'{sample_type} default', key=f'tnt_layout_start_{sample_type}',
                      on_click=_tnt_start_from, args=(sample_type,))

    # A save or reset ends the started draft
    submitted = st.session_state.get(TNT_EDITOR_KEY)
    if submitted and submitted.get('nonce') != st.session_state.get(f'{TNT_EDITOR_KEY}_handled'):
        st.session_state.pop(_TNT_START_KEY, None)
    start = st.session_state.get(_TNT_START_KEY)

    if layout_editor(file_manager, 'flashtnt_layout', TNT_COMPONENTS, TNT_DEFAULT_LAYOUT,
                     key=TNT_EDITOR_KEY,
                     start_layout=[TNT_DEFAULT_LAYOUTS[start]] if start else None):
        st.rerun()


def layout_button(file_manager, tool):
    """Draw the "Edit layout" button of a Viewer; ``tool`` is 'flashdeconv' or 'flashtnt'."""
    dialog = _deconv_dialog if tool == 'flashdeconv' else _tnt_dialog
    if st.button('Edit layout', icon=':material/dashboard_customize:',
                 key=f'{tool}_edit_layout', width='content',
                 help='Choose and arrange the components shown on this page'):
        dialog(file_manager)
