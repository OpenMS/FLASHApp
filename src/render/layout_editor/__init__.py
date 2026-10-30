"""Drag-and-drop Layout Manager.

The editor runs entirely in the browser (``frontend/index.html``, plain
JavaScript with no build step), so arranging components never reruns the
script. Python only hears from it when the user presses Save or Reset.
"""

from pathlib import Path

import streamlit as st
import streamlit.components.v1 as st_components

from src.render.layout_editor.catalog import MAX_COLUMNS, MAX_EXPERIMENTS, validate_layout

_editor = st_components.declare_component(
    'flash_layout_editor', path=str(Path(__file__).parent / 'frontend')
)


def get_saved_layout(file_manager, dataset_id):
    """Return ``(layout, side_by_side)`` stored under ``dataset_id``, or ``None``."""
    if not file_manager.result_exists(dataset_id, 'layout'):
        return None
    saved = file_manager.get_results(dataset_id, 'layout')['layout']
    return saved['layout'], saved['side_by_side']


def layout_editor(file_manager, dataset_id, components, default_layout, key):
    """Render the editor and persist the layout it submits.

    The saved format is unchanged from the previous Layout Manager:
    ``{'layout': [[[name, ...], ...], ...], 'side_by_side': bool}`` under
    ``file_manager.store_data(dataset_id, 'layout', ...)``.
    """
    revision_key = '%s_revision' % key
    handled_key = '%s_handled' % key
    st.session_state.setdefault(revision_key, 0)

    # A submitted value is in session state at the start of the rerun it
    # triggers, so persist it before rendering and the editor gets fresh args.
    submitted = st.session_state.get(key)
    error = ''
    if submitted and submitted.get('nonce') != st.session_state.get(handled_key):
        st.session_state[handled_key] = submitted.get('nonce')
        if submitted.get('action') == 'save':
            error = validate_layout(submitted.get('layout'), components)
            if not error:
                layout = submitted['layout']
                file_manager.store_data(dataset_id, 'layout', {
                    'layout': layout,
                    'side_by_side': bool(submitted.get('side_by_side')) and len(layout) == 2,
                })
        elif submitted.get('action') == 'reset':
            if file_manager.result_exists(dataset_id, 'layout'):
                file_manager.remove_results(dataset_id)
            st.session_state[revision_key] += 1

    saved = get_saved_layout(file_manager, dataset_id)
    _editor(
        components=components,
        default_layout=default_layout,
        saved_layout=saved[0] if saved else None,
        side_by_side=saved[1] if saved else False,
        max_experiments=MAX_EXPERIMENTS,
        max_columns=MAX_COLUMNS,
        revision=st.session_state[revision_key],
        server_error=error,
        key=key,
        default=None,
    )
