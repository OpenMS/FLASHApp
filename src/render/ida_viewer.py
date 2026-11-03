"""Experiment picker and fixed grid of the two FLASHIda viewers."""
from pathlib import Path

import streamlit as st

from src.render.render import render_grid
from src.workflow.FileManager import FileManager


def show(workflow_dir, result, layout, tool, empty_message):
    file_manager = FileManager(
        st.session_state['workspace'],
        Path(st.session_state['workspace'], workflow_dir, 'cache'),
    )
    results = file_manager.get_results_list([result])
    if not results:
        st.info(empty_message)
        return
    names = {file_manager.get_display_name(r): r for r in results}
    key = f'{tool}-experiment'
    if st.session_state.get(key) not in names:
        st.session_state[key] = list(names)[-1]
    choice = st.selectbox('choose experiment', list(names), key=key)
    render_grid(names[choice], layout, file_manager, tool, f'{tool}_experiment0', f'{tool}_grid')
