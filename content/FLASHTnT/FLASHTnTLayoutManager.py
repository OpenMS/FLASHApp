from pathlib import Path

import streamlit as st

from src.common.common import page_setup, save_params
from src.workflow.FileManager import FileManager
from src.render.layout_editor import layout_editor
from src.render.layout_editor.catalog import TNT_COMPONENTS, TNT_DEFAULT_LAYOUT

# page_setup() first: it initialises st.session_state['workspace'], which the
# cache access below needs when this page is opened by direct link.
params = page_setup()

# Setup cache access
file_manager = FileManager(
    st.session_state["workspace"],
    Path(st.session_state['workspace'], 'cache')
)

st.title("Layout Manager")
st.caption(
    "Arrange the components shown on the **👀 Viewer** page. Drag components "
    "into the grid (up to 3 per row), add experiments to compare several files, "
    "then press **Save**."
)

layout_editor(file_manager, 'flashtnt_layout', TNT_COMPONENTS, TNT_DEFAULT_LAYOUT, key='tnt_layout_editor')

save_params(params)
