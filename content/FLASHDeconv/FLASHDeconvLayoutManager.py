from pathlib import Path

import streamlit as st

from src.common.common import page_setup, save_params
from src.workflow.FileManager import FileManager
from src.render.layout_editor import layout_editor
from src.render.layout_editor.catalog import (
    DECONV_COMPONENTS, DECONV_SEQUENCE_COMPONENTS, DECONV_DEFAULT_LAYOUT
)

# page_setup() first: it initialises st.session_state['workspace'], which the
# cache access below needs when this page is opened by direct link.
params = page_setup()

# Setup cache access
file_manager = FileManager(
    st.session_state["workspace"],
    Path(st.session_state['workspace'], 'cache')
)

components = DECONV_COMPONENTS
default_layout = DECONV_DEFAULT_LAYOUT
# "Sequence view" and "Internal fragment map" need a submitted sequence
if file_manager.result_exists('sequence', 'sequence'):
    components = components + DECONV_SEQUENCE_COMPONENTS
    default_layout = default_layout + [['sequence_view']]

st.title("Layout Manager")
st.caption(
    "Arrange the components shown on the **👀 Viewer** page. Drag components "
    "into the grid (up to 3 per row), add experiments to compare several files, "
    "then press **Save**. Submit a sequence on the Sequence Input page to unlock "
    "the sequence components."
)

layout_editor(file_manager, 'layout', components, default_layout, key='deconv_layout_editor')

save_params(params)
