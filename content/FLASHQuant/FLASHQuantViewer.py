import os

import streamlit as st

from pathlib import Path

from src.workflow.FileManager import FileManager
from src.common.common import page_setup, save_params
# from src.render.components import flash_viewer_grid_component, FlashViewerComponent, FLASHQuant
from src.render.render import render_grid

# Migration flag (shared across workflows): render with OpenMS-Insight engine.
USE_OPENMS_INSIGHT = os.environ.get("FLASHAPP_USE_OPENMS_INSIGHT", "").lower() in (
    "1", "true", "yes",
)

# page initialization
params = page_setup()


# Get available results
file_manager = FileManager(
    st.session_state["workspace"],
    Path(st.session_state['workspace'], 'flashquant', 'cache')
)
results = file_manager.get_results_list(
    ['quant_dfs']
)

### if no input file is given, show blank page
if len(results) == 0:
    st.error('No results to show yet. Please run a workflow first!')
    st.stop()

# Map names to index
name_to_index = {n : i for i, n in enumerate(results)}


# for only single experiment on one view
st.selectbox("choose experiment", results, key="selected_experiment0_quant")
selected_exp0 = st.session_state.selected_experiment0_quant

if USE_OPENMS_INSIGHT:
    from src.render_oi import render_experiment_quant

    render_experiment_quant(
        st.session_state.selected_experiment0_quant, file_manager,
        panel_key='selected_experiment0_quant',
    )
else:
    render_grid(
        st.session_state.selected_experiment0_quant, [['quant_visualization']],
        file_manager, 'flashquant', 'selected_experiment0_quant'
    )

# # Get data
# quant_df = file_manager.get_results(selected_exp0, 'quant_dfs')['quant_dfs']

# component = [[FlashViewerComponent(FLASHQuant())]]
# flash_viewer_grid_component(components=component, data={'quant_data': quant_df, 'dataset': selected_exp0}, component_key='flash_viewer_grid')

save_params(params)
