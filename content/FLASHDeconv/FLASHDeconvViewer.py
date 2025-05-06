import streamlit as st

from pathlib import Path

from src.common.common import page_setup, save_params
from src.workflow.FileManager import FileManager
from src.render.render import render_grid

DEFAULT_LAYOUT = [['ms1_deconv_heat_map'], ['scan_table', 'mass_table'],
                  ['anno_spectrum', 'deconv_spectrum'], ['3D_SN_plot']]
if 'input_sequence' in st.session_state and st.session_state.input_sequence:
    DEFAULT_LAYOUT = DEFAULT_LAYOUT + [['sequence_view']]

def select_experiment():
    st.session_state.selected_experiment0 = st.session_state.selected_experiment_dropdown
    if "saved_layout_setting" in st.session_state and len(st.session_state["saved_layout_setting"]) > 1:
        for exp_index in range(1, len(st.session_state["saved_layout_setting"])):
            if st.session_state[f'selected_experiment_dropdown_{exp_index}'] is None:
                continue
            st.session_state[f"selected_experiment{exp_index}"] = st.session_state[f'selected_experiment_dropdown_{exp_index}']

# page initialization
params = page_setup()
st.title("FLASHViewer")

# Get available results
file_manager = FileManager(
    st.session_state["workspace"],
    Path(st.session_state['workspace'], 'flashdeconv', 'cache')
)
results = file_manager.get_results_list(['deconv_dfs', 'anno_dfs'])

### if no input file is given, show blank page
if len(results) == 0:
    st.error('No results to show yet. Please run a workflow first!')
    st.stop()

# Map names to index
name_to_index = {n : i for i, n in enumerate(results)}

### for only single experiment on one view
st.selectbox(
    "choose experiment", results, 
    key="selected_experiment_dropdown", 
    index=name_to_index[st.session_state.selected_experiment0] if 'selected_experiment0' in st.session_state else None,
    on_change=select_experiment
)

if 'selected_experiment0' in st.session_state:
    layout_info = DEFAULT_LAYOUT
    if "saved_layout_setting" in st.session_state:  # when layout manager was used
        layout_info = st.session_state["saved_layout_setting"][0]
    render_grid(st.session_state.selected_experiment0, layout_info, file_manager, 'flashdeconv', 'selected_experiment0')


### for multiple experiments on one view
if "saved_layout_setting" in st.session_state and len(st.session_state["saved_layout_setting"]) > 1:

    for exp_index, exp_layout in enumerate(st.session_state["saved_layout_setting"]):
        if exp_index == 0: continue  # skip the first experiment

        st.divider()  # horizontal line

        st.selectbox(
            "choose experiment", results, 
            key=f'selected_experiment_dropdown_{exp_index}',
            index = name_to_index[st.session_state[f'selected_experiment{exp_index}']] if f'selected_experiment{exp_index}' in st.session_state else None,
            on_change=select_experiment
        )
        # if #experiment input files are less than #layouts, all the pre-selection will be the first experiment
        if f"selected_experiment{exp_index}" in st.session_state:
            layout_info = st.session_state["saved_layout_setting"][exp_index]
            render_grid(st.session_state["selected_experiment%d" % exp_index], layout_info, file_manager, 'flashdeconv', "selected_experiment%d" % exp_index, 'flash_viewer_grid_%d' % exp_index,)

save_params(params)
