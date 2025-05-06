import streamlit as st

from pathlib import Path

from src.common.common import page_setup, save_params
from src.workflow.FileManager import FileManager
from src.render.render import render_grid


DEFAULT_LAYOUT = [
    ['protein_table'], 
    ['sequence_view'], 
    ['tag_table'],
    ['combined_spectrum']
]


def select_experiment():
    st.session_state.selected_experiment0_tagger = st.session_state.selected_experiment_dropdown_tagger
    if "saved_layout_setting_tagger" in st.session_state and len(st.session_state["saved_layout_setting_tagger"]) > 1:
        for exp_index in range(1, len(st.session_state["saved_layout_setting_tagger"])):
            if st.session_state[f'selected_experiment_dropdown_{exp_index}_tagger'] is None:
                continue
            st.session_state[f"selected_experiment{exp_index}_tagger"] = st.session_state[f'selected_experiment_dropdown_{exp_index}_tagger']


# page initialization
params = page_setup("TaggerViewer")
st.title('FLASHViewer')

# Get available results
file_manager = FileManager(
    st.session_state["workspace"],
    Path(st.session_state['workspace'], 'flashtnt', 'cache')
)
results = file_manager.get_results_list(
    ['deconv_dfs', 'anno_dfs', 'tag_dfs', 'protein_dfs']
)

### if no input file is given, show blank page
if len(results) == 0:
    st.error('No results to show yet. Please run a workflow first!')
    st.stop()

# Map names to index
name_to_index = {n : i for i, n in enumerate(results)}

### for only single experiment on one view
st.selectbox(
    "choose experiment", results, 
    key="selected_experiment_dropdown_tagger", 
    index=name_to_index[st.session_state.selected_experiment0_tagger] if 'selected_experiment0_tagger' in st.session_state else None,
    on_change=select_experiment
)

if 'selected_experiment0_tagger' in st.session_state:
    layout_info = DEFAULT_LAYOUT
    if "saved_layout_setting_tagger" in st.session_state:  # when layout manager was used
        layout_info = st.session_state["saved_layout_setting_tagger"][0]
    render_grid(st.session_state.selected_experiment0_tagger, layout_info, file_manager, 'flashtnt', 'selected_experiment0_tagger')


### for multiple experiments on one view
if "saved_layout_setting_tagger" in st.session_state and len(st.session_state["saved_layout_setting_tagger"]) > 1:

    for exp_index, exp_layout in enumerate(st.session_state["saved_layout_setting_tagger"]):
        if exp_index == 0: continue  # skip the first experiment

        st.divider() # horizontal line

        st.selectbox(
            "choose experiment", results, 
            key=f'selected_experiment_dropdown_{exp_index}_tagger',
            index = name_to_index[st.session_state[f'selected_experiment{exp_index}_tagger']] if f'selected_experiment{exp_index}_tagger' in st.session_state else None,
            on_change=select_experiment
        )

        # if #experiment input files are less than #layouts, all the pre-selection will be the first experiment
        if f"selected_experiment{exp_index}_tagger" in st.session_state:
            layout_info = st.session_state["saved_layout_setting_tagger"][exp_index]
            with st.spinner('Loading component...'):
                render_grid(st.session_state["selected_experiment%d_tagger" % exp_index], layout_info, file_manager, 'flashtnt', f"selected_experiment{exp_index}_tagger", 'flash_viewer_grid_%d' % exp_index)

save_params(params)