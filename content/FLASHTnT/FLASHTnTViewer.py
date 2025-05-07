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
    if len(layout) > 1:
        for exp_index in range(1, len(layout)):
            if st.session_state[f'selected_experiment_dropdown_{exp_index}_tagger'] is None:
                continue
            st.session_state[f"selected_experiment{exp_index}_tagger"] = st.session_state[f'selected_experiment_dropdown_{exp_index}_tagger']


# page initialization
params = page_setup("TaggerViewer")

# Get available results
file_manager = FileManager(
    st.session_state["workspace"],
    Path(st.session_state['workspace'], 'flashtnt', 'cache')
)
results = file_manager.get_results_list(
    ['deconv_dfs', 'anno_dfs', 'tag_dfs', 'protein_dfs']
)

if file_manager.result_exists('layout', 'layout'):
    layout = file_manager.get_results('layout', 'layout')['layout']
    side_by_side = layout['side_by_side']
    layout = layout['layout']
    
else:
    layout = [DEFAULT_LAYOUT]
    side_by_side = False

### if no input file is given, show blank page
if len(results) == 0:
    st.error('No results to show yet. Please run a workflow first!')
    st.stop()

# Map names to index
name_to_index = {n : i for i, n in enumerate(results)}

if len(layout) == 2 and side_by_side:
    c1, c2 = st.columns(2)
    with c1:
        st.selectbox(
            "choose experiment", results, 
            key="selected_experiment_dropdown_tagger", 
            index=name_to_index[st.session_state.selected_experiment0_tagger] if 'selected_experiment0_tagger' in st.session_state else None,
            on_change=select_experiment
        )
        if 'selected_experiment0_tagger' in st.session_state:
            render_grid(st.session_state.selected_experiment0_tagger, layout[0], file_manager, 'flashtnt', 'selected_experiment0_tagger')
    with c2:
        st.selectbox(
            "choose experiment", results, 
            key=f'selected_experiment_dropdown_1_tagger',
            index = name_to_index[st.session_state[f'selected_experiment1_tagger']] if f'selected_experiment1_tagger' in st.session_state else None,
            on_change=select_experiment
        )
        if f"selected_experiment1_tagger" in st.session_state:
            render_grid(st.session_state.selected_experiment1_tagger, layout[1], file_manager, 'flashtnt', 'selected_experiment1_tagger', 'flash_viewer_grid_1')


else:
    ### for only single experiment on one view
    st.selectbox(
        "choose experiment", results, 
        key="selected_experiment_dropdown_tagger", 
        index=name_to_index[st.session_state.selected_experiment0_tagger] if 'selected_experiment0_tagger' in st.session_state else None,
        on_change=select_experiment
    )

    if 'selected_experiment0_tagger' in st.session_state:
        render_grid(st.session_state.selected_experiment0_tagger, layout[0], file_manager, 'flashtnt', 'selected_experiment0_tagger')

    ### for multiple experiments on one view
    if len(layout) > 1:

        for exp_index, exp_layout in enumerate(layout):
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
                render_grid(st.session_state["selected_experiment%d_tagger" % exp_index], layout[exp_index], file_manager, 'flashtnt', f"selected_experiment{exp_index}_tagger", 'flash_viewer_grid_%d' % exp_index)

save_params(params)