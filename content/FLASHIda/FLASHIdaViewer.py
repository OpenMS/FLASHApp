import streamlit as st

from pathlib import Path

from src.common.common import page_setup, save_params
from src.workflow.FileManager import FileManager
from src.render.render import render_grid

DEFAULT_LAYOUT = [['ms1_deconv_heat_map']]

def select_experiment():
    st.session_state.selected_experiment0_ida = st.session_state.selected_experiment_dropdown_ida
    print(st.session_state.selected_experiment0_ida)
    if len(layout) > 1:
        for exp_index in range(1, len(layout)):
            if st.session_state[f'selected_experiment_dropdown_{exp_index}_ida'] is None:
                continue
            st.session_state[f"selected_experiment{exp_index}_ida"] = st.session_state[f'selected_experiment_dropdown_{exp_index}_ida']

def validate_selected_index(file_manager, selected_experiment):
    results = file_manager.get_results_list(['simulation_dfs'])
    if selected_experiment in st.session_state:
        if st.session_state[selected_experiment] in results:
            return name_to_index[st.session_state[selected_experiment]]
        else:
            del st.session_state[selected_experiment]
    return None

# page initialization
params = page_setup()

# Get available results
file_manager = FileManager(
    st.session_state["workspace"],
    Path(st.session_state['workspace'], 'flashidasimulator', 'cache')
)

results = file_manager.get_results_list(['simulation_dfs'])

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
            key="selected_experiment_dropdown_ida", 
            index=validate_selected_index(file_manager, 'selected_experiment0_ida'),
            on_change=select_experiment
        )
        if 'selected_experiment0_ida' in st.session_state:
            render_grid(
                st.session_state.selected_experiment0_ida, layout[0], file_manager, 
                'flashidasimulator', "selected_experiment0_ida", 'flash_viewer_grid_0_ida'
            )
    with c2:
        st.selectbox(
            "choose experiment", results, 
            key=f'selected_experiment_dropdown_1_ida',
            index=validate_selected_index(file_manager, 'selected_experiment1_ida'),
            on_change=select_experiment
        )
        if f"selected_experiment1_ida" in st.session_state:
            with st.spinner('Loading component...'):
                render_grid(
                     st.session_state["selected_experiment1_ida"], layout[1], 
                     file_manager, 'flashidasimulator', 'selected_experiment1_ida', 
                     'flash_viewer_grid_1_ida'
                )

else:
    ### for only single experiment on one view
    st.selectbox(
        "choose experiment", results, 
        key="selected_experiment_dropdown_ida", 
        index=validate_selected_index(file_manager, 'selected_experiment0_ida'),
        on_change=select_experiment
    )


    if 'selected_experiment0_ida' in st.session_state:
        print('Lets go!')
        render_grid(
            st.session_state.selected_experiment0_ida, layout[0], file_manager, 
            'flashidasimulator', 'selected_experiment0_ida'
        )

    ### for multiple experiments on one view
    if len(layout) > 1:

        for exp_index, exp_layout in enumerate(layout):
            if exp_index == 0: continue  # skip the first experiment

            st.divider()  # horizontal line

            st.selectbox(
                "choose experiment", results, 
                key=f'selected_experiment_dropdown_{exp_index}_ida',
                index=validate_selected_index(file_manager, f'selected_experiment{exp_index}_ida'),
                on_change=select_experiment
            )
            # if #experiment input files are less than #layouts, all the pre-selection will be the first experiment
            if f"selected_experiment{exp_index}_ida" in st.session_state:
                render_grid(
                     st.session_state["selected_experiment%d_ida" % exp_index], 
                     layout[exp_index], file_manager, 'flashidasimulator', 
                     "selected_experiment%d_ida" % exp_index, 
                     'flash_viewer_grid_%d_ida' % exp_index
                )

save_params(params)
