import streamlit as st

from pathlib import Path

from src.common.common import page_setup, save_params
from src.workflow.FileManager import FileManager
from src.render.render import render_grid
from src.render.layout_editor.dialogs import layout_button
from src.render import quant_view

DEFAULT_LAYOUT = [['ms1_deconv_heat_map'], ['scan_table', 'mass_table'],
                  ['anno_spectrum', 'deconv_spectrum'], ['3D_SN_plot']]

def select_experiment():
    # Map display name back to experiment ID
    # A cleared selectbox reports None; show no experiment instead of failing.
    if st.session_state.selected_experiment_dropdown is None:
        st.session_state.pop('selected_experiment0', None)
    else:
        st.session_state.selected_experiment0 = display_name_to_id[st.session_state.selected_experiment_dropdown]
    if len(layout) > 1:
        for exp_index in range(1, len(layout)):
            # Not drawn while a comparison is shown
            if st.session_state.get(f'selected_experiment_dropdown_{exp_index}') is None:
                continue
            st.session_state[f"selected_experiment{exp_index}"] = display_name_to_id[st.session_state[f'selected_experiment_dropdown_{exp_index}']]

def validate_selected_index(file_manager, selected_experiment, names=None):
    """Index of the stored selection in the selector (`names`, default every
    entry), or None after dropping a selection that no longer exists."""
    names = display_names if names is None else names
    if selected_experiment in st.session_state:
        display_name = display_id_to_name.get(st.session_state[selected_experiment])
        if display_name in names:
            return names.index(display_name)
        del st.session_state[selected_experiment]
    return None


def comparison_selected():
    """True while the first selector shows a comparison of runs."""
    return st.session_state.get('selected_experiment0') in comparison_ids

# page initialization
params = page_setup()

# Get available results
file_manager = FileManager(
    st.session_state["workspace"],
    Path(st.session_state['workspace'], 'cache')
)

def get_sequence():
    # Check if layout has been set
    if not file_manager.result_exists('sequence', 'sequence'):
        return None
    # fetch layout from cache
    sequence = file_manager.get_results('sequence', 'sequence')['sequence']

    return sequence['input_sequence'], sequence['fixed_mod_cysteine'], sequence['fixed_mod_methionine'] 

if get_sequence() is not None:
    DEFAULT_LAYOUT = DEFAULT_LAYOUT + [['sequence_view']]

results = file_manager.get_results_list(['threedim_SN_plot'])

if file_manager.result_exists('layout', 'layout'):
    layout = file_manager.get_results('layout', 'layout')['layout']
    side_by_side = layout['side_by_side']
    layout = layout['layout']
    
else:
    layout = [DEFAULT_LAYOUT]
    side_by_side = False

### if no input file is given, show blank page
comparison_names = quant_view.comparison_names(file_manager, 'flashdeconv')
comparison_ids = set(comparison_names.values())

if len(results) == 0 and not comparison_names:
    layout_button(file_manager, 'flashdeconv')
    st.error('No results to show yet. Please run a workflow first!')
    st.stop()

# Create display names and mappings: the runs, then the comparisons of runs
run_names = [file_manager.get_display_name(exp_id) for exp_id in results]
display_names = run_names + list(comparison_names)
display_name_to_id = {file_manager.get_display_name(exp_id): exp_id for exp_id in results}
display_name_to_id.update(comparison_names)
display_id_to_name = {exp_id: name for name, exp_id in display_name_to_id.items()}

if len(layout) == 2 and side_by_side and not comparison_selected():
    c1, c2 = st.columns(2)
    with c1:
        st.selectbox(
            "choose experiment", display_names,
            key="selected_experiment_dropdown",
            index=validate_selected_index(file_manager, 'selected_experiment0'),
            on_change=select_experiment
        )
        if 'selected_experiment0' in st.session_state:
            render_grid(
                st.session_state.selected_experiment0, layout[0], file_manager, 
                'flashdeconv', "selected_experiment0", 'flash_viewer_grid_0'
            )
    with c2:
        with st.container(horizontal=True, vertical_alignment='bottom'):
            st.selectbox(
                "choose experiment", run_names,
                key=f'selected_experiment_dropdown_1',
                index=validate_selected_index(file_manager, 'selected_experiment1', run_names),
                on_change=select_experiment
            )
            layout_button(file_manager, 'flashdeconv')
        if f"selected_experiment1" in st.session_state:
            with st.spinner('Loading component...'):
                render_grid(
                     st.session_state["selected_experiment1"], layout[1], 
                     file_manager, 'flashdeconv', 'selected_experiment1', 
                     'flash_viewer_grid_1'
                )

else:
    ### for only single experiment on one view
    with st.container(horizontal=True, vertical_alignment='bottom'):
        st.selectbox(
            "choose experiment", display_names,
            key="selected_experiment_dropdown",
            index=validate_selected_index(file_manager, 'selected_experiment0'),
            on_change=select_experiment
        )
        if comparison_selected():
            quant_view.buttons(file_manager, st.session_state.selected_experiment0, 'flashdeconv')
        else:
            layout_button(file_manager, 'flashdeconv')

    if comparison_selected():
        # A comparison spans its runs, so it fills the page on its own
        quant_view.render(file_manager, st.session_state.selected_experiment0, 'flashdeconv')

    elif 'selected_experiment0' in st.session_state:
        render_grid(
            st.session_state.selected_experiment0, layout[0], file_manager, 
            'flashdeconv', 'selected_experiment0'
        )

    ### for multiple experiments on one view
    if len(layout) > 1 and not comparison_selected():

        for exp_index, exp_layout in enumerate(layout):
            if exp_index == 0: continue  # skip the first experiment

            st.divider()  # horizontal line

            st.selectbox(
                "choose experiment", run_names,
                key=f'selected_experiment_dropdown_{exp_index}',
                index=validate_selected_index(file_manager, f'selected_experiment{exp_index}', run_names),
                on_change=select_experiment
            )
            # if #experiment input files are less than #layouts, all the pre-selection will be the first experiment
            if f"selected_experiment{exp_index}" in st.session_state:
                render_grid(
                     st.session_state["selected_experiment%d" % exp_index], 
                     layout[exp_index], file_manager, 'flashdeconv', 
                     "selected_experiment%d" % exp_index, 
                     'flash_viewer_grid_%d' % exp_index
                )

save_params(params)
