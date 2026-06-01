import os

import streamlit as st

from pathlib import Path

from src.common.common import page_setup, save_params
from src.workflow.FileManager import FileManager
from src.render.render import render_grid

# Migration flag (shared with FLASHDeconv): when truthy, render each experiment
# panel with the OpenMS-Insight engine instead of the legacy flash_viewer_grid.
USE_OPENMS_INSIGHT = os.environ.get("FLASHAPP_USE_OPENMS_INSIGHT", "").lower() in (
    "1", "true", "yes",
)


DEFAULT_LAYOUT = [
    ['protein_table'],
    ['sequence_view'],
    ['tag_table'],
    ['combined_spectrum']
]


def render_panel_tnt(dataset_id, layout_rows, file_manager, exp_key, grid_key=None):
    """Render one FLASHTnT experiment panel with the selected engine."""
    if USE_OPENMS_INSIGHT:
        from src.render_oi import render_experiment_tnt

        render_experiment_tnt(
            dataset_id, layout_rows, file_manager,
            panel_key=(grid_key or exp_key),
        )
    elif grid_key is not None:
        render_grid(dataset_id, layout_rows, file_manager, 'flashtnt', exp_key, grid_key)
    else:
        render_grid(dataset_id, layout_rows, file_manager, 'flashtnt', exp_key)


def select_experiment():
    # Map display name back to experiment ID
    st.session_state.selected_experiment0_tagger = display_name_to_id[st.session_state.selected_experiment_dropdown_tagger]
    if len(layout) > 1:
        for exp_index in range(1, len(layout)):
            if st.session_state[f'selected_experiment_dropdown_{exp_index}_tagger'] is None:
                continue
            st.session_state[f"selected_experiment{exp_index}_tagger"] = display_name_to_id[st.session_state[f'selected_experiment_dropdown_{exp_index}_tagger']]

def validate_selected_index(file_manager, selected_experiment):
    results = file_manager.get_results_list(
        ['deconv_dfs', 'anno_dfs', 'tag_dfs', 'protein_dfs']
    )
    if selected_experiment in st.session_state:
        if st.session_state[selected_experiment] in results:
            # Map experiment ID to display name for the dropdown index
            exp_id = st.session_state[selected_experiment]
            display_name = file_manager.get_display_name(exp_id)
            return display_name_to_index[display_name]
        else:
            del st.session_state[selected_experiment]
    return None

# page initialization
params = page_setup("TaggerViewer")

# Get available results
file_manager = FileManager(
    st.session_state["workspace"],
    Path(st.session_state['workspace'], 'cache')
)
results = file_manager.get_results_list(
    ['protein_dfs']
)

if file_manager.result_exists('flashtnt_layout', 'layout'):
    layout = file_manager.get_results('flashtnt_layout', 'layout')['layout']
    side_by_side = layout['side_by_side']
    layout = layout['layout']

else:
    layout = [DEFAULT_LAYOUT]
    side_by_side = False

### if no input file is given, show blank page
if len(results) == 0:
    st.error('No results to show yet. Please run a workflow first!')
    st.stop()

# Create display names and mappings
display_names = [file_manager.get_display_name(exp_id) for exp_id in results]
display_name_to_id = {file_manager.get_display_name(exp_id): exp_id for exp_id in results}
display_name_to_index = {n : i for i, n in enumerate(display_names)}
# Keep backward compatibility mapping for experiment IDs
name_to_index = {n : i for i, n in enumerate(results)}

if len(layout) == 2 and side_by_side:
    c1, c2 = st.columns(2)
    with c1:
        st.selectbox(
            "choose experiment", display_names,
            key="selected_experiment_dropdown_tagger",
            index=validate_selected_index(file_manager, 'selected_experiment0_tagger'),
            on_change=select_experiment
        )
        if 'selected_experiment0_tagger' in st.session_state:
            render_panel_tnt(st.session_state.selected_experiment0_tagger, layout[0], file_manager, 'selected_experiment0_tagger')
    with c2:
        st.selectbox(
            "choose experiment", display_names,
            key=f'selected_experiment_dropdown_1_tagger',
            index=validate_selected_index(file_manager, 'selected_experiment1_tagger'),
            on_change=select_experiment
        )
        if f"selected_experiment1_tagger" in st.session_state:
            render_panel_tnt(st.session_state.selected_experiment1_tagger, layout[1], file_manager, 'selected_experiment1_tagger', 'flash_viewer_grid_1')


else:
    ### for only single experiment on one view
    st.selectbox(
        "choose experiment", display_names,
        key="selected_experiment_dropdown_tagger",
        index=validate_selected_index(file_manager, 'selected_experiment0_tagger'),
        on_change=select_experiment
    )

    if 'selected_experiment0_tagger' in st.session_state:
        render_panel_tnt(st.session_state.selected_experiment0_tagger, layout[0], file_manager, 'selected_experiment0_tagger')

    ### for multiple experiments on one view
    if len(layout) > 1:

        for exp_index, exp_layout in enumerate(layout):
            if exp_index == 0: continue  # skip the first experiment

            st.divider() # horizontal line

            st.selectbox(
                "choose experiment", display_names,
                key=f'selected_experiment_dropdown_{exp_index}_tagger',
                index=validate_selected_index(file_manager, f'selected_experiment{exp_index}_tagger'),
                on_change=select_experiment
            )

            # if #experiment input files are less than #layouts, all the pre-selection will be the first experiment
            if f"selected_experiment{exp_index}_tagger" in st.session_state:
                render_panel_tnt(st.session_state["selected_experiment%d_tagger" % exp_index], layout[exp_index], file_manager, f"selected_experiment{exp_index}_tagger", 'flash_viewer_grid_%d' % exp_index)

save_params(params)