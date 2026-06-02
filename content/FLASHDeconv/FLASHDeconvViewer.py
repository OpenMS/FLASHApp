import streamlit as st

from pathlib import Path

from src.common.common import page_setup, save_params
from src.workflow.FileManager import FileManager
# Legacy bespoke-grid render path (kept importable until OI integration is verified).
from src.render.render import render_grid
# The OpenMS-Insight viewer (Stage B) is imported lazily inside render_panel (see
# below) so an import failure (e.g. a missing openms-insight install) falls back
# to the legacy grid instead of breaking the whole page.


def _use_oi_viewer():
    return st.session_state.get("settings", {}).get(
        "use_openms_insight_viewer", True
    )


def render_panel(experiment_id, layout_info_per_exp, file_manager, identifier,
                 grid_key, panel_index):
    """Render one experiment panel via the configured viewer.

    Routes to the new OpenMS-Insight viewer when enabled, else the legacy grid.
    The OI viewer is imported lazily and guarded so an import failure falls back
    to the legacy grid rather than breaking the page.
    """
    if _use_oi_viewer():
        try:
            from content.FLASHDeconv.FLASHDeconvViewerOI import (
                render_experiment_panel,
            )
        except Exception as exc:  # noqa: BLE001 - OI viewer unavailable
            st.warning(
                f"OpenMS-Insight viewer unavailable ({exc}); using legacy grid."
            )
        else:
            render_experiment_panel(
                experiment_id, layout_info_per_exp, file_manager, panel_index
            )
            return
    render_grid(
        experiment_id, layout_info_per_exp, file_manager,
        'flashdeconv', identifier, grid_key
    )

DEFAULT_LAYOUT = [['ms1_deconv_heat_map'], ['scan_table', 'mass_table'],
                  ['anno_spectrum', 'deconv_spectrum'], ['3D_SN_plot']]

def select_experiment():
    # Map display name back to experiment ID
    st.session_state.selected_experiment0 = display_name_to_id[st.session_state.selected_experiment_dropdown]
    if len(layout) > 1:
        for exp_index in range(1, len(layout)):
            if st.session_state[f'selected_experiment_dropdown_{exp_index}'] is None:
                continue
            st.session_state[f"selected_experiment{exp_index}"] = display_name_to_id[st.session_state[f'selected_experiment_dropdown_{exp_index}']]

def validate_selected_index(file_manager, selected_experiment):
    results = file_manager.get_results_list(['deconv_dfs', 'anno_dfs'])
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
            key="selected_experiment_dropdown",
            index=validate_selected_index(file_manager, 'selected_experiment0'),
            on_change=select_experiment
        )
        if 'selected_experiment0' in st.session_state:
            render_panel(
                st.session_state.selected_experiment0, layout[0], file_manager,
                "selected_experiment0", 'flash_viewer_grid_0', panel_index=0
            )
    with c2:
        st.selectbox(
            "choose experiment", display_names,
            key=f'selected_experiment_dropdown_1',
            index=validate_selected_index(file_manager, 'selected_experiment1'),
            on_change=select_experiment
        )
        if f"selected_experiment1" in st.session_state:
            with st.spinner('Loading component...'):
                render_panel(
                     st.session_state["selected_experiment1"], layout[1],
                     file_manager, 'selected_experiment1',
                     'flash_viewer_grid_1', panel_index=1
                )

else:
    ### for only single experiment on one view
    st.selectbox(
        "choose experiment", display_names,
        key="selected_experiment_dropdown",
        index=validate_selected_index(file_manager, 'selected_experiment0'),
        on_change=select_experiment
    )


    if 'selected_experiment0' in st.session_state:
        render_panel(
            st.session_state.selected_experiment0, layout[0], file_manager,
            'selected_experiment0', 'flash_viewer_grid', panel_index=0
        )

    ### for multiple experiments on one view
    if len(layout) > 1:

        for exp_index, exp_layout in enumerate(layout):
            if exp_index == 0: continue  # skip the first experiment

            st.divider()  # horizontal line

            st.selectbox(
                "choose experiment", display_names,
                key=f'selected_experiment_dropdown_{exp_index}',
                index=validate_selected_index(file_manager, f'selected_experiment{exp_index}'),
                on_change=select_experiment
            )
            # if #experiment input files are less than #layouts, all the pre-selection will be the first experiment
            if f"selected_experiment{exp_index}" in st.session_state:
                render_panel(
                     st.session_state["selected_experiment%d" % exp_index],
                     layout[exp_index], file_manager,
                     "selected_experiment%d" % exp_index,
                     'flash_viewer_grid_%d' % exp_index, panel_index=exp_index
                )

save_params(params)
