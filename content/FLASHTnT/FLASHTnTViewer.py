import streamlit as st

from pathlib import Path

from src.common.common import page_setup, save_params
from src.workflow.FileManager import FileManager
# Legacy bespoke-grid render path (kept importable until OI integration is verified).
from src.render.render import render_grid
# The OpenMS-Insight viewer (Stage C) is imported lazily inside render_panel (see
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
            from content.FLASHTnT.FLASHTnTViewerOI import (
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
        'flashtnt', identifier, grid_key
    )


DEFAULT_LAYOUT = [
    ['protein_table'],
    ['sequence_view'],
    ['tag_table'],
    ['combined_spectrum']
]


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
            render_panel(
                st.session_state.selected_experiment0_tagger, layout[0],
                file_manager, 'selected_experiment0_tagger',
                'flash_viewer_grid_0', panel_index=0
            )
    with c2:
        st.selectbox(
            "choose experiment", display_names,
            key=f'selected_experiment_dropdown_1_tagger',
            index=validate_selected_index(file_manager, 'selected_experiment1_tagger'),
            on_change=select_experiment
        )
        if f"selected_experiment1_tagger" in st.session_state:
            with st.spinner('Loading component...'):
                render_panel(
                    st.session_state.selected_experiment1_tagger, layout[1],
                    file_manager, 'selected_experiment1_tagger',
                    'flash_viewer_grid_1', panel_index=1
                )


else:
    ### for only single experiment on one view
    st.selectbox(
        "choose experiment", display_names,
        key="selected_experiment_dropdown_tagger",
        index=validate_selected_index(file_manager, 'selected_experiment0_tagger'),
        on_change=select_experiment
    )

    if 'selected_experiment0_tagger' in st.session_state:
        render_panel(
            st.session_state.selected_experiment0_tagger, layout[0],
            file_manager, 'selected_experiment0_tagger',
            'flash_viewer_grid', panel_index=0
        )

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
                render_panel(
                    st.session_state["selected_experiment%d_tagger" % exp_index],
                    layout[exp_index], file_manager,
                    f"selected_experiment{exp_index}_tagger",
                    'flash_viewer_grid_%d' % exp_index, panel_index=exp_index
                )

save_params(params)