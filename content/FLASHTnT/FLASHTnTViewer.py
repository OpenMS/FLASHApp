import streamlit as st

from pathlib import Path

from src.common.common import page_setup, save_params
from src.workflow.FileManager import FileManager
from src import sample_type
from src.render.render import render_grid
from src.render.layout_editor.catalog import TNT_DEFAULT_LAYOUTS
from src.render.layout_editor.dialogs import layout_button
from src.render import quant_view


def experiment_layout(saved_layout, dataset_id):
    """The saved layout of an experiment, or the default for its sample type."""
    if saved_layout is not None:
        return saved_layout
    kind = sample_type.of_result(file_manager, dataset_id) or sample_type.NATIVE
    st.caption(
        f"{kind} default layout: this experiment looks like a "
        f"{'native or targeted run' if kind == sample_type.NATIVE else 'lysate LC-MS/MS run'}. "
        "Press **Edit layout** to arrange your own."
    )
    return TNT_DEFAULT_LAYOUTS[kind]


def select_experiment():
    # Map display name back to experiment ID
    # A cleared selectbox reports None; show no experiment instead of failing.
    if st.session_state.selected_experiment_dropdown_tagger is None:
        st.session_state.pop('selected_experiment0_tagger', None)
    else:
        st.session_state.selected_experiment0_tagger = display_name_to_id[st.session_state.selected_experiment_dropdown_tagger]
    if len(layout) > 1:
        for exp_index in range(1, len(layout)):
            # Not drawn while a comparison is shown
            if st.session_state.get(f'selected_experiment_dropdown_{exp_index}_tagger') is None:
                continue
            st.session_state[f"selected_experiment{exp_index}_tagger"] = display_name_to_id[st.session_state[f'selected_experiment_dropdown_{exp_index}_tagger']]

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
    return st.session_state.get('selected_experiment0_tagger') in comparison_ids

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
    # No saved layout: each experiment gets the default for its sample type
    layout = [None]
    side_by_side = False

### if no input file is given, show blank page
comparison_names = quant_view.comparison_names(file_manager, 'flashtnt')
comparison_ids = set(comparison_names.values())

if len(results) == 0 and not comparison_names:
    layout_button(file_manager, 'flashtnt')
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
            key="selected_experiment_dropdown_tagger",
            index=validate_selected_index(file_manager, 'selected_experiment0_tagger'),
            on_change=select_experiment
        )
        if 'selected_experiment0_tagger' in st.session_state:
            render_grid(st.session_state.selected_experiment0_tagger, experiment_layout(layout[0], st.session_state.selected_experiment0_tagger), file_manager, 'flashtnt', 'selected_experiment0_tagger')
    with c2:
        with st.container(horizontal=True, vertical_alignment='bottom'):
            st.selectbox(
                "choose experiment", run_names,
                key=f'selected_experiment_dropdown_1_tagger',
                index=validate_selected_index(file_manager, 'selected_experiment1_tagger', run_names),
                on_change=select_experiment
            )
            layout_button(file_manager, 'flashtnt')
        if f"selected_experiment1_tagger" in st.session_state:
            render_grid(st.session_state.selected_experiment1_tagger, experiment_layout(layout[1], st.session_state.selected_experiment1_tagger), file_manager, 'flashtnt', 'selected_experiment1_tagger', 'flash_viewer_grid_1')


else:
    ### for only single experiment on one view
    with st.container(horizontal=True, vertical_alignment='bottom'):
        st.selectbox(
            "choose experiment", display_names,
            key="selected_experiment_dropdown_tagger",
            index=validate_selected_index(file_manager, 'selected_experiment0_tagger'),
            on_change=select_experiment
        )
        if comparison_selected():
            quant_view.buttons(file_manager, st.session_state.selected_experiment0_tagger, 'flashtnt')
        else:
            layout_button(file_manager, 'flashtnt')

    if comparison_selected():
        # A comparison spans its runs, so it fills the page on its own
        quant_view.render(file_manager, st.session_state.selected_experiment0_tagger, 'flashtnt')

    elif 'selected_experiment0_tagger' in st.session_state:
        render_grid(st.session_state.selected_experiment0_tagger, experiment_layout(layout[0], st.session_state.selected_experiment0_tagger), file_manager, 'flashtnt', 'selected_experiment0_tagger')

    ### for multiple experiments on one view
    if len(layout) > 1 and not comparison_selected():

        for exp_index, exp_layout in enumerate(layout):
            if exp_index == 0: continue  # skip the first experiment

            st.divider() # horizontal line

            st.selectbox(
                "choose experiment", run_names,
                key=f'selected_experiment_dropdown_{exp_index}_tagger',
                index=validate_selected_index(file_manager, f'selected_experiment{exp_index}_tagger', run_names),
                on_change=select_experiment
            )

            # if #experiment input files are less than #layouts, all the pre-selection will be the first experiment
            if f"selected_experiment{exp_index}_tagger" in st.session_state:
                render_grid(st.session_state["selected_experiment%d_tagger" % exp_index], experiment_layout(exp_layout, st.session_state["selected_experiment%d_tagger" % exp_index]), file_manager, 'flashtnt', f"selected_experiment{exp_index}_tagger", 'flash_viewer_grid_%d' % exp_index)

save_params(params)