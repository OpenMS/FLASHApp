"""The Viewer of the merged workflow: every run, deconvolved and perhaps
identified, and every comparison of runs.

An experiment shows the identification components when its run was
identified and the deconvolution components otherwise; an identified run can
be switched to its deconvolution. Each kind keeps its own layout (Layout
Manager under ``flashtnt_layout`` and ``layout``), so the layouts saved before
the workflows were merged stay as they were.
"""

import streamlit as st

from pathlib import Path

from src.common.common import page_setup, save_params
from src.workflow.FileManager import FileManager
from src import sample_type
from src.render.render import render_grid
from src.render.layout_editor import get_saved_layout
from src.render.layout_editor.catalog import TNT_DEFAULT_LAYOUTS
from src.render.layout_editor.dialogs import deconv_components, layout_button
from src.render import quant_view

IDENTIFICATION, DECONVOLUTION = 'Identification', 'Deconvolution'
TOOL = {IDENTIFICATION: 'flashtnt', DECONVOLUTION: 'flashdeconv'}
LAYOUT_DATASET = {'flashtnt': 'flashtnt_layout', 'flashdeconv': 'layout'}


def selection_key(index):
    return f'selected_experiment{index}'


def view_of(index):
    """'flashtnt' or 'flashdeconv': what experiment `index` shows."""
    dataset_id = st.session_state.get(selection_key(index))
    if dataset_id not in identified:
        return 'flashdeconv'
    return TOOL[st.session_state.get(f'experiment_view_{index}', IDENTIFICATION)]


def saved_layout(tool):
    """(experiments, side_by_side) saved for `tool`; experiments are None
    while nothing is saved (the default is chosen per experiment)."""
    saved = get_saved_layout(file_manager, LAYOUT_DATASET[tool])
    return saved if saved else ([None], False)


def experiment_layout(tool, index, dataset_id):
    """The rows experiment `index` shows: its saved layout, else the default."""
    experiments, _ = saved_layout(tool)
    rows = experiments[index] if index < len(experiments) else experiments[0]
    if rows is not None:
        return rows
    if tool == 'flashdeconv':
        return deconv_components(file_manager)[1]
    kind = sample_type.of_result(file_manager, dataset_id) or sample_type.NATIVE
    st.caption(
        f"{kind} default layout: this experiment looks like a "
        f"{'native or targeted run' if kind == sample_type.NATIVE else 'lysate LC-MS/MS run'}. "
        "Press **Edit layout** to arrange your own."
    )
    return TNT_DEFAULT_LAYOUTS[kind]


def select_experiment(index):
    # A cleared selectbox reports None; show no experiment instead of failing.
    name = st.session_state.get(f'selected_experiment_dropdown_{index}')
    if name is None:
        st.session_state.pop(selection_key(index), None)
    else:
        st.session_state[selection_key(index)] = display_name_to_id[name]


def validate_selected_index(index, names):
    """Index of the stored selection in `names`, or None after dropping a
    selection that no longer exists."""
    key = selection_key(index)
    if key in st.session_state:
        display_name = display_id_to_name.get(st.session_state[key])
        if display_name in names:
            return names.index(display_name)
        del st.session_state[key]
    return None


def comparison_selected():
    """True while the first selector shows a comparison of runs."""
    return st.session_state.get(selection_key(0)) in comparison_ids


def selector(index, names, with_button=None):
    """The experiment selector of experiment `index`, the switch between its
    identification and its deconvolution, and (`with_button`) a button."""
    with st.container(horizontal=True, vertical_alignment='bottom'):
        st.selectbox(
            'choose experiment', names, key=f'selected_experiment_dropdown_{index}',
            index=validate_selected_index(index, names),
            on_change=select_experiment, args=(index,),
        )
        if st.session_state.get(selection_key(index)) in identified:
            st.segmented_control(
                'Show', [IDENTIFICATION, DECONVOLUTION], key=f'experiment_view_{index}',
                default=IDENTIFICATION, width='content',
                help='The proteoforms identified in this run, or its deconvolved masses',
            )
            # A segmented control can be cleared: keep showing something
            if st.session_state.get(f'experiment_view_{index}') is None:
                st.session_state[f'experiment_view_{index}'] = IDENTIFICATION
        if with_button is not None:
            with_button()


def experiment(index):
    """Draw experiment `index` below its selector."""
    dataset_id = st.session_state.get(selection_key(index))
    if dataset_id is None:
        return
    tool = view_of(index)
    render_grid(
        dataset_id, experiment_layout(tool, index, dataset_id), file_manager, tool,
        selection_key(index), f'flash_viewer_grid_{index}',
    )


# page initialization
params = page_setup()

file_manager = FileManager(
    st.session_state["workspace"],
    Path(st.session_state['workspace'], 'cache')
)
quant_view.migrate_layout(file_manager)

# Every deconvolved run, of which some were identified, then the comparisons
results = file_manager.get_results_list(['threedim_SN_plot'])
identified = set(file_manager.get_results_list(['protein_dfs']))
comparison_names = quant_view.comparison_names(file_manager)
comparison_ids = set(comparison_names.values())

if len(results) == 0 and not comparison_names:
    layout_button(file_manager, 'flashdeconv')
    st.error('No results to show yet. Please run a workflow first!')
    st.stop()

run_names = [file_manager.get_display_name(exp_id) for exp_id in results]
display_names = run_names + list(comparison_names)
display_name_to_id = {file_manager.get_display_name(exp_id): exp_id for exp_id in results}
display_name_to_id.update(comparison_names)
display_id_to_name = {exp_id: name for name, exp_id in display_name_to_id.items()}

# The number of experiments (and side by side) follows the layout of what the
# first experiment shows; its Edit layout edits that layout.
first_tool = view_of(0)
experiments, side_by_side = saved_layout(first_tool)
edit_layout = lambda: layout_button(file_manager, first_tool)

if comparison_selected():
    selector(0, display_names, lambda: quant_view.buttons(file_manager, st.session_state[selection_key(0)]))
    # A comparison spans its runs, so it fills the page on its own
    quant_view.render(file_manager, st.session_state[selection_key(0)])

elif len(experiments) == 2 and side_by_side:
    c1, c2 = st.columns(2)
    with c1:
        selector(0, display_names)
        experiment(0)
    with c2:
        selector(1, run_names, edit_layout)
        experiment(1)

else:
    selector(0, display_names, edit_layout)
    experiment(0)
    for exp_index in range(1, len(experiments)):
        st.divider()
        selector(exp_index, run_names)
        experiment(exp_index)

save_params(params)
