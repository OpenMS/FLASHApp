"""The Viewer of the merged workflow: every run in one combined view of its
deconvolution and identification, and every comparison of runs.

All runs share one layout (Edit layout, ``topdown_layout``). Until one is
saved, each experiment starts from the default of its run: the
deconvolution default for a run deconvolved only, the Native or Lysate
identification default for an identified run, as it was searched
(src/render/topdown_layout.py).
"""

import streamlit as st

from pathlib import Path

from src.common.common import page_setup, save_params
from src.workflow.FileManager import FileManager
from src import analysis, former_examples, sample_type
from src.render.render import render_grid
from src.render import topdown_layout
from src.render.layout_editor.dialogs import layout_button
from src.render import quant_view
from src.render.tool_results import tool_results_button
from src.common import result_download

TOOL = 'topdown'


def selection_key(index):
    return f'selected_experiment{index}'


def experiment_layout(index, dataset_id):
    """The rows experiment `index` shows: the saved layout, else the default of its run."""
    identified = dataset_id in identified_runs
    labelled = dataset_id in labelled_runs
    kind = (sample_type.of_result(file_manager, dataset_id) or sample_type.NATIVE) if identified else None
    if topdown_layout.saved(file_manager) is None:
        if topdown_layout.has_legacy(file_manager, identified):
            name = 'identified runs' if identified else 'deconvolved runs'
            text = f"The layout saved earlier for {name}."
        elif identified:
            sample = analysis.SAMPLE_LABELS[kind].lower()
            text = f"Default layout for {sample}, as this run was searched."
        else:
            text = "Default layout for a run deconvolved only."
        if labelled:
            text += " The reporter ions of its labels are shown below the spectra it selects."
        st.caption(text + " Press **Edit layout** to arrange your own.")
    return topdown_layout.rows(file_manager, identified, kind, index, labelled)


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
    """The experiment selector of experiment `index` and (`with_button`) a button."""
    with st.container(horizontal=True, vertical_alignment='bottom'):
        st.selectbox(
            'choose experiment', names, key=f'selected_experiment_dropdown_{index}',
            index=validate_selected_index(index, names),
            on_change=select_experiment, args=(index,),
        )
        if with_button is not None:
            with_button()


def experiment(index):
    """Draw experiment `index` below its selector."""
    dataset_id = st.session_state.get(selection_key(index))
    if dataset_id is None:
        return
    render_grid(
        dataset_id, experiment_layout(index, dataset_id), file_manager, TOOL,
        selection_key(index), f'flash_viewer_grid_{index}',
    )


# page initialization
params = page_setup()

# A workspace from an earlier release gets today's example runs (once)
former_examples.add_example_results(Path(st.session_state['workspace'], 'cache'))
file_manager = FileManager(
    st.session_state["workspace"],
    Path(st.session_state['workspace'], 'cache')
)
quant_view.migrate_layout(file_manager)

# Every deconvolved run, of which some were identified, then the comparisons
results = file_manager.get_results_list(['threedim_SN_plot'])
identified_runs = set(file_manager.get_results_list(['protein_dfs']))
labelled_runs = set(file_manager.get_results_list(['isobaric_dfs']))
comparison_names = quant_view.comparison_names(file_manager)
comparison_ids = set(comparison_names.values())
downloadable = set(result_download.downloadable_runs(file_manager))

if len(results) == 0 and not comparison_names:
    with st.container(horizontal=True):
        layout_button(file_manager)
        tool_results_button()
    st.error('No results to show yet. Run the analysis first, or add the output of the tools with **Add tool results**.')
    st.stop()

run_names = [file_manager.get_display_name(exp_id) for exp_id in results]
display_names = run_names + list(comparison_names)
display_name_to_id = {file_manager.get_display_name(exp_id): exp_id for exp_id in results}
display_name_to_id.update(comparison_names)
display_id_to_name = {exp_id: name for name, exp_id in display_name_to_id.items()}

# The saved layout sets the number of experiments (and side by side); without
# one, a single experiment
experiments, side_by_side = topdown_layout.saved(file_manager) or ([None], False)


def edit_layout(index=0):
    layout_button(file_manager)
    tool_results_button()
    shown = st.session_state.get(selection_key(index))
    if downloadable:
        result_download.viewer_button(file_manager, shown if shown in downloadable else None)


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
        selector(1, run_names, lambda: edit_layout(1))
        experiment(1)

else:
    selector(0, display_names, edit_layout)
    experiment(0)
    for exp_index in range(1, len(experiments)):
        st.divider()
        selector(exp_index, run_names)
        experiment(exp_index)

save_params(params)
