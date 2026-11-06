"""The Viewer of the merged workflow: every run in one combined view of its
deconvolution and identification, and every comparison of runs.

The view bar (the template's layout manager) picks the view all runs show:
Deconvolution, Native, Lysate or one of the user's own (+ and ✏️). Until a
view is picked, each experiment shows the default of its run: the
deconvolution default for a run deconvolved only, the Native or Lysate
identification default for an identified run, as it was searched
(src/render/topdown_layout.py). Every component of an experiment shares one
selection, so a scan or proteoform stays selected when the view changes.
"""

import streamlit as st

from pathlib import Path

from src.common.common import page_setup, save_params
from src.workflow.FileManager import FileManager
from src import analysis, former_examples, sample_type
from src.render.render import render_grid
from src.render import topdown_layout
from src.render.layout_editor.catalog import MAX_EXPERIMENTS
from src.layout_manager import view_bar
from src.render import quant_view
from src.render.tool_results import tool_results_button
from src.common import result_download

TOOL = 'topdown'


def selection_key(index):
    return f'selected_experiment{index}'


def run_kind(dataset_id):
    """(identified, sample type) of a run."""
    identified = dataset_id in identified_runs
    kind = (sample_type.of_result(file_manager, dataset_id) or sample_type.NATIVE) if identified else None
    return identified, kind


def experiment_layout(index, dataset_id):
    """The rows experiment `index` shows: the picked view, else the default of its run."""
    identified, kind = run_kind(dataset_id)
    labelled = dataset_id in labelled_runs
    picked = topdown_layout.view_rows(file_manager, workspace, identified, kind, index, labelled)
    if picked is not None:
        return picked
    if topdown_layout.has_legacy(file_manager, identified):
        name = 'identified runs' if identified else 'runs deconvolved only'
        text = f"The layout saved earlier for {name} (the view **{topdown_layout.LEGACY_VIEWS[identified][1]}**)."
    elif identified:
        sample = analysis.SAMPLE_LABELS[kind].lower()
        text = f"Default layout for {sample}, as this run was searched."
    else:
        text = "Default layout for a run deconvolved only."
    if labelled:
        text += " The reporter ions of its labels are shown below the spectra it selects."
    st.caption(text + " Pick a view above, or make your own with **+**.")
    return topdown_layout.default_rows(file_manager, identified, kind, index, labelled)


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
workspace = Path(st.session_state['workspace'])
former_examples.add_example_results(Path(st.session_state['workspace'], 'cache'))
file_manager = FileManager(
    st.session_state["workspace"],
    Path(st.session_state['workspace'], 'cache')
)
quant_view.migrate_layout(file_manager)
topdown_layout.adopt_saved_layout(file_manager, workspace)

# Every deconvolved run, of which some were identified, then the comparisons
results = file_manager.get_results_list(['threedim_SN_plot'])
identified_runs = set(file_manager.get_results_list(['protein_dfs']))
labelled_runs = set(file_manager.get_results_list(['isobaric_dfs']))
comparison_names = quant_view.comparison_names(file_manager)
comparison_ids = set(comparison_names.values())
downloadable = set(result_download.downloadable_runs(file_manager))

if len(results) == 0 and not comparison_names:
    with st.container(horizontal=True):
        tool_results_button()
    st.error('No results to show yet. Run the analysis first, or add the output of the tools with **Add tool results**.')
    st.stop()

run_names = [file_manager.get_display_name(exp_id) for exp_id in results]
display_names = run_names + list(comparison_names)
display_name_to_id = {file_manager.get_display_name(exp_id): exp_id for exp_id in results}
display_name_to_id.update(comparison_names)
display_id_to_name = {exp_id: name for name, exp_id in display_name_to_id.items()}



def view_selector():
    """The view bar: the views as many as fit, + for a view of one's own, More and ✏️."""
    first = st.session_state.get(selection_key(0))
    identified, kind = run_kind(first) if first in display_id_to_name else (False, None)
    view_bar(
        topdown_layout.views_path(workspace), topdown_layout.components(file_manager),
        topdown_layout.presets(file_manager), topdown_layout.deconvolution_default(file_manager),
        key='topdown_view', default_view=topdown_layout.default_view(identified, kind, file_manager),
        max_experiments=MAX_EXPERIMENTS, target='the Viewer',
        help_text=(
            "Drag a view to order the bar, untick it to keep it under More. Your own views take up "
            "to 3 components per row; ready-made views can be duplicated. Add experiments to compare "
            "runs side by side."
            + ('' if topdown_layout.has_sequence(file_manager) else
               " For a run deconvolved only, add the **Sequence view** and set a sequence there with "
               "**Change sequence**.")
        ),
    )


def edit_layout(index=0):
    tool_results_button()
    shown = st.session_state.get(selection_key(index))
    if downloadable:
        result_download.viewer_button(file_manager, shown if shown in downloadable else None)


if not comparison_selected():
    view_selector()
# The views file sets the number of experiments (and side by side); without
# one, a single experiment
count, side_by_side = topdown_layout.experiments_shown(workspace) or (1, False)
experiments = [None] * count

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
