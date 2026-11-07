"""The Viewer of the merged workflow: every run in one combined view of its
deconvolution and identification, and every comparison of runs.

The view bar (the template's layout manager) lists every view at once:
Deconvolution, Native, Lysate, the comparison views Label-free and TMT, iTRAQ
(src/render/quant_view.py), and the user's own (+ and ✏️). A view the shown
experiment cannot fill (identifications of a run deconvolved only, a
comparison view for one run, a run view for a comparison) stays listed,
locked, and offers what fills it. Until a view is picked, each experiment
shows the default of its run: the deconvolution default for a run
deconvolved only, the Native or Lysate identification default for an
identified run, as it was searched, the view of its kind for a comparison
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
from src.layout_manager import request_analysis, view_bar
from src.render import missing_analysis
from src.render import quant_view
from src.render.tool_results import tool_results_button
from src.common import result_download
from src.workflow.StreamlitUI import latest_results_key

TOOL = 'topdown'
# TopDownWorkflow's folder in the workspace
WORKFLOW_DIR = 'flashtnt'


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


VIEW_KEY = 'topdown_view'
MISSING_INDEX = 'topdown_missing_experiment'


def add_missing(index):
    """The button of a component whose run lacks its analysis: the view bar's dialog."""
    def on_missing(missing, component):
        st.session_state[MISSING_INDEX] = index
        request_analysis(VIEW_KEY, missing, component=component)
    return on_missing


def missing_options(missing, view):
    """What the dialog offers for the run of the experiment that asked."""
    index = st.session_state.get(MISSING_INDEX, 0) if view is None else 0
    dataset_id = st.session_state.get(selection_key(index))
    if dataset_id not in results and dataset_id not in comparison_ids:
        return []
    return missing_analysis.options(workspace, file_manager, dataset_id, missing, selection_key(index))


def experiment(index):
    """Draw experiment `index` below its selector."""
    dataset_id = st.session_state.get(selection_key(index))
    if dataset_id is None:
        return
    render_grid(
        dataset_id, experiment_layout(index, dataset_id), file_manager, TOOL,
        selection_key(index), f'flash_viewer_grid_{index}', on_missing=add_missing(index),
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
quant_view.adopt_saved_layout(file_manager, workspace)
topdown_layout.adopt_saved_layout(file_manager, workspace)
quant_view.merge_views(workspace)

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

# After Process more files, the comparison of the earlier runs with the added
# ones is shown (once), in place of the runs one by one
added = quant_view.added_comparison(file_manager, Path(workspace, missing_analysis.WORKFLOW_DIR))
if added and st.session_state.get('_viewer_added_comparison') != added:
    st.session_state['_viewer_added_comparison'] = added
    st.session_state[selection_key(0)] = added



def view_selector():
    """The view bar: every view, as many as fit, + for a view of one's own, More and ✏️."""
    first = st.session_state.get(selection_key(0))
    shown = first if first in display_id_to_name else None
    if shown in comparison_ids:
        default = quant_view.default_view(file_manager, shown)
    else:
        identified, kind = run_kind(shown) if shown is not None else (False, None)
        default = topdown_layout.default_view(identified, kind, file_manager)
    view_bar(
        topdown_layout.views_path(workspace), topdown_layout.components(file_manager),
        topdown_layout.presets(file_manager), topdown_layout.deconvolution_default(file_manager),
        key=VIEW_KEY, default_view=default,
        available=missing_analysis.available(file_manager, shown),
        analyses=missing_analysis.ANALYSES, unavailable_options=missing_options,
        dataset_label='comparison' if shown in comparison_ids else 'run',
        max_experiments=MAX_EXPERIMENTS, target='the Viewer',
        help_text=(
            "Drag a view to order the bar, untick it to keep it under More. Your own views take up "
            "to 3 components per row; ready-made views can be duplicated. Add experiments to compare "
            "runs side by side. Label-free and TMT, iTRAQ show a comparison of runs."
            + ('' if topdown_layout.has_sequence(file_manager) else
               " For a run deconvolved only, add the **Sequence view** and set a sequence there with "
               "**Change sequence**.")
        ),
    )


def add_files_button(dataset_id, runs=None):
    """Add files: process more files with the settings of the shown run (or
    comparison) and compare them with it."""
    if dataset_id is None:
        return
    if st.button('Add files', icon=':material/library_add:', key=f'add_files_{dataset_id}', width='content',
                 help='Process more files with these settings and compare them with this experiment'):
        missing_analysis.add_files(workspace, file_manager, dataset_id, runs=runs)


def comparison_buttons():
    shown = st.session_state[selection_key(0)]
    quant_view.buttons(file_manager, shown)
    add_files_button(shown, runs=quant_view.comparison_runs(file_manager, shown))


def edit_layout(index=0):
    tool_results_button()
    shown = st.session_state.get(selection_key(index))
    if index == 0 and shown in results:
        add_files_button(shown)
    if downloadable:
        result_download.viewer_button(file_manager, shown if shown in downloadable else None)


# Opened with Explore results on Run analysis: show what that start made
if st.session_state.pop(latest_results_key(Path(workspace, WORKFLOW_DIR)), False):
    newest = analysis.newest_result(list(results) + list(comparison_ids))
    if newest is not None:
        st.session_state[selection_key(0)] = newest

view_selector()
# The views file sets the number of experiments (and side by side); without
# one, a single experiment
count, side_by_side = topdown_layout.experiments_shown(workspace) or (1, False)
experiments = [None] * count

if comparison_selected():
    selector(0, display_names, comparison_buttons)
    # A comparison spans its runs, so it fills the page on its own
    shown = st.session_state[selection_key(0)]
    quant_view.render(file_manager, shown, quant_view.view_rows(file_manager, workspace, shown),
                      on_missing=add_missing(0))

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
