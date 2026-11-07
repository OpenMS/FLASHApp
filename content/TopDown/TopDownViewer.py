"""The Viewer of the merged workflow: every run in one combined view of its
deconvolution and identification, and every comparison of runs.

The view bar (the template's layout manager) picks the view all runs show:
Deconvolution, Native, Lysate or one of the user's own (+ and ✏️). Until a
view is picked, each experiment shows the default of its run: the
deconvolution default for a run deconvolved only, the Native or Lysate
identification default for an identified run, as it was searched
(src/render/topdown_layout.py). Every component of an experiment shares one
selection, so a scan or proteoform stays selected when the view changes.

Above the view bar, **Single run | Comparison** picks what is looked at: a
run of the experiment, with the views above, or a comparison of its runs,
with its own views (Label-free, TMT, iTRAQ and the user's own comparison
views, src/render/quant_view.py). Each mode lists only its own datasets and
remembers the one picked in it.
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


def select_experiment(index, widget_key):
    # A cleared selectbox reports None; show no experiment instead of failing.
    name = st.session_state.get(widget_key)
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


def selector(index, names, with_button=None, label='choose experiment', widget_key=None):
    """The experiment selector of experiment `index` and (`with_button`) a button.
    Each mode has a selectbox of its own (`widget_key`), so switching modes
    never leaves one showing a dataset of the other."""
    widget_key = widget_key or f'selected_experiment_dropdown_{index}'
    with st.container(horizontal=True, vertical_alignment='bottom'):
        st.selectbox(
            label, names, key=widget_key,
            index=validate_selected_index(index, names),
            on_change=select_experiment, args=(index, widget_key),
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
    if dataset_id not in results:
        return []
    return missing_analysis.options(workspace, file_manager, dataset_id, missing)


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
    """The view bar: the views as many as fit, + for a view of one's own, More and ✏️."""
    first = st.session_state.get(selection_key(0))
    identified, kind = run_kind(first) if first in display_id_to_name else (False, None)
    view_bar(
        topdown_layout.views_path(workspace), topdown_layout.components(file_manager),
        topdown_layout.presets(file_manager), topdown_layout.deconvolution_default(file_manager),
        key=VIEW_KEY, default_view=topdown_layout.default_view(identified, kind, file_manager),
        available=missing_analysis.available(file_manager, first if first in display_id_to_name else None),
        analyses=missing_analysis.ANALYSES, unavailable_options=missing_options,
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


RUN_MODE, COMPARISON_MODE = 'Single run', 'Comparison'
MODE_KEY = 'viewer_mode'
_LAST_SHOWN = {RUN_MODE: '_viewer_last_run', COMPARISON_MODE: '_viewer_last_comparison'}
_SHOWN_MODE = '_viewer_mode_shown'


def _mode_of(dataset_id):
    if dataset_id in comparison_ids:
        return COMPARISON_MODE
    if dataset_id in results:
        return RUN_MODE
    return None


def _switch_mode():
    """Single run | Comparison changed: show what was last picked in that mode."""
    before = st.session_state.get(_SHOWN_MODE, RUN_MODE)
    mode = st.session_state.get(MODE_KEY)
    if mode is None:  # a segmented control can be cleared: keep the mode
        st.session_state[MODE_KEY] = before
        return
    shown = st.session_state.get(selection_key(0))
    if shown is not None:
        st.session_state[_LAST_SHOWN[before]] = shown
    last = st.session_state.get(_LAST_SHOWN[mode])
    if _mode_of(last) != mode:
        pool = list(comparison_ids) if mode == COMPARISON_MODE else list(results)
        last = (analysis.newest_result(pool) or pool[0]) if pool else None
    if last is None:
        st.session_state.pop(selection_key(0), None)
    else:
        st.session_state[selection_key(0)] = last
    st.session_state[_SHOWN_MODE] = mode


def mode_switch():
    """Single run | Comparison: look at one run of the experiment or at a comparison of its runs."""
    # A selection made elsewhere (Explore results, a comparison just made) sets the mode
    mode = _mode_of(st.session_state.get(selection_key(0)))
    if mode is not None:
        st.session_state[MODE_KEY] = mode
    elif st.session_state.get(MODE_KEY) not in (RUN_MODE, COMPARISON_MODE):
        st.session_state[MODE_KEY] = RUN_MODE if results else COMPARISON_MODE
    st.session_state[_SHOWN_MODE] = st.session_state[MODE_KEY]
    st.segmented_control(
        'Show', [RUN_MODE, COMPARISON_MODE], key=MODE_KEY, on_change=_switch_mode,
        label_visibility='collapsed',
        help='**Single run**: the spectra, masses and proteoforms of one run. '
             '**Comparison**: the quantities of several runs side by side, from a start with '
             'quantification switched on.',
    )


def no_comparison():
    """Comparison mode without any comparison: say how one is made."""
    st.info(
        "No comparison of runs yet. A comparison is made when several files are processed "
        "with **Quantification** (label-free or isobaric labels) switched on in Describe experiment."
    )
    newest = (analysis.newest_result(list(results)) or list(results)[0]) if results else None
    if newest is not None and st.button(
        f'Compare {file_manager.get_display_name(newest)} with more files', icon=':material/library_add:',
        help='Opens Describe experiment with the settings of this run and label-free quantification '
             'switched on, to pick the files to compare it with',
    ):
        missing_analysis.add_files(workspace, file_manager, newest)


# Opened with Explore results on Run analysis: show what that start made
if st.session_state.pop(latest_results_key(Path(workspace, WORKFLOW_DIR)), False):
    newest = analysis.newest_result(list(results) + list(comparison_ids))
    if newest is not None:
        st.session_state[selection_key(0)] = newest

mode_switch()

if st.session_state[MODE_KEY] == COMPARISON_MODE and not comparison_names:
    no_comparison()
    save_params(params)
    st.stop()
if st.session_state[MODE_KEY] == RUN_MODE and not run_names:
    st.info('Every result here is a comparison of runs. Switch to **Comparison** to see it.')
    save_params(params)
    st.stop()

if comparison_selected():
    quant_view.view_selector(file_manager, workspace, st.session_state[selection_key(0)])
else:
    view_selector()
# The views file sets the number of experiments (and side by side); without
# one, a single experiment
count, side_by_side = topdown_layout.experiments_shown(workspace) or (1, False)
experiments = [None] * count

if st.session_state[MODE_KEY] == COMPARISON_MODE:
    selector(0, list(comparison_names), comparison_buttons, label='choose comparison',
             widget_key='selected_comparison_dropdown')
    if not comparison_selected():
        save_params(params)
        st.stop()
    # A comparison spans its runs, so it fills the page on its own
    shown = st.session_state[selection_key(0)]
    quant_view.render(file_manager, shown, quant_view.view_rows(file_manager, workspace, shown))

elif len(experiments) == 2 and side_by_side:
    c1, c2 = st.columns(2)
    with c1:
        selector(0, run_names)
        experiment(0)
    with c2:
        selector(1, run_names, lambda: edit_layout(1))
        experiment(1)

else:
    selector(0, run_names, edit_layout)
    experiment(0)
    for exp_index in range(1, len(experiments)):
        st.divider()
        selector(exp_index, run_names)
        experiment(exp_index)

save_params(params)
