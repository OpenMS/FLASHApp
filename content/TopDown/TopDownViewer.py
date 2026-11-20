"""The Viewer of the merged workflow: every run in one combined view of its
deconvolution and identification, and every comparison of runs.

The view bar (the template's layout manager) picks the view all runs show:
Deconvolution, Native, Lysate or one of the user's own (+ and ✏️). Until a
view is picked, each experiment shows the default of its run: the
deconvolution default for a run deconvolved only, the Native or Lysate
identification default for an identified run, as it was searched
(src/render/topdown_layout.py). Every component of an experiment shares one
selection, so a scan or proteoform stays selected when the view changes.

**choose experiment** lists every comparison of runs and every run that
belongs to none. For a comparison, **Single run | Comparison** picks what is
looked at: one of its runs, with the run views, or the comparison, with its
own views (Identification, Quantification, QC: label-free, QC: TMT, iTRAQ
and the user's own comparison views,
src/render/quant_view.py) and **Edit comparison**, whose **Compared** column
leaves runs out (``QuantSession.hidden``). For a run of no comparison,
Comparison is locked and offers to process more files. The view bar is below
this row. **History** shows how the result came about, its settings and its
processing report. **Add analysis** processes more files with the settings
shown, or adds tool output made elsewhere.
"""


import streamlit as st

from pathlib import Path

from src.common.common import page_setup, save_params
from src.workflow.FileManager import FileManager
from src import analysis, former_examples, sample_type
from openms_insight import StateManager
from src.render.render import _state_manager, render_grid
from src.render import topdown_layout
from src.render.layout_editor.catalog import MAX_EXPERIMENTS
from src.layout_manager import request_analysis, view_bar
from src.render import missing_analysis
from src.render import quant_view
from src.render import inspect_scan
from src.render.widget_state import follow
from src.render.tool_results import tool_results_button, tool_results_dialog
from src.render.history_ui import history_button
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
        on_inspect=request_inspect,
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


MATCHER_PAGE = 'content/Tools/FragmentMatcher.py'
MATCHER_STATE = 'fragment_matcher_state'


INSPECT_REQUEST = '_viewer_inspect'


def request_inspect(dataset_id, scan_idx, proteoform_index):
    """The inspect button of a table row: open its spectrum on the next run
    (a page cannot be switched while a component is drawn)."""
    if scan_idx is None:
        return False
    st.session_state[INSPECT_REQUEST] = (dataset_id, scan_idx, proteoform_index)
    return True


def inspect_requested():
    """Open the spectrum of the row whose inspect button was clicked on Match
    fragments: its deconvolved masses, and the residues of the row's
    proteoform when it was identified in that spectrum."""
    request = st.session_state.pop(INSPECT_REQUEST, None)
    if request is None:
        return
    dataset_id, scan_idx, proteoform_index = request
    if dataset_id not in results:
        return
    data = inspect_scan.scan_input(
        file_manager, dataset_id, scan_idx, proteoform_index,
        run_name=file_manager.get_display_name(dataset_id),
    )
    if data is None:
        st.toast('This spectrum has no deconvolved masses.')
        return
    inspect_scan.store(workspace, data)
    # The matcher shows this proteoform, or keeps its own sequence
    st.session_state.pop('matcher_sequence_seen', None)
    if data['sequence']:
        StateManager(session_key=MATCHER_STATE).set_selection('sequence', data['sequence'])
    st.switch_page(MATCHER_PAGE)


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
        # A view of one's own keeps each experiment's zoom and selection (📌)
        states=lambda i: _state_manager(TOOL, selection_key(i), st.session_state.get(selection_key(i))),
        datasets=lambda i: st.session_state.get(selection_key(i)),
        help_text=(
            "Drag a view to order the bar, untick it to keep it under More. Your own views take up "
            "to 3 components per row; ready-made views can be duplicated. Add experiments to compare "
            "runs side by side."
            + ('' if topdown_layout.has_sequence(file_manager) else
               " For a run deconvolved only, add the **Sequence view** and set a sequence there with "
               "**Change sequence**.")
        ),
    )


# -- experiments: the comparisons and the runs of none ----------------------

EXPERIMENT_KEY = 'viewer_experiment'
SHOW_KEY = 'viewer_show'
RUN_KEY = 'viewer_run'
PENDING = '_viewer_pending'
LAST_RUN = '_viewer_last_run'
SHOW_COMPARISON, SHOW_RUN = 'Comparison', 'Single run'
LOCKED_KEY = 'viewer_show_locked'
LOCKED_COMPARISON = 'Comparison :material/add_circle:'


def stored_runs(comparison_id):
    """The runs of this workspace a comparison was made from (matched by
    their input file, so the example runs count too)."""
    manifest = file_manager.get_results(comparison_id, ['quant_manifest'])['quant_manifest']
    names = list((manifest.get('sources') or {}).get('runs') or []) + list(manifest.get('runs') or [])
    by_input = {missing_analysis._base_name(r): r for r in results}
    found = [r if r in results else by_input.get(missing_analysis._base_name(r)) for r in names]
    return [r for r in dict.fromkeys(found) if r is not None]


def experiment_groups():
    """Experiment name -> (the comparison or run it opens on, the runs in it):
    every comparison, then every run that belongs to none."""
    groups, grouped = {}, set()
    for name, comparison_id in comparison_names.items():
        members = stored_runs(comparison_id)
        groups[name] = (comparison_id, members)
        grouped.update(members)
    for run_id in results:
        if run_id not in grouped:
            groups[file_manager.get_display_name(run_id)] = (run_id, [run_id])
    return groups


def experiment_of(shown, groups):
    """The experiment that holds `shown` (a comparison, or a run)."""
    for name, (opens_on, members) in groups.items():
        if shown == opens_on or (opens_on in comparison_ids and shown in members):
            return name
    return None


def _on_experiment(groups):
    name = st.session_state.get(EXPERIMENT_KEY)
    if name is not None:
        st.session_state[selection_key(0)] = groups[name][0]


# The callbacks below only note what was picked: applied in the script run
# (apply_pending), whose file manager is the current one.
def _on_show(comparison_id, members):
    st.session_state[PENDING] = ('show', comparison_id, members, st.session_state.get(SHOW_KEY))


def _on_locked_comparison():
    """Comparison of a run compared with none: the dialog to process more files."""
    picked = st.session_state.get(LOCKED_KEY)
    st.session_state[LOCKED_KEY] = SHOW_RUN
    if picked != LOCKED_COMPARISON:
        return  # sent again by the browser after a switch, not a click on Comparison
    st.session_state[MISSING_INDEX] = 0
    request_analysis(VIEW_KEY, [missing_analysis.COMPARISON], component='Comparison')


def _on_run(comparison_id, members):
    st.session_state[PENDING] = ('run', comparison_id, members, st.session_state.get(RUN_KEY))


def apply_pending():
    """What was picked in Single run | Comparison, or the run."""
    pending = st.session_state.pop(PENDING, None)
    if pending is None:
        return
    what, comparison_id, members, value = pending
    shown = st.session_state.get(selection_key(0))
    if shown != comparison_id and shown not in members:
        # Picked on a comparison that is no longer shown: the browser can still
        # send the control of the experiment shown before a switch, and it
        # would take the page back there.
        return
    if what == 'show':
        if value == SHOW_RUN and members:
            last = st.session_state.get(LAST_RUN)
            st.session_state[selection_key(0)] = last if last in members else members[0]
        elif value == SHOW_COMPARISON:
            st.session_state[selection_key(0)] = comparison_id
        # a cleared segmented control keeps what is shown
    elif what == 'run' and value in display_name_to_id:
        st.session_state[selection_key(0)] = st.session_state[LAST_RUN] = display_name_to_id[value]


def add_files_menu(dataset_id, runs=None):
    """Add analysis: process more files with the settings shown, or add the
    output of FLASHDeconv and FLASHTnT run elsewhere."""
    with st.popover('Add analysis', icon=':material/library_add:', width='content'):
        if dataset_id is not None and st.button(
            'Process more files', icon=':material/play_circle:', key=f'add_files_{dataset_id}', width='stretch',
            help='Process more files with these settings and compare them with this experiment',
        ):
            missing_analysis.add_files(workspace, file_manager, dataset_id, runs=runs)
        if st.button('Add tool results', icon=':material/upload_file:', key='add_tool_results', width='stretch',
                     help='Show the output of FLASHDeconv and FLASHTnT run outside the app'):
            tool_results_dialog()


def experiment_bar():
    """choose experiment; for a comparison, Single run | Comparison with its
    run, or Edit comparison."""
    groups = experiment_groups()
    shown = st.session_state.get(selection_key(0))
    name = experiment_of(shown, groups)
    if name is None:
        name = next(iter(groups))
        shown = st.session_state[selection_key(0)] = groups[name][0]
    follow(EXPERIMENT_KEY, name)
    opens_on, members = groups[name]
    comparison = opens_on in comparison_ids
    with st.container(horizontal=True, vertical_alignment='bottom'):
        st.selectbox('choose experiment', list(groups), key=EXPERIMENT_KEY,
                     on_change=_on_experiment, args=(groups,))
        if comparison:
            follow(SHOW_KEY, SHOW_COMPARISON if shown == opens_on else SHOW_RUN)
            st.segmented_control(
                'show', [SHOW_RUN, SHOW_COMPARISON], key=SHOW_KEY, on_change=_on_show, args=(opens_on, members),
                help='**Single run**: the spectra, masses and proteoforms of one of its runs. '
                     '**Comparison**: the quantities of its runs side by side.',
            )
            if shown != opens_on:
                names = [file_manager.get_display_name(r) for r in members]
                follow(RUN_KEY, display_id_to_name.get(shown))
                st.selectbox('run', names, key=RUN_KEY, on_change=_on_run, args=(opens_on, members))
        else:
            follow(LOCKED_KEY, SHOW_RUN)
            st.segmented_control(
                'show', [SHOW_RUN, LOCKED_COMPARISON], key=LOCKED_KEY, on_change=_on_locked_comparison,
                help='This run is compared with no other run: process more files with its settings '
                     'to compare them.',
            )
        if comparison and shown == opens_on:
            quant_view.buttons(file_manager, opens_on)
        history_button(file_manager, shown, workspace)
        add_files_menu(shown, runs=members if shown == opens_on and comparison else None)
        if shown not in comparison_ids and downloadable:
            result_download.viewer_button(file_manager, shown if shown in downloadable else None)


# Opened with Explore results on Run analysis: show what that start made
if st.session_state.pop(latest_results_key(Path(workspace, WORKFLOW_DIR)), False):
    newest = analysis.newest_result(list(results) + list(comparison_ids))
    if newest is not None:
        st.session_state[selection_key(0)] = newest

apply_pending()
inspect_requested()
# Nothing (or something gone) shown: the first experiment
if st.session_state.get(selection_key(0)) not in display_id_to_name:
    first_groups = experiment_groups()
    st.session_state[selection_key(0)] = next(iter(first_groups.values()))[0]

# The views file sets the number of experiments (and side by side); without
# one, a single experiment
count, side_by_side = topdown_layout.experiments_shown(workspace) or (1, False)
experiments = [None] * count

# The experiment first, then the views of what it shows below it
experiment_bar()
if comparison_selected():
    quant_view.view_selector(file_manager, workspace, st.session_state[selection_key(0)])
else:
    view_selector()

if comparison_selected():
    # A comparison spans its runs, so it fills the page on its own
    shown = st.session_state[selection_key(0)]
    quant_view.render(file_manager, shown, quant_view.view_rows(file_manager, workspace, shown))

elif len(experiments) == 2 and side_by_side:
    c1, c2 = st.columns(2)
    with c1:
        experiment(0)
    with c2:
        selector(1, run_names)
        experiment(1)

else:
    experiment(0)
    for exp_index in range(1, len(experiments)):
        st.divider()
        selector(exp_index, run_names)
        experiment(exp_index)

save_params(params)
