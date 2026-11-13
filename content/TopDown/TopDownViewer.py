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
from src.layout_manager import DatasetKind, ViewerSpec, show_dataset, show_viewer
from src.render import missing_analysis
from src.render import quant_view
from src.render import inspect_scan
from src.render.tool_results import tool_results_button, tool_results_dialog
from src.render.history_ui import history_button
from src.common import result_download
from src.workflow.StreamlitUI import latest_results_key

TOOL = 'topdown'
# TopDownWorkflow's folder in the workspace
WORKFLOW_DIR = 'flashtnt'
KEY = 'topdown_viewer'
VIEW_KEY = 'topdown_view'


def selection_key(index):
    """The identifier of experiment `index`'s selections (its Insight state)."""
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
    st.session_state.pop('matcher_seen_sequence', None)
    if data['sequence']:
        StateManager(session_key=MATCHER_STATE).set_selection('sequence', data['sequence'])
    st.switch_page(MATCHER_PAGE)


# -- experiments: the comparisons and the runs of none ----------------------

def stored_runs(comparison_id):
    """The runs of this workspace a comparison was made from (matched by
    their input file, so the example runs count too)."""
    manifest = file_manager.get_results(comparison_id, ['quant_manifest'])['quant_manifest']
    names = list((manifest.get('sources') or {}).get('runs') or []) + list(manifest.get('runs') or [])
    by_input = {missing_analysis._base_name(r): r for r in results}
    found = [r if r in results else by_input.get(missing_analysis._base_name(r)) for r in names]
    return [r for r in dict.fromkeys(found) if r is not None]


def experiments():
    """Experiment id -> name, and experiment id -> its runs: every comparison,
    then every run that belongs to none."""
    names, members, grouped = {}, {}, set()
    for name, comparison_id in comparison_names.items():
        names[comparison_id] = name
        members[comparison_id] = stored_runs(comparison_id)
        grouped.update(members[comparison_id])
    for run_id in results:
        if run_id not in grouped:
            names[run_id] = file_manager.get_display_name(run_id)
    return names, members


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


def top_buttons(ctx):
    """What acts on the experiment shown, beside its selector: Edit comparison,
    History, Add analysis and Export."""
    comparison = ctx.dataset in comparison_ids
    if comparison:
        quant_view.buttons(file_manager, ctx.dataset)
    history_button(file_manager, ctx.dataset, workspace)
    add_files_menu(ctx.dataset, runs=ctx.members if comparison else None)
    if not comparison and downloadable:
        result_download.viewer_button(file_manager, ctx.dataset if ctx.dataset in downloadable else None)


def no_results():
    with st.container(horizontal=True):
        tool_results_button()
    st.error('No results to show yet. Run the analysis first, or add the output of the tools with **Add tool results**.')


RUNS = DatasetKind(
    name='run', components=topdown_layout.components(file_manager), presets=topdown_layout.presets(file_manager),
    default_layout=topdown_layout.deconvolution_default(file_manager),
    views_path=topdown_layout.views_path(workspace), view_key=VIEW_KEY,
    # Until a view is picked, each experiment shows the default of its run
    default_view=lambda ctx: topdown_layout.default_view(*run_kind(ctx.dataset), file_manager),
    rows=lambda ctx, rows, saved: experiment_layout(ctx.index, ctx.dataset),
    render=lambda ctx, rows: render_grid(
        ctx.dataset, rows, file_manager, TOOL, selection_key(ctx.index), f'flash_viewer_grid_{ctx.index}',
        on_missing=lambda missing, component, c=ctx: c.request_analysis(missing, component),
        on_inspect=request_inspect,
    ),
    # A view of one's own keeps each experiment's zoom and selection (📌)
    states=lambda ctx: _state_manager(TOOL, selection_key(ctx.index), ctx.dataset),
    available=lambda ctx: missing_analysis.available(file_manager, ctx.dataset),
    analyses=missing_analysis.ANALYSES,
    unavailable_options=lambda missing, view, ctx: (
        missing_analysis.options(workspace, file_manager, ctx.dataset, missing) if ctx.dataset in results else []),
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


def _before_comparison(ctx):
    quant_view.migrate_views(quant_view.views_path(workspace))
    quant_view.fit_current(quant_view.views_path(workspace), quant_view.comparison_presets(file_manager, ctx.dataset))


# A comparison spans its runs, so it fills the page on its own
COMPARISONS = DatasetKind(
    name='comparison', components=quant_view.QUANT_COMPONENTS,
    # The QC of its kind only (label-free, or TMT, iTRAQ)
    presets=lambda ctx: quant_view.comparison_presets(file_manager, ctx.dataset),
    default_layout=quant_view.QUANT_QUANTIFICATION_LAYOUT,
    views_path=quant_view.views_path(workspace), view_key=quant_view.VIEW_KEY, before=_before_comparison,
    default_view=lambda ctx: quant_view.default_view(file_manager, ctx.dataset),
    rows=lambda ctx, rows, saved: quant_view.view_rows(file_manager, workspace, ctx.dataset),
    render=lambda ctx, rows: quant_view.render(file_manager, ctx.dataset, rows),
    states=None, max_experiments=1,
    available=lambda ctx: quant_view.available(file_manager, ctx.dataset),
    analyses=missing_analysis.ANALYSES,
    unavailable_options=lambda missing, view, ctx: quant_view._options(workspace, file_manager, ctx.dataset)(
        missing, view),
    dataset_label='comparison', target='comparisons of runs',
    help_text=(
        "Drag a view to order the bar, untick it to keep it under More. These views are used for "
        "every comparison of runs; your own take up to 3 components per row, ready-made views can "
        "be duplicated."
    ),
)

names, members = experiments()

# After Process more files, the comparison of the earlier runs with the added
# ones is shown (once), in place of the runs one by one
added = quant_view.added_comparison(file_manager, Path(workspace, missing_analysis.WORKFLOW_DIR))
if added and st.session_state.get('_viewer_added_comparison') != added:
    st.session_state['_viewer_added_comparison'] = added
    show_dataset(KEY, added)

# Opened with Explore results on Run analysis: show what that start made
if st.session_state.pop(latest_results_key(Path(workspace, WORKFLOW_DIR)), False):
    newest = analysis.newest_result(list(results) + list(comparison_ids))
    if newest is not None:
        show_dataset(KEY, newest)

inspect_requested()

show_viewer(ViewerSpec(
    key=KEY, kinds=[RUNS, COMPARISONS], datasets=names,
    kind_of=lambda dataset: 'comparison' if dataset in comparison_ids else 'run',
    # A comparison is shown as a whole or one of its runs (Single run | Comparison)
    member_mode='switch', members=lambda experiment: members.get(experiment, []),
    whole=lambda experiment: experiment in comparison_ids,
    member_format=file_manager.get_display_name,
    switch_help='**Single run**: the spectra, masses and proteoforms of one of its runs. '
                '**Comparison**: the quantities of its runs side by side.',
    # A run compared with no other run: process more files with its settings
    lone_whole=lambda ctx: ctx.request_analysis([missing_analysis.COMPARISON], 'Comparison'),
    lone_help='This run is compared with no other run: process more files with its settings to compare them.',
    label='choose experiment',
    # Later experiments show runs, side by side or one below the other
    later=lambda kind: {run: file_manager.get_display_name(run) for run in results},
    buttons=top_buttons, empty=no_results,
))

save_params(params)
