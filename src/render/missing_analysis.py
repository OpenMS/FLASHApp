"""What the Viewer offers when a view needs an analysis the shown run was not processed with.

A run deconvolved only has no identifications, a run without isobaric labels no
reporter ions. Views needing them are locked in the view bar (the components'
``needs`` in ``layout_editor/catalog.py``), and a click, or the **Add …** button
of such a component, opens a dialog that offers one thing: process more files
with this run's settings and the missing analysis switched on. Only the added
files are processed, and this run is compared with them in one comparison,
computed in the Viewer, as it is if it was processed with what the new start
does (else its file is processed again with them, so picking no new file runs
this run again with the analysis).

The bar lists the views of a comparison of runs too: a run shown alone offers
the comparisons it is part of and to process more files to compare it with,
and a comparison shown offers its runs for the views of one run.
"""

import json
import re
from pathlib import Path

import streamlit as st

from src import analysis
from src.render.layout_editor.catalog import COMPARISON, IDENTIFICATION, ISOBARIC, RUN

CONFIGURE_PAGE = Path('content', 'TopDown', 'TopDownConfigure.py')
WORKFLOW_DIR = 'flashtnt'

ANALYSES = {
    RUN: {
        'label': 'one run',
        'missing': 'This view shows one run, and a comparison spans several: show one of its runs.',
        'action': 'Show a run',
    },
    COMPARISON: {
        'label': 'a comparison of runs',
        'action': 'Compare runs',
        'help': ('FLASHQuantConsensus compares the features of two or more runs (**Quantification: '
                 'Label-free** or **Isobaric labels** on the Describe experiment page).'),
    },
    IDENTIFICATION: {
        'label': 'identification',
        'help': ('FLASHTnT searches a protein database for the proteoforms of the deconvolved masses '
                 '(**Identify proteoforms** on the Describe experiment page).'),
    },
    ISOBARIC: {
        'label': 'isobaric labels',
        'help': ('FLASHDeconv measures the reporter ions of TMT or iTRAQ labels '
                 '(**Quantification: Isobaric labels** on the Describe experiment page).'),
    },
}

# What switching an analysis on changes in the workflow's parameters
_SWITCH_ON = {
    IDENTIFICATION: {analysis.IDENTIFY_KEY: True},
    ISOBARIC: {analysis.QUANT_KEY: analysis.ISOBARIC},
}


def is_comparison(file_manager, dataset_id):
    return file_manager.result_exists(dataset_id, 'quant_manifest')


def available(file_manager, dataset_id):
    """What experiment ``dataset_id`` can show: ``run`` and the analyses a run
    was processed with, ``comparison`` for a comparison of runs (None: nothing shown)."""
    if dataset_id is None:
        return None
    if is_comparison(file_manager, dataset_id):
        return {COMPARISON}
    has = {RUN}
    if file_manager.result_exists(dataset_id, 'protein_dfs'):
        has.add(IDENTIFICATION)
    if file_manager.result_exists(dataset_id, 'isobaric_dfs'):
        has.add(ISOBARIC)
    return has


# The example runs of the demo workspace and the input file each was made from
# (``EXAMPLES`` in ``example-data/build_examples.py``)
EXAMPLE_INPUTS = {
    'example_deconvolution': 'Deconvolution-only_ProteinMix_LC-MS',
    'example_native_identified': 'Native-identified_GapA_HCD',
    'example_lysate_identified': 'Label-free_SpikeIn-1x_R1',
    'example_isobaric': 'Isobaric_TMT6_Ecoli-Ratio',
}

_TIME_SUFFIX = re.compile(r'_\d{8}-\d{6}$')


def _base_name(dataset_id):
    """The input file's name without extension: dataset ids are ``<name>_<time>``."""
    if dataset_id in EXAMPLE_INPUTS:
        return EXAMPLE_INPUTS[dataset_id]
    return _TIME_SUFFIX.sub('', dataset_id)


def input_file(workspace, dataset_id):
    """The workflow's input file run ``dataset_id`` was made from, if it is still there."""
    folder = Path(workspace, WORKFLOW_DIR, 'input-files', 'mzML-files')
    candidates = []
    if folder.is_dir():
        candidates = [p for p in folder.iterdir() if p.name != 'external_files.txt']
        external = Path(folder, 'external_files.txt')
        if external.exists():
            candidates += [Path(line.strip()) for line in external.read_text().splitlines() if line.strip()]
    base = _base_name(dataset_id)
    for path in candidates:
        if Path(str(path).rstrip('/\\')).stem == base:
            return str(path)
    return None


def _current_files(workspace):
    try:
        params = json.loads(Path(workspace, WORKFLOW_DIR, 'params.json').read_text(encoding='utf-8'))
    except (OSError, ValueError):
        return []
    files = params.get('mzML-files') or []
    return [files] if isinstance(files, str) else list(files)


def run_settings(workspace, file_manager, dataset_id):
    """The workflow's parameters as run ``dataset_id`` was made with them.

    The workflow's current settings (database, samples), with the FLASHDeconv
    and FLASHTnT parameters stored with the run and its input file.
    """
    try:
        params = json.loads(Path(workspace, WORKFLOW_DIR, 'params.json').read_text(encoding='utf-8'))
    except (OSError, ValueError):
        params = {}
    for tool, name in (('FLASHDeconv', 'FD_parameters'), ('FLASHTnT', 'FTnT_parameters')):
        if file_manager.result_exists(dataset_id, name):
            params[tool] = dict(file_manager.get_results(dataset_id, [name])[name] or {})
    params[analysis.IDENTIFY_KEY] = IDENTIFICATION in (available(file_manager, dataset_id) or set())
    source = input_file(workspace, dataset_id)
    if source is not None:
        params['mzML-files'] = [source]
    return params, source is not None


def more_files_settings(params, changes, reusable):
    """The settings, changes and runs added to of Process more files.

    The added files are compared with this run, so the start quantifies
    (label-free unless the run has labels). ``reusable(identify, quant, plex)``
    says if the run can join that comparison as it is: then only the files
    picked now are processed; else its own file is processed again with them.
    Returns (params, changes, adds_to).
    """
    changes = dict(changes)
    merged = {**params, **changes}
    if analysis.quantification(merged) == analysis.NO_QUANT:
        changes[analysis.QUANT_KEY] = analysis.LABEL_FREE
        merged[analysis.QUANT_KEY] = analysis.LABEL_FREE
    params = dict(params)
    if reusable(merged.get(analysis.IDENTIFY_KEY, True), analysis.quantification(merged), analysis.plex(merged)):
        params['mzML-files'] = []
        return params, changes, True
    return params, changes, False


def _follow_up(workspace, file_manager, dataset_id, changes, note, more_files=False, runs=None):
    """Open Describe experiment with the settings of ``dataset_id`` and ``changes``.

    ``more_files`` adds the files picked there to its experiment: the runs
    ``runs`` (default the run itself; the runs of a comparison) join the new
    comparison as they are when they fit the settings, else their files are
    processed again with the new ones.
    """
    from src.Workflow import TopDownWorkflow  # the page set up the workspace already
    from src.quant.run import reusable_runs

    runs = list(runs or [dataset_id])
    params, _ = run_settings(workspace, file_manager, runs[0])
    name = file_manager.get_display_name(dataset_id)
    if not more_files:
        TopDownWorkflow().start_from(params, changes, name=name, page=CONFIGURE_PAGE, note=note)
        return
    own = []
    for run in runs:
        source = input_file(workspace, run)
        if source is not None:
            own.append(source)
    own = list(dict.fromkeys(own))
    params, changes, adds = more_files_settings(
        params, changes,
        lambda identify, quant, plex: len(reusable_runs(file_manager, runs, identify, quant, plex)) == len(runs),
    )
    if not adds:
        # The experiment's files are processed again with the files added in the Browse dialog
        params['mzML-files'] = list(dict.fromkeys(_current_files(workspace) + own))
    TopDownWorkflow().start_from(
        params, changes, name=name, page=CONFIGURE_PAGE, note=note,
        adds_to=runs if adds else None, files=own if adds else None,
        files_key='mzML-files' if adds else None,
    )


def add_files(workspace, file_manager, dataset_id, runs=None):
    """Process more files with the settings of ``dataset_id`` (a run, or a
    comparison of ``runs``) and compare them with it."""
    _follow_up(workspace, file_manager, dataset_id, {}, 'to add more files', more_files=True, runs=runs)


def comparison_members(file_manager, comparison):
    """The runs of the Viewer (dataset ids) that comparison ``comparison`` compares.

    The manifest names the runs FLASHQuantConsensus read (``sources``; the
    fractions of a labelled set); an example run is found by its input file.
    """
    manifest = file_manager.get_results(comparison, ['quant_manifest'])['quant_manifest']
    sources = manifest.get('sources') if isinstance(manifest.get('sources'), dict) else {}
    compared = list(sources.get('runs') or manifest.get('fractions') or manifest.get('runs') or [])
    runs = file_manager.get_results_list(['threedim_SN_plot'])
    out = []
    for source in compared:
        base = _TIME_SUFFIX.sub('', str(source))
        for run in runs:
            if run not in out and (run == source or _base_name(run) == base):
                out.append(run)
    return out


def _show(selection_key, dataset_id):
    def show():
        st.session_state[selection_key] = dataset_id
    return show


def options(workspace, file_manager, dataset_id, missing, selection_key=None):
    """The dialog's choices for experiment ``dataset_id`` lacking the analyses ``missing``.

    ``selection_key`` is the session key of the experiment's selector: the
    choices that show another experiment (a run of the comparison, a
    comparison the run is part of) set it.
    """
    from src.render import quant_view
    if RUN in missing:
        # A comparison: the views of one run show one of its runs
        out = []
        for run in comparison_members(file_manager, dataset_id) if selection_key is not None else []:
            out.append({'label': f'Show {file_manager.get_display_name(run)}', 'icon': ':material/visibility:',
                        'on_click': _show(selection_key, run)})
        return out
    out = []
    if COMPARISON in missing:
        if selection_key is not None:
            for name, comparison in quant_view.comparison_names(file_manager).items():
                if dataset_id in comparison_members(file_manager, comparison):
                    out.append({'label': f'Show {name}', 'icon': ':material/visibility:',
                                'help': f'A comparison {file_manager.get_display_name(dataset_id)} is part of.',
                                'on_click': _show(selection_key, comparison)})
        missing = [m for m in missing if m != COMPARISON]
        if not missing:
            # Process more files compares them with this run (label-free unless it has labels)
            name = file_manager.get_display_name(dataset_id)
            _, has_file = run_settings(workspace, file_manager, dataset_id)
            out.append({
                'label': 'Process more files with these settings',
                'icon': ':material/library_add:',
                'help': (f"Opens Describe experiment with the settings of {name}, to pick files to compare "
                         f"with it. The Viewer then shows {name} and the added files in one comparison."
                         + ('' if has_file else " Its input file is no longer in the workspace: pick it again.")),
                'on_click': lambda: add_files(workspace, file_manager, dataset_id),
            })
            return out
    return out + _process_more(workspace, file_manager, dataset_id, missing)


def _process_more(workspace, file_manager, dataset_id, missing):
    """Process more files with the run's settings and ``missing`` switched on."""
    labels = ' and '.join(ANALYSES[m]['label'].lower() for m in missing)
    changes = {}
    for name in missing:
        changes.update(_SWITCH_ON.get(name, {}))
    _, has_file = run_settings(workspace, file_manager, dataset_id)
    name = file_manager.get_display_name(dataset_id)
    return [{
        'label': 'Process more files with these settings',
        'icon': ':material/library_add:',
        'help': (f"Opens Describe experiment with the settings of {name} and {labels} switched on, to pick "
                 f"files to add. The Viewer then compares them with {name} in one comparison, computed "
                 "here; a run processed without these settings is processed again with them."
                 + ('' if has_file else " Its input file is no longer in the workspace: pick it again.")
                 + (" Check the protein database there." if IDENTIFICATION in missing else "")),
        'on_click': lambda: _follow_up(workspace, file_manager, dataset_id, changes,
                                       f'with {labels} switched on, for more files', more_files=True),
    }]
