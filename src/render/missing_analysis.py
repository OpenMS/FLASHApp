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
"""

import json
import re
from pathlib import Path

import streamlit as st

from src import analysis
from src.render.layout_editor.catalog import IDENTIFICATION, ISOBARIC

CONFIGURE_PAGE = Path('content', 'TopDown', 'TopDownConfigure.py')
WORKFLOW_DIR = 'flashtnt'

ANALYSES = {
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


def available(file_manager, dataset_id):
    """The analyses run ``dataset_id`` was processed with (None: not a run)."""
    if dataset_id is None:
        return None
    has = set()
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


def options(workspace, file_manager, dataset_id, missing):
    """The dialog's choices for run ``dataset_id`` lacking the analyses ``missing``."""
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
