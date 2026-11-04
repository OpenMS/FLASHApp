"""The comparison step of a workflow run (Streamlit-free: it runs in the
workflow process or the queue worker).

After FLASHDeconv (with -quant) and, for FLASHTnT, the search have run on
every file, FLASHQuantConsensus links the features of the runs and, with the
PrSMs, says which proteoform each linked feature is. The result is stored as
a comparison dataset next to the runs, which the Viewer lists.
"""

import time
from shutil import rmtree
from os.path import basename, join, splitext
from pathlib import Path

from src.quant import contract, isobaric, session
from src.workflow import ExperimentalDesign

TOOL = 'FLASHQuantConsensus'
# Written by the workflow itself
CONSENSUS_EXCLUDE = ['in', 'in_feature', 'in_id', 'out', 'out_id', 'out_mbr_detail']


def run_name(path):
    """A run's name in FLASHQuantConsensus's output: the base name of its mzML."""
    return splitext(basename(str(path)))[0]


def writes_identifications(parameter_manager):
    """True when the installed FLASHQuantConsensus has -out_id."""
    if not parameter_manager.create_ini(TOOL):
        return False
    ini = Path(parameter_manager.ini_dir, f'{TOOL}.ini')
    return 'name="out_id"' in ini.read_text(errors='replace')


def design_of(workflow_dir, selected_files, mzmls):
    """run -> (condition, replicate) from the workflow's Samples table. The
    table lists the selected files (perhaps .raw or .d); runs are matched by
    base name."""
    table = ExperimentalDesign.load_design(workflow_dir, [str(f) for f in selected_files])
    by_name = {
        run_name(row.file): (row.condition, int(row.replicate)) for row in table.itertuples()
    }
    return {run_name(m): by_name.get(run_name(m), ('', i + 1)) for i, m in enumerate(mzmls)}


def database_of(workflow):
    """The FASTA the runs were searched against (None without one)."""
    try:
        return str(workflow.file_manager.get_files(workflow.params['fasta-file'])[0])
    except (KeyError, ValueError, IndexError, TypeError):
        return None


def compare_runs(workflow, tool, dataset_ids, mzmls, selected_files, prsm_key=None):
    """Run FLASHQuantConsensus over the stored runs `dataset_ids` (with their
    `mzmls`, in the same order) and store the comparison. `prsm_key` names the
    stored PrSM tsv of each run (FLASHTnT) to add proteoforms and proteins.

    Returns the comparison's dataset id, or None when the step failed.
    """
    file_manager, executor, logger = workflow.file_manager, workflow.executor, workflow.logger
    features = [
        str(file_manager.get_results(dataset_id, ['out_tsv'])['out_tsv']) for dataset_id in dataset_ids
    ]
    dataset_id = f"comparison_{time.strftime('%Y%m%d-%H%M%S')}"
    folder = Path(workflow.workflow_dir, 'comparison', dataset_id)
    folder.mkdir(parents=True, exist_ok=True)
    out = join(folder, 'consensus.tsv')
    io = {'in': [[str(m) for m in mzmls]], 'in_feature': [features], 'out': [out]}

    out_id = None
    if prsm_key is not None:
        if writes_identifications(executor.parameter_manager):
            io['in_id'] = [[
                str(file_manager.get_results(d, [prsm_key])[prsm_key]) for d in dataset_ids
            ]]
            out_id = join(folder, 'identifications.tsv')
            io['out_id'] = [out_id]
        else:
            logger.log(f'This {TOOL} cannot assign proteoforms (no -out_id): comparing features only.')

    logger.log(f'-> Linking the features of {len(mzmls)} runs ({TOOL})...')
    if not executor.run_topp(TOOL, input_output=io):
        return None

    logger.log('-> Storing the comparison...')
    levels, runs = contract.build(out, out_id, database_of(workflow))
    file_manager.store_file(dataset_id, 'consensus_tsv', Path(out), file_name='consensus.tsv')
    if out_id is not None and Path(out_id).exists():
        file_manager.store_file(dataset_id, 'identifications_tsv', Path(out_id), file_name='identifications.tsv')
    session.store(
        file_manager, dataset_id, levels, runs,
        design_of(workflow.workflow_dir, selected_files, mzmls),
        {'runs': list(dataset_ids), 'tool': TOOL}, tool,
    )
    rmtree(folder, ignore_errors=True)
    return dataset_id


def compare_isobaric(workflow, tool, dataset_ids, plex):
    """Store the comparison of the labelled samples of the runs `dataset_ids`
    (fractions of one labelled set): the reporter ions of every spectrum
    (stored by the workflow as `isobaric_dfs`), with the PrSMs of an
    identified run. Returns the comparison's dataset id, or None when no
    spectrum carried reporter ions."""
    file_manager, logger = workflow.file_manager, workflow.logger
    channels = isobaric.load_channels(workflow.workflow_dir, plex)
    runs = []
    for dataset_id in dataset_ids:
        quant = file_manager.get_results(dataset_id, ['isobaric_dfs'], use_polars=True)['isobaric_dfs'].collect()
        prsm_path = file_manager.get_results(dataset_id, ['prsms_tsv'], partial=True).get('prsms_tsv')
        runs.append((dataset_id, quant, isobaric.read_prsms(prsm_path) if prsm_path else None))
    logger.log(f'-> Comparing the {len(channels)} labelled samples ({isobaric.label(plex)})...')
    levels = isobaric.build(runs, channels, database_of(workflow))
    if not levels:
        logger.log('No spectrum carried reporter ions with a precursor mass: no comparison made.')
        return None
    dataset_id = f"comparison_{time.strftime('%Y%m%d-%H%M%S')}"
    samples = list(channels['sample'])
    session.store(
        file_manager, dataset_id, levels, samples, isobaric.design(channels),
        {'runs': list(dataset_ids), 'tool': 'FLASHDeconv', 'plex': plex}, tool,
        extra={'kind': 'isobaric', 'plex': plex, 'fractions': list(dataset_ids),
               'channels': list(channels['channel'])},
    )
    return dataset_id
