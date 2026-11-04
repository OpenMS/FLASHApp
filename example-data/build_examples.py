"""Build the example results of the demo workspace.

Each example shows one default layout of the Viewer, so a new workspace opens
with one run per kind of analysis, named after it. The script runs the
Top-down analysis on the example inputs exactly as the Run page does
(tasks.execute_workflow) and stores the results in the demo workspace's
shared cache (example-data/workspaces/default/cache), under a fixed dataset id
and a display name that says what the example shows.

The label-free example uses four runs of the FLASHQuant SpikeIn dataset
(MassIVE MSV000091923 / PXD042169: the six-protein mix spiked into an E. coli
lysate at relative amounts 1 (A) and 1/10 (F), two replicates each). A missing
input is downloaded from MassIVE and cut to a short retention-time window with
FileFilter, so the repository stays small.

It needs the TOPP tools (FileFilter, FLASHDeconv, DecoyDatabase, FLASHTnT,
FLASHAnnotator, FLASHQuantConsensus) on PATH, so run it inside the FLASHApp image from the repository root:

    docker run --rm -v "$PWD:/repo" -w /repo ghcr.io/openms/flashapp:latest \
        mamba run -n streamlit-env python example-data/build_examples.py
"""

import json
import shutil
import subprocess
import urllib.error
import urllib.request
import sqlite3
import sys
import tempfile
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from src.Workflow import FLASHANNOTATOR_FLAGS, FLASHDECONV_FLAGS, FLASHTNT_FLAGS  # noqa: E402

DEMO = REPO / 'example-data' / 'workspaces' / 'default'
INPUTS = DEMO / 'flashtnt' / 'input-files'
DATABASE = INPUTS / 'fasta-file' / 'example_database.fasta'

MASSIVE = 'https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?forceDownload=true&file=f.MSV000091923/peak/'
# Retention time window kept of the SpikeIn runs (seconds): spike-in proteins
# and E. coli background elute here
SPIKEIN_RT = (3990, 4140)
# Example input -> run of MSV000091923 it is cut from
SOURCES = {
    'Label-free_SpikeIn-1x_R1.mzML': '20220414_Ecoli_ProtMix_A_R1.mzML',
    'Label-free_SpikeIn-1x_R2.mzML': '20220414_Ecoli_ProtMix_A_R2.mzML',
    'Label-free_SpikeIn-0.1x_R1.mzML': '20220414_Ecoli_ProtMix_F_R1.mzML',
    'Label-free_SpikeIn-0.1x_R2.mzML': '20220414_Ecoli_ProtMix_F_R2.mzML',
}

# Each example: the input files, its analysis parameters, and the dataset id
# and display name of every result it stores (one per run, then the
# comparison of a label-free example).
EXAMPLES = [
    {
        'files': ['Deconvolution-only_ProteinMix_LC-MS.mzML'],
        'params': {'identify': False, 'quantification': 'none'},
        'runs': [('example_deconvolution', 'Example: deconvolution only (protein mix, LC-MS)')],
    },
    {
        'files': ['Native-identified_GapA_HCD.mzML'],
        # The Native preset (presets.json)
        'params': {'identify': True, 'quantification': 'none', 'FLASHTnT': {'native': True}, 'few_proteins': True},
        'runs': [('example_native_identified', 'Example: native, identified (GapA)')],
    },
    {
        'files': list(SOURCES),
        # The Lysate preset (presets.json)
        'params': {'identify': True, 'quantification': 'label-free', 'FLASHTnT': {'native': False},
                   'few_proteins': False},
        'conditions': ['Spike-in 1x', 'Spike-in 1x', 'Spike-in 0.1x', 'Spike-in 0.1x'],
        'replicates': [1, 2, 1, 2],
        # The first run is the lysate example. The comparison keeps its own
        # tables, so the other runs' results are not shipped (None).
        'runs': [('example_lysate_identified', 'Example: lysate, identified (E. coli + spike-in)'), None, None, None],
        'comparison': ('example_label_free', 'Example: label-free comparison (spike-in 1x vs 0.1x)'),
    },
]

# Valueless switches, as the Configure page records them (_flag_params)
FLAGS = {
    'FLASHDeconv': FLASHDECONV_FLAGS, 'FLASHTnT': FLASHTNT_FLAGS, 'FLASHAnnotator': FLASHANNOTATOR_FLAGS,
}


def download(url, path, attempts=6):
    """Download `url` to `path`, waiting longer after each refusal (MassIVE rate-limits)."""
    request = urllib.request.Request(url, headers={'User-Agent': 'FLASHApp-examples/1.0'})
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(request, timeout=600) as response, open(path, 'wb') as out:
                shutil.copyfileobj(response, out)
            return
        except (urllib.error.URLError, TimeoutError) as error:
            if attempt == attempts - 1:
                raise
            wait = 60 * 2 ** attempt
            print(f'{error}: retrying in {wait} s', flush=True)
            time.sleep(wait)


def trim(source, target):
    """Keep the spectra within SPIKEIN_RT (FileFilter); `source` may be `target`."""
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp, 'trimmed.mzML')
        subprocess.run(
            ['FileFilter', '-in', str(source), '-out', str(out), '-rt', f'{SPIKEIN_RT[0]}:{SPIKEIN_RT[1]}'],
            check=True,
        )
        shutil.move(out, target)


def fetch(name):
    """Download the SpikeIn run behind example input `name` and keep SPIKEIN_RT."""
    target = INPUTS / 'mzML-files' / name
    if target.exists():
        trim(target, target)
        return
    with tempfile.TemporaryDirectory() as tmp:
        full = Path(tmp, SOURCES[name])
        print(f'Downloading {SOURCES[name]}', flush=True)
        download(MASSIVE + SOURCES[name], full)
        trim(full, target)
    print(f'{name}: {target.stat().st_size / 1e6:.1f} MB', flush=True)


def run(example, workspace):
    """Run one example in `workspace`; return the stored run ids (in the order
    of its files) and the comparison id (None without one)."""
    from src.workflow import ExperimentalDesign
    from src.workflow.FileManager import FileManager
    from src.workflow.tasks import execute_workflow

    workflow_dir = workspace / 'flashtnt'
    shutil.rmtree(workflow_dir, ignore_errors=True)
    workflow_dir.mkdir(parents=True)
    files = [str(INPUTS / 'mzML-files' / name) for name in example['files']]
    params = dict(example['params'], **{'mzML-files': files, '_flag_params': FLAGS})
    if params.get('identify'):
        params['fasta-file'] = str(example.get('database', DATABASE))
    (workflow_dir / 'params.json').write_text(json.dumps(params, indent=4))
    if 'conditions' in example:
        design = ExperimentalDesign.default_design(files)
        design['condition'] = example['conditions']
        design['replicate'] = example['replicates']
        ExperimentalDesign.save_design(workflow_dir, design, files)

    cache = workspace / 'cache'

    def stored():
        file_manager = FileManager(workflow_dir, cache)
        return (set(file_manager.get_results_list(['threedim_SN_plot'])),
                set(file_manager.get_results_list(['quant_manifest'])))

    runs_before, comparisons_before = stored()
    try:
        execute_workflow(str(workflow_dir), 'TopDownWorkflow', 'src.Workflow')
    except Exception:
        print((workflow_dir / 'logs' / 'all.log').read_text())
        raise
    runs_after, comparisons_after = stored()
    new_runs, new_comparisons = runs_after - runs_before, comparisons_after - comparisons_before
    run_ids = []
    for name in example['files']:
        stem = Path(name).stem
        found = [d for d in new_runs if d.rsplit('_', 1)[0] == stem]
        if len(found) != 1:
            print((workflow_dir / 'logs' / 'all.log').read_text())
            raise SystemExit(f'{name}: expected one new result, found {sorted(new_runs)}')
        run_ids.append(found[0])
    expected = 1 if 'comparison' in example else 0
    if len(new_comparisons) != expected:
        print((workflow_dir / 'logs' / 'all.log').read_text())
        raise SystemExit(f'expected {expected} new comparison, found {sorted(new_comparisons)}')
    return run_ids, (new_comparisons.pop() if new_comparisons else None)


def rename(cache, old_id, new_id, display_name):
    """Store result `old_id` as `new_id` (folder and index) with a display name."""
    folder = Path(cache, 'files', old_id)
    shutil.move(folder, Path(cache, 'files', new_id))
    db = sqlite3.connect(Path(cache, 'cache.db'))
    for table in ('stored_data', 'stored_files'):
        columns = [row[1] for row in db.execute(f'PRAGMA table_info({table})')]
        if 'display_name' not in columns:
            db.execute(f'ALTER TABLE {table} ADD COLUMN display_name TEXT')
            columns.append('display_name')
        row = db.execute(f'SELECT * FROM {table} WHERE id = ?', (old_id,)).fetchone()
        if row is None:
            continue
        values = dict(zip(columns, row))
        values['id'] = new_id
        values['display_name'] = display_name
        for key, value in values.items():
            if isinstance(value, str) and value.startswith(f'files/{old_id}/'):
                values[key] = f'files/{new_id}/' + value[len(f'files/{old_id}/'):]
        db.execute(f'DELETE FROM {table} WHERE id = ?', (old_id,))
        db.execute(
            f'INSERT INTO {table} ({", ".join(values)}) VALUES ({", ".join("?" * len(values))})',
            list(values.values()),
        )
    db.commit()
    db.close()


def main():
    from src.workflow.FileManager import FileManager

    kept = []
    with tempfile.TemporaryDirectory() as tmp:
        workspace = Path(tmp, 'workspace')
        cache = workspace / 'cache'
        for example in EXAMPLES:
            for name in example['files']:
                if name in SOURCES:
                    fetch(name)
            named = [entry for entry in example['runs'] if entry] + [example.get('comparison')]
            print(f"== {', '.join(name for name in (e[1] for e in named if e))}", flush=True)
            run_ids, comparison = run(example, workspace)
            for old_id, entry in zip(run_ids, example['runs']):
                if entry is None:
                    FileManager(workspace / 'flashtnt', cache).remove_results(old_id)
                    continue
                rename(cache, old_id, *entry)
                kept.append(entry[0])
            if comparison is not None:
                rename(cache, comparison, *example['comparison'])
                kept.append(example['comparison'][0])
        target = DEMO / 'cache'
        shutil.rmtree(target, ignore_errors=True)
        target.mkdir()
        shutil.copy(cache / 'cache.db', target / 'cache.db')
        for dataset_id in kept:
            shutil.copytree(cache / 'files' / dataset_id, target / 'files' / dataset_id)
    print(f'Stored {len(kept)} example results in {target.relative_to(REPO)}')


if __name__ == '__main__':
    main()
