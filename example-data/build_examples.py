"""Build the example results of the demo workspace.

Each example shows one default layout of the Viewer, so a new workspace opens
with one run per kind of analysis, named after it. The script runs the
Top-down analysis on the example inputs exactly as the Run page does
(tasks.execute_workflow) and stores the results in the demo workspace's
shared cache (example-data/workspaces/default/cache), under a fixed dataset id
and a display name that says what the example shows.

It needs the TOPP tools (FLASHDeconv, DecoyDatabase, FLASHTnT, FLASHAnnotator)
on PATH, so run it inside the FLASHApp image from the repository root:

    docker run --rm -v "$PWD:/repo" -w /repo ghcr.io/openms/flashapp:latest \
        mamba run -n streamlit-env python example-data/build_examples.py
"""

import json
import shutil
import sqlite3
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

DEMO = REPO / 'example-data' / 'workspaces' / 'default'
INPUTS = DEMO / 'flashtnt' / 'input-files'
DATABASE = INPUTS / 'fasta-file' / 'example_database.fasta'

# dataset id: (input file, display name, analysis parameters)
EXAMPLES = {
    'example_deconvolution': (
        'Deconvolution-only_ProteinMix_LC-MS.mzML',
        'Example: deconvolution only (protein mix, LC-MS)',
        {'identify': False, 'quantification': 'none'},
    ),
    'example_native_identified': (
        'Native-identified_GapA_HCD.mzML',
        'Example: native, identified (GapA)',
        # The Native preset (presets.json)
        {'identify': True, 'quantification': 'none', 'FLASHTnT': {'native': True}, 'few_proteins': True},
    ),
}


def run(dataset_id, file_name, params, workspace):
    """Run one example in `workspace` and return the dataset id it was stored under."""
    from src.workflow.tasks import execute_workflow
    from src.workflow.FileManager import FileManager

    workflow_dir = workspace / 'flashtnt'
    shutil.rmtree(workflow_dir, ignore_errors=True)
    workflow_dir.mkdir(parents=True)
    params = dict(params, **{'mzML-files': [str(INPUTS / 'mzML-files' / file_name)]})
    if params.get('identify'):
        params['fasta-file'] = str(DATABASE)
    (workflow_dir / 'params.json').write_text(json.dumps(params, indent=4))

    cache = workspace / 'cache'
    before = set(FileManager(workflow_dir, cache).get_results_list(['threedim_SN_plot']))
    execute_workflow(str(workflow_dir), 'TopDownWorkflow', 'src.Workflow')
    stored = set(FileManager(workflow_dir, cache).get_results_list(['threedim_SN_plot'])) - before
    if len(stored) != 1:
        print((workflow_dir / 'logs' / 'all.log').read_text())
        raise SystemExit(f'{dataset_id}: expected one new result, found {sorted(stored)}')
    return stored.pop()


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
    with tempfile.TemporaryDirectory() as tmp:
        workspace = Path(tmp, 'workspace')
        for dataset_id, (file_name, display_name, params) in EXAMPLES.items():
            print(f'== {display_name}', flush=True)
            stored = run(dataset_id, file_name, params, workspace)
            rename(workspace / 'cache', stored, dataset_id, display_name)
        target = DEMO / 'cache'
        shutil.rmtree(target, ignore_errors=True)
        target.mkdir()
        shutil.copy(workspace / 'cache' / 'cache.db', target / 'cache.db')
        for dataset_id in EXAMPLES:
            shutil.copytree(workspace / 'cache' / 'files' / dataset_id, target / 'files' / dataset_id)
    print(f'Stored {len(EXAMPLES)} examples in {target.relative_to(REPO)}')


if __name__ == '__main__':
    main()
