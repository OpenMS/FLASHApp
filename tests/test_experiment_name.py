"""The name typed beside Run analysis names the runs a start processed and its
comparison, as the Viewer lists them."""

import os
import sys
from pathlib import Path
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.workflow.FileManager import FileManager  # noqa: E402


def _workflow(tmp_path, name):
    # Imported here: see test_merged_workflow.py
    from src.Workflow import TopDownWorkflow
    from src.workflow.RunManager import RUN_NAME_FILE

    wf = object.__new__(TopDownWorkflow)
    wf.workflow_dir = Path(tmp_path, 'flashtnt')
    wf.workflow_dir.mkdir()
    if name:
        Path(wf.workflow_dir, RUN_NAME_FILE).write_text(name)
    cache = Path(tmp_path, 'cache')
    cache.mkdir(exist_ok=True)
    wf.file_manager = FileManager(cache, cache)
    wf.logger = mock.MagicMock()
    return wf


def _store(wf, *dataset_ids):
    for dataset_id in dataset_ids:
        wf.file_manager.store_data(dataset_id, 'settings', {'x': 1})


def test_one_run_takes_the_name(tmp_path):
    wf = _workflow(tmp_path, 'Liver, day 3')
    _store(wf, 'liver_20261108-120000')
    wf.name_results(['liver_20261108-120000'])
    assert wf.file_manager.get_display_name('liver_20261108-120000') == 'Liver, day 3'


def test_a_comparison_takes_the_name_and_its_runs_their_files(tmp_path):
    wf = _workflow(tmp_path, 'Liver')
    runs = ['a_20261108-120000', 'b_20261108-120001']
    _store(wf, *runs, 'comparison_20261108-120100')
    wf.name_results(runs, 'comparison_20261108-120100')
    assert wf.file_manager.get_display_name('comparison_20261108-120100') == 'Liver'
    assert [wf.file_manager.get_display_name(r) for r in runs] == ['Liver: a', 'Liver: b']


def test_a_name_given_twice_stays_apart(tmp_path):
    wf = _workflow(tmp_path, 'Liver')
    _store(wf, 'a_20261108-120000', 'a_20261108-130000')
    wf.name_results(['a_20261108-120000'])
    wf.name_results(['a_20261108-130000'])
    assert wf.file_manager.get_display_name('a_20261108-130000') == 'Liver (2)'


def test_without_a_name_the_runs_keep_theirs(tmp_path):
    wf = _workflow(tmp_path, '')
    _store(wf, 'a_20261108-120000')
    wf.name_results(['a_20261108-120000'])
    assert wf.file_manager.get_display_name('a_20261108-120000') == 'a_20261108-120000'
