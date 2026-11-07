"""Views a run was not processed for: what the Viewer locks and what it offers instead."""

import importlib.util
import json
from pathlib import Path

from src.render import missing_analysis
from src.render.layout_editor.catalog import IDENTIFICATION, ISOBARIC, TOPDOWN_COMPONENTS
from src.workflow.FileManager import FileManager

DEMO = Path('example-data', 'workspaces', 'default')


def _file_manager(tmp_path):
    return FileManager(tmp_path, DEMO / 'cache')


def test_components_name_the_analysis_they_need():
    needs = {c['name']: c.get('needs') for c in TOPDOWN_COMPONENTS}
    assert needs['protein_table'] == IDENTIFICATION
    assert needs['reporter_table'] == ISOBARIC
    # Deconvolution components and the sequence view (a sequence can be set by hand) need nothing
    assert needs['ms1_deconv_heat_map'] is None
    assert needs['sequence_view'] is None


def test_available_reads_what_each_run_was_processed_with(tmp_path):
    file_manager = _file_manager(tmp_path)
    assert missing_analysis.available(file_manager, 'example_deconvolution') == set()
    assert missing_analysis.available(file_manager, 'example_native_identified') == {IDENTIFICATION}
    assert missing_analysis.available(file_manager, 'example_isobaric') == {ISOBARIC}
    assert missing_analysis.available(file_manager, None) is None


def test_example_runs_name_the_input_files_they_were_made_from():
    spec = importlib.util.spec_from_file_location('build_examples', Path('example-data', 'build_examples.py'))
    build = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(build)
    made_from = {run[0]: Path(example['files'][i]).stem
                 for example in build.EXAMPLES for i, run in enumerate(example['runs']) if run}
    assert missing_analysis.EXAMPLE_INPUTS == made_from


def test_input_file_of_a_run_and_of_an_example(tmp_path):
    folder = Path(tmp_path, 'flashtnt', 'input-files', 'mzML-files')
    folder.mkdir(parents=True)
    Path(folder, 'my_sample_1.mzML').touch()
    Path(folder, 'Deconvolution-only_ProteinMix_LC-MS.mzML').touch()
    elsewhere = Path(tmp_path, 'elsewhere', 'other.raw')
    Path(folder, 'external_files.txt').write_text(f'{elsewhere}\n')

    # The time stamp goes, underscores of the file name stay
    assert missing_analysis.input_file(tmp_path, 'my_sample_1_20260101-120000').endswith('my_sample_1.mzML')
    assert missing_analysis.input_file(tmp_path, 'other_20260101-120000') == str(elsewhere)
    assert missing_analysis.input_file(tmp_path, 'example_deconvolution').endswith(
        'Deconvolution-only_ProteinMix_LC-MS.mzML')
    assert missing_analysis.input_file(tmp_path, 'gone_20260101-120000') is None


def test_run_settings_are_the_ones_the_run_was_made_with(tmp_path):
    file_manager = _file_manager(tmp_path)
    Path(tmp_path, 'flashtnt').mkdir()
    Path(tmp_path, 'flashtnt', 'params.json').write_text(json.dumps(
        {'fasta-file': 'db.fasta', 'identify': True, 'mzML-files': ['a.mzML']}))

    params, has_file = missing_analysis.run_settings(tmp_path, file_manager, 'example_deconvolution')

    assert params['fasta-file'] == 'db.fasta'  # the workflow's current settings stay
    assert params['identify'] is False  # the run was deconvolved only
    assert params['FLASHDeconv'] == file_manager.get_results(
        'example_deconvolution', ['FD_parameters'])['FD_parameters']
    assert not has_file  # its input file is not in this workspace


def test_options_add_identification_to_the_run_and_process_more_files_for_labels(tmp_path):
    file_manager = _file_manager(tmp_path)
    name = file_manager.get_display_name('example_deconvolution')
    labels = [o['label'] for o in missing_analysis.options(
        tmp_path, file_manager, 'example_deconvolution', [IDENTIFICATION])]
    assert labels == [f'Add identification to {name}']
    labels = [o['label'] for o in missing_analysis.options(
        tmp_path, file_manager, 'example_deconvolution', [ISOBARIC])]
    assert labels == ['Process more files with these settings']


def _capture_start(monkeypatch):
    started = {}

    class Workflow:
        def start_from(self, params, changes=None, **kwargs):
            started.update(params=params, changes=changes, **kwargs)

    import src.Workflow
    monkeypatch.setattr(src.Workflow, 'TopDownWorkflow', Workflow)
    return started


def test_add_identification_locks_the_runs_file_and_reuses_its_deconvolution(tmp_path, monkeypatch):
    from src import analysis

    file_manager = _file_manager(tmp_path)
    folder = Path(tmp_path, 'flashtnt', 'input-files', 'mzML-files')
    folder.mkdir(parents=True)
    source = Path(folder, 'Deconvolution-only_ProteinMix_LC-MS.mzML')
    source.touch()
    Path(folder, 'other.mzML').touch()
    Path(tmp_path, 'flashtnt', 'params.json').write_text(json.dumps(
        {'fasta-file': 'db.fasta', 'mzML-files': [str(Path(folder, 'other.mzML'))],
         analysis.QUANT_KEY: analysis.LABEL_FREE}))
    started = _capture_start(monkeypatch)

    missing_analysis.add_identification(tmp_path, file_manager, 'example_deconvolution')

    assert started['same_files'] == 'mzML-files' and started['files'] == [str(source)]
    assert started['keep'] == ['FLASHDeconv'] and started['reuses'] == ['example_deconvolution']
    # Identification on, the quantification as the run was processed (none)
    assert started['changes'][analysis.IDENTIFY_KEY] is True
    assert started['changes'][analysis.QUANT_KEY] == analysis.NO_QUANT
    assert started['params']['fasta-file'] == 'db.fasta'


def test_add_identification_without_the_runs_file_picks_files_again(tmp_path, monkeypatch):
    started = _capture_start(monkeypatch)
    missing_analysis.add_identification(tmp_path, _file_manager(tmp_path), 'example_deconvolution')
    assert 'same_files' not in started and 'keep' not in started


def test_execution_reuses_the_runs_deconvolution(tmp_path):
    import shutil
    from src.Workflow import TopDownWorkflow

    cache = Path(tmp_path, 'cache')
    shutil.copytree(DEMO / 'cache', cache)
    workflow = object.__new__(TopDownWorkflow)  # as tasks.execute_workflow builds it
    workflow.workflow_dir = Path(tmp_path, 'flashtnt')
    workflow.workflow_dir.mkdir()
    workflow.file_manager = FileManager(tmp_path, cache)
    assert workflow._reused_deconvolution(1) == [None]

    Path(workflow.workflow_dir, 'last-follow-up.json').write_text(json.dumps(
        {'name': 'x', 'same_files': 'mzML-files', 'files': ['a.mzML'], 'kept': ['FLASHDeconv'],
         'reuses': ['example_deconvolution']}))
    assert workflow._reused_deconvolution(1) == ['example_deconvolution']
    assert workflow._reused_deconvolution(2) == [None, None]  # files and runs do not pair up

    out = Path(tmp_path, 'out')
    out.mkdir()
    assert workflow._copy_deconvolution('example_deconvolution', out)
    assert {'out_deconv.mzML', 'anno_annotated.mzML', 'spec1.tsv', 'out.tsv'} <= {p.name for p in out.iterdir()}
    assert not workflow._copy_deconvolution('unknown_run', Path(tmp_path, 'out2'))


def test_more_files_compare_the_added_files_with_the_run(tmp_path):
    from src import analysis
    from src.quant.run import matches, reusable_runs, run_record

    mzml = Path(tmp_path, 'a.mzML')
    mzml.write_text('')
    file_manager = FileManager(tmp_path, Path(tmp_path, 'cache'))
    record = run_record(True, analysis.LABEL_FREE, None, mzml, ('control', 1))
    file_manager.store_data('a_20260101-000000', 'analysis_settings', record)

    # A run processed with what the start needs joins the comparison as it is
    assert reusable_runs(file_manager, ['a_20260101-000000', 'unknown'], True, analysis.LABEL_FREE, None) \
        == [('a_20260101-000000', record)]
    # Not when the start needs what the run was not processed with
    assert not matches(record, True, analysis.ISOBARIC, 'tmt6plex')
    assert matches(run_record(False, analysis.LABEL_FREE, None, mzml), False, analysis.LABEL_FREE, None)
    assert not matches(run_record(False, analysis.LABEL_FREE, None, mzml), True, analysis.LABEL_FREE, None)

    # Processing more files quantifies: the files are compared with the run
    params, changes, adds = missing_analysis.more_files_settings(
        {'mzML-files': ['a.mzML'], analysis.IDENTIFY_KEY: True}, {}, lambda *a: True)
    assert adds and params['mzML-files'] == []
    assert changes == {analysis.QUANT_KEY: analysis.LABEL_FREE}
    # A run without the analysis switched on now is processed again with the new files
    params, changes, adds = missing_analysis.more_files_settings(
        {'mzML-files': ['a.mzML'], analysis.IDENTIFY_KEY: False}, {analysis.IDENTIFY_KEY: True}, lambda *a: False)
    assert not adds and params['mzML-files'] == ['a.mzML']


def test_the_viewer_finds_the_comparison_of_added_files(tmp_path, monkeypatch):
    from src.render import quant_view

    workflow_dir = Path(tmp_path, 'flashtnt')
    workflow_dir.mkdir()
    file_manager = type('FM', (), {'get_results': lambda self, d, k: {'quant_manifest': {
        'comparison_2': {'runs': ['a', 'b', 'c']}, 'comparison_1': {'runs': ['b', 'c']}}[d]}})()
    monkeypatch.setattr(quant_view, 'comparisons', lambda fm, tool=None: ['comparison_2', 'comparison_1'])
    assert quant_view.added_comparison(file_manager, workflow_dir) is None
    Path(workflow_dir, 'last-follow-up.json').write_text(json.dumps({'adds_to': ['a']}))
    assert quant_view.added_comparison(file_manager, workflow_dir) == 'comparison_2'
