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


def test_options_only_process_more_files(tmp_path):
    file_manager = _file_manager(tmp_path)
    labels = [o['label'] for o in missing_analysis.options(
        tmp_path, file_manager, 'example_deconvolution', [IDENTIFICATION])]
    assert labels == ['Process more files with these settings']


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
