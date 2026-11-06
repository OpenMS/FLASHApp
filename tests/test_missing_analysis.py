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


def test_options_rerun_add_files_and_bring_in_a_run_that_has_it(tmp_path):
    file_manager = _file_manager(tmp_path)
    runs = ['example_deconvolution', 'example_native_identified', 'example_isobaric']
    labels = [o['label'] for o in missing_analysis.options(
        tmp_path, file_manager, 'example_deconvolution', [IDENTIFICATION], {'name': 'native'},
        runs=runs, selection_key=lambda i: f'sel{i}', view_args=None)]
    assert labels == [
        'Run Example: deconvolution only (protein mix, LC-MS) again with identification',
        'Process more files with these settings',
        'Show Example: native, identified (GapA)',
        'Compare side by side with Example: native, identified (GapA)',
    ]
    # In the second experiment of a comparison there is nothing to put beside it
    second = missing_analysis.options(
        tmp_path, file_manager, 'example_deconvolution', [IDENTIFICATION], None,
        runs=runs, selection_key=lambda i: f'sel{i}', view_args=None, index=1)
    assert not any(o['label'].startswith('Compare') for o in second)
