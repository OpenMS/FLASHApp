"""The merged FLASHDeconv & FLASHTnT workflow: what carries over from the two
former workflows of a workspace (their uploads and the comparison layout)."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pathlib import Path  # noqa: E402

from src.workflow.FileManager import FileManager  # noqa: E402


def _external(workflow_dir):
    listed = Path(workflow_dir, 'input-files', 'mzML-files', 'external_files.txt')
    return listed.read_text().splitlines() if listed.exists() else []


# src.Workflow is imported inside the tests: importing it here, at collection,
# would load src.workflow.ParameterManager with the real streamlit before
# test_parameter_defaults.py swaps in its mock.


def test_former_workflows_name_the_merged_one():
    from src.Workflow import DeconvWorkflow, TagWorkflow, TopDownWorkflow
    assert TagWorkflow is TopDownWorkflow and DeconvWorkflow is TopDownWorkflow


def test_deconvolution_uploads_are_listed_once_by_reference(tmp_path):
    from src.Workflow import LEGACY_DECONV_INPUTS, adopt_legacy_inputs
    legacy = Path(tmp_path, LEGACY_DECONV_INPUTS)
    legacy.mkdir(parents=True)
    Path(legacy, 'run_a.mzML').write_text('a')
    Path(legacy, 'shared.mzML').write_text('b')
    Path(legacy, 'external_files.txt').write_text('/data/run_b.raw\n')
    workflow_dir = Path(tmp_path, 'flashtnt')
    target = Path(workflow_dir, 'input-files', 'mzML-files')
    target.mkdir(parents=True)
    Path(target, 'shared.mzML').write_text('b')  # already among the identification runs

    adopted = adopt_legacy_inputs(tmp_path, workflow_dir)

    assert sorted(p.name for p in adopted) == ['run_a.mzML', 'run_b.raw']
    assert sorted(Path(line).name for line in _external(workflow_dir)) == ['run_a.mzML', 'run_b.raw']
    # Referenced, not copied
    assert sorted(p.name for p in target.iterdir()) == ['external_files.txt', 'shared.mzML']
    # A file removed later from the merged workflow does not come back
    Path(target, 'external_files.txt').write_text('')
    assert adopt_legacy_inputs(tmp_path, workflow_dir) == []
    assert _external(workflow_dir) == []


def test_workspace_without_deconvolution_uploads(tmp_path):
    from src.Workflow import adopt_legacy_inputs
    assert adopt_legacy_inputs(tmp_path, Path(tmp_path, 'flashtnt')) == []
    assert _external(Path(tmp_path, 'flashtnt')) == []


def _layout(rows):
    return {'layout': [rows], 'side_by_side': False}


def test_comparison_layout_of_identified_runs_carries_over(tmp_path):
    from src.render import quant_view

    file_manager = FileManager(tmp_path, Path(tmp_path, 'cache'))
    file_manager.store_data('flashtnt_quant_layout', 'layout', _layout([['quant_table']]))
    quant_view.migrate_layout(file_manager)
    assert quant_view.layout(file_manager) == [['quant_table']]
    assert not file_manager.result_exists('flashtnt_quant_layout', 'layout')


def test_saved_comparison_layout_wins(tmp_path):
    from src.render import quant_view

    file_manager = FileManager(tmp_path, Path(tmp_path, 'cache'))
    file_manager.store_data('quant_layout', 'layout', _layout([['quant_volcano']]))
    file_manager.store_data('flashtnt_quant_layout', 'layout', _layout([['quant_table']]))
    quant_view.migrate_layout(file_manager)
    assert quant_view.layout(file_manager) == [['quant_volcano']]
    assert not file_manager.result_exists('flashtnt_quant_layout', 'layout')


def test_one_catalog_holds_both_kinds():
    from src.render.layout_editor.catalog import (
        DECONV_DEFAULT_LAYOUT, TNT_DEFAULT_LAYOUTS, TOPDOWN_COMPONENTS, TOPDOWN_SEQUENCE_COMPONENTS,
        requirements, validate_layout,
    )
    catalog = TOPDOWN_COMPONENTS + TOPDOWN_SEQUENCE_COMPONENTS
    names = [c['name'] for c in catalog]
    assert len(names) == len(set(names))
    assert all(r in names for c in catalog for r in requirements(c))
    for layout in [DECONV_DEFAULT_LAYOUT, DECONV_DEFAULT_LAYOUT + [['sequence_view']], *TNT_DEFAULT_LAYOUTS.values()]:
        assert validate_layout([layout], catalog) == ''
    # Spectra follow an identification as well as the scan table
    assert validate_layout([[['protein_table'], ['deconv_spectrum', 'mass_table']]], catalog) == ''
    assert 'needs' in validate_layout([[['deconv_spectrum']]], catalog)


def test_each_experiment_starts_from_the_default_of_its_run(tmp_path):
    from src.render import topdown_layout
    from src.render.layout_editor.catalog import DECONV_DEFAULT_LAYOUT, TNT_DEFAULT_LAYOUTS

    file_manager = FileManager(tmp_path, Path(tmp_path, 'cache'))
    # A run deconvolved only ends with the sequence view, which asks for a sequence
    deconvolution = DECONV_DEFAULT_LAYOUT + [['sequence_view']]
    assert topdown_layout.rows(file_manager, False, None) == deconvolution
    assert topdown_layout.rows(file_manager, True, 'Native') == TNT_DEFAULT_LAYOUTS['Native']
    assert topdown_layout.rows(file_manager, True, 'Lysate') == TNT_DEFAULT_LAYOUTS['Lysate']
    # A layout the former identification Viewer saved stays the default of identified runs
    file_manager.store_data('flashtnt_layout', 'layout', _layout([['protein_table']]))
    assert topdown_layout.rows(file_manager, True, 'Lysate') == [['protein_table']]
    assert topdown_layout.rows(file_manager, False, None) == deconvolution
    # A saved layout is shown for every run
    file_manager.store_data('topdown_layout', 'layout', _layout([['scan_table']]))
    assert topdown_layout.rows(file_manager, True, 'Native') == [['scan_table']]
    assert topdown_layout.rows(file_manager, False, None) == [['scan_table']]


def test_sequence_input_saves_and_reads_back(tmp_path):
    from src.render import sequence_input

    file_manager = FileManager(tmp_path, Path(tmp_path, 'cache'))
    assert sequence_input.saved(file_manager) is None
    assert sequence_input.is_valid(' mkv laa\n') and sequence_input.is_valid('ACXD')
    assert not sequence_input.is_valid('ACBD') and not sequence_input.is_valid('  ')
    sequence_input.save(file_manager, ' mkv laa\n', ['Carbamidomethyl (C)'])
    assert sequence_input.saved(file_manager) == ('MKVLAA', ['Carbamidomethyl (C)'])
    sequence_input.clear(file_manager)
    assert sequence_input.saved(file_manager) is None
