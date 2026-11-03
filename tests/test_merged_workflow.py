"""The merged FLASHDeconv & FLASHTnT workflow: what carries over from the two
former workflows of a workspace (their uploads and the comparison layout)."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pathlib import Path  # noqa: E402

from src.Workflow import (  # noqa: E402
    LEGACY_DECONV_INPUTS, DeconvWorkflow, TagWorkflow, TopDownWorkflow, adopt_legacy_inputs,
)
from src.workflow.FileManager import FileManager  # noqa: E402


def _external(workflow_dir):
    listed = Path(workflow_dir, 'input-files', 'mzML-files', 'external_files.txt')
    return listed.read_text().splitlines() if listed.exists() else []


def test_former_workflows_name_the_merged_one():
    assert TagWorkflow is TopDownWorkflow and DeconvWorkflow is TopDownWorkflow


def test_deconvolution_uploads_are_listed_once_by_reference(tmp_path):
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
