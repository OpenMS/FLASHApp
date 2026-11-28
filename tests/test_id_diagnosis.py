"""Identification diagnostics (src/render/id_diagnosis.py).

Hand-made step counts, plus the identified lysate example of the demo
workspace shipped in example-data/.
"""
import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.render import id_diagnosis  # noqa: E402
from src.workflow.FileManager import FileManager  # noqa: E402

WORKSPACE = Path(__file__).resolve().parent.parent / 'example-data' / 'workspaces' / 'default'
LYSATE = 'example_lysate_identified'


def _counts(**changes):
    counts = {
        'ms1': 40, 'ms2': 1000, 'deconvolved': 950, 'precursor': 900, 'tagged': 800,
        'matched': 300, 'kept': 250, 'q_estimated': True, 'calibration_ppm': 1.0,
        'tolerance_ppm': 10.0,
    }
    counts.update(changes)
    return counts


def _causes(table):
    return dict(zip(table['step'], table['cause']))


def test_a_good_run_names_no_cause():
    table = id_diagnosis.diagnose(_counts())
    assert table['step'].to_list() == list(id_diagnosis.STEP_LABELS.values())
    assert table['spectra'].to_list() == [1000, 950, 900, 800, 300, 250]
    assert table['share'].to_list() == pytest.approx([100.0, 95.0, 90.0, 80.0, 30.0, 25.0])
    assert not any(table['cause'])
    assert 'No step loses' in id_diagnosis.headline(_counts())


def test_no_ms2_spectra():
    counts = _counts(ms2=0, deconvolved=0, precursor=0, tagged=0, matched=0, kept=0)
    causes = _causes(id_diagnosis.diagnose(counts))
    assert 'no MS2 spectra' in causes['MS2 spectra']
    assert sum(bool(c) for c in causes.values()) == 1
    assert id_diagnosis.headline(counts) == 'This run has no MS2 spectra to identify.'


def test_failed_deconvolution_and_few_tags():
    causes = _causes(id_diagnosis.diagnose(_counts(deconvolved=300, tagged=100, matched=50, kept=40)))
    assert 'only 30 %' in causes['with deconvolved masses']
    assert 'tag' in causes['with sequence tags']


def test_database_that_does_not_fit():
    causes = _causes(id_diagnosis.diagnose(_counts(matched=20, kept=15)))
    assert 'database' in causes['matched to a proteoform']


def test_no_ms1_spectra_explains_the_missing_precursor_masses():
    causes = _causes(id_diagnosis.diagnose(_counts(ms1=0, precursor=0)))
    assert 'no MS1 spectra' in causes['with a precursor mass']


def test_no_q_values_and_most_hits_failing_the_fdr():
    causes = _causes(id_diagnosis.diagnose(_counts(q_estimated=False, kept=0)))
    assert 'no q-values' in causes['PrSMs at q ≤ 0.01']
    causes = _causes(id_diagnosis.diagnose(_counts(kept=50)))
    assert 'decoys' in causes['PrSMs at q ≤ 0.01']


def test_calibration_offset():
    causes = _causes(id_diagnosis.diagnose(_counts(calibration_ppm=-7.5)))
    assert '-7.5 ppm' in causes['matched to a proteoform']


def test_at_most_three_causes_the_most_telling_first():
    counts = _counts(ms1=10, deconvolved=400, precursor=100, tagged=50, matched=1, kept=0)
    table = id_diagnosis.diagnose(counts)
    named = [s for s, c in zip(table['step'], table['cause']) if c]
    assert len(named) == id_diagnosis.MAX_CAUSES
    assert 'with deconvolved masses' in named and 'with a precursor mass' not in named
    assert '5 steps' in id_diagnosis.headline(counts)


def test_steps_of_a_run_deconvolved_only_leave_out_the_identification_steps():
    counts = _counts(tagged=None, matched=None, kept=None, q_estimated=None)
    assert id_diagnosis.diagnose(counts)['step'].to_list() == [
        'MS2 spectra', 'with deconvolved masses', 'with a precursor mass']


def test_lysate_example(tmp_path):
    counts = id_diagnosis.steps(FileManager(tmp_path, cache_path=WORKSPACE / 'cache'), LYSATE)
    assert counts['ms2'] > 0
    assert counts['ms2'] >= counts['deconvolved'] >= counts['tagged'] >= counts['matched'] >= counts['kept'] > 0
    assert counts['q_estimated'] is True
    assert counts['tolerance_ppm'] == 10.0
    table = id_diagnosis.diagnose(counts)
    assert table.height == len(id_diagnosis.STEP_LABELS)
