"""Fragmentation method of a run, for the Sequence View's ion series."""
from pathlib import Path

import pyopenms as oms
import pytest

from src.render.fragmentation import method_of, mzml_method, run_method

AM = oms.Precursor.ActivationMethod
DEMO_MZML = Path(
    'example-data/workspaces/default/flashdeconv/cache/files/example_fd/out_deconv.mzML'
)


@pytest.mark.parametrize('methods, expected', [
    ({AM.CID}, 'CID/HCD'),
    ({AM.HCD}, 'CID/HCD'),
    ({AM.ETD}, 'ETD/ECD'),
    ({AM.ECD}, 'ETD/ECD'),
    ({AM.EThcD}, 'EThcD'),
    ({AM.ETciD}, 'EThcD'),
    ({AM.ETD, AM.HCD}, 'EThcD'),
    ({AM.ETD, AM.CID}, 'EThcD'),
    ({AM.PD}, None),
    (set(), None),
])
def test_method_of(methods, expected):
    assert method_of(methods) == expected


def _experiment(per_spectrum):
    experiment = oms.MSExperiment()
    ms1 = oms.MSSpectrum()
    ms1.setMSLevel(1)
    experiment.addSpectrum(ms1)
    for methods in per_spectrum:
        spectrum = oms.MSSpectrum()
        spectrum.setMSLevel(2)
        precursor = oms.Precursor()
        precursor.setActivationMethods(set(methods))
        spectrum.setPrecursors([precursor])
        experiment.addSpectrum(spectrum)
    return experiment


def test_run_method_most_common():
    experiment = _experiment([{AM.ETD}, {AM.ETD}, {AM.HCD}, {AM.PD}])
    assert run_method(experiment) == 'ETD/ECD'


def test_run_method_tie_prefers_collisional():
    assert run_method(_experiment([{AM.ETD}, {AM.HCD}])) == 'CID/HCD'


def test_run_method_unknown():
    assert run_method(_experiment([{AM.PD}, set()])) is None
    assert run_method(_experiment([])) is None


def test_mzml_method(tmp_path):
    path = tmp_path / 'etd.mzML'
    oms.MzMLFile().store(str(path), _experiment([{AM.ETD}, {AM.ETD}]))
    assert mzml_method(path) == 'ETD/ECD'
    assert mzml_method(tmp_path / 'missing.mzML') is None


@pytest.mark.skipif(not DEMO_MZML.exists(), reason='demo workspace not present')
def test_demo_run_is_collisional():
    assert mzml_method(DEMO_MZML) == 'CID/HCD'
