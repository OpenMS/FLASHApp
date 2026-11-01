"""Fragmentation method per scan, for the Sequence View's ion series."""
from pathlib import Path

import polars as pl
import pyopenms as oms
import pytest

from src.render.fragmentation import method_of, mzml_scan_methods, spectrum_method
from src.render.insight_data import with_fragmentation

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


def _spectrum(level, methods=()):
    spectrum = oms.MSSpectrum()
    spectrum.setMSLevel(level)
    if level > 1:
        precursor = oms.Precursor()
        precursor.setActivationMethods(set(methods))
        spectrum.setPrecursors([precursor])
    return spectrum


def test_spectrum_method():
    assert spectrum_method(_spectrum(1)) is None
    assert spectrum_method(_spectrum(2, {AM.ETD})) == 'ETD/ECD'
    assert spectrum_method(_spectrum(2, {AM.PD})) is None
    assert spectrum_method(_spectrum(2)) is None


def test_mzml_scan_methods_mixed_run(tmp_path):
    """Alternating HCD and ETD scans each get their own method."""
    experiment = oms.MSExperiment()
    for spectrum in [_spectrum(1), _spectrum(2, {AM.HCD}), _spectrum(2, {AM.ETD}),
                     _spectrum(1), _spectrum(2, {AM.ETD, AM.HCD})]:
        experiment.addSpectrum(spectrum)
    path = tmp_path / 'mixed.mzML'
    oms.MzMLFile().store(str(path), experiment)
    assert mzml_scan_methods(path) == (None, 'CID/HCD', 'ETD/ECD', None, 'EThcD')
    assert mzml_scan_methods(tmp_path / 'missing.mzML') == ()


def test_with_fragmentation_joins_by_scan_index():
    lf = pl.LazyFrame({'scan_idx': [2, 0, 7, None], 'x': [1, 2, 3, 4]},
                      schema={'scan_idx': pl.Int64, 'x': pl.Int64})
    out = with_fragmentation(lf, (None, 'CID/HCD', 'ETD/ECD'), 'scan_idx').collect()
    assert out['x'].to_list() == [1, 2, 3, 4]
    assert out['fragmentation'].to_list() == ['ETD/ECD', None, None, None]
    empty = with_fragmentation(lf, (), 'scan_idx').collect()
    assert empty['fragmentation'].null_count() == 4


@pytest.mark.skipif(not DEMO_MZML.exists(), reason='demo workspace not present')
def test_demo_scans_are_collisional():
    methods = mzml_scan_methods(DEMO_MZML)
    assert len(methods) == 317
    assert set(methods) == {None, 'CID/HCD'}
