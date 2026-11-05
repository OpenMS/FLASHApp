"""Theoretical isotope envelope of the selected mass (src/render/envelope_overlay.py)."""

import shutil
from pathlib import Path

import pytest

from src.render import envelope_overlay
from src.render.validation_data import ISOTOPE_SPACING, PROTON_MASS

MASS = 10000.0
ABUNDANCES = [0.5, 1.0, 0.8, 0.3, 0.005]


def _mz(isotope, charge):
    return (MASS + isotope * ISOTOPE_SPACING) / charge + PROTON_MASS


def _signal(charge, intensities):
    return [[0.0, _mz(i, charge), value, float(charge)] for i, value in intensities.items()]


def test_scaled_to_each_charge():
    peaks = _signal(10, {0: 50.0, 1: 100.0, 2: 80.0}) + _signal(11, {1: 20.0, 2: 16.0})
    result = envelope_overlay.theoretical_peaks(MASS, peaks, ABUNDANCES)

    by_charge = {}
    for peak in result:
        by_charge.setdefault(peak['charge'], []).append(peak)
    assert sorted(by_charge) == [10, 11]
    # Isotope 4 is below 1 % of the highest isotope: not drawn
    assert [p['isotope'] for p in by_charge[10]] == [0, 1, 2, 3]
    assert by_charge[10][1]['intensity'] == pytest.approx(100.0)
    assert by_charge[10][3]['intensity'] == pytest.approx(30.0)
    assert by_charge[11][1]['intensity'] == pytest.approx(20.0)
    assert by_charge[11][0]['mz'] == pytest.approx(_mz(0, 11))


def test_abundances_are_normalised():
    peaks = _signal(5, {1: 200.0})
    result = envelope_overlay.theoretical_peaks(MASS, peaks, [x * 7 for x in ABUNDANCES])
    assert result[1]['intensity'] == pytest.approx(200.0)


def test_nothing_to_draw():
    assert envelope_overlay.theoretical_peaks(MASS, [], ABUNDANCES) == []
    assert envelope_overlay.theoretical_peaks(MASS, None, ABUNDANCES) == []
    assert envelope_overlay.theoretical_peaks(None, _signal(5, {0: 1.0}), ABUNDANCES) == []
    assert envelope_overlay.theoretical_peaks(MASS, _signal(5, {0: 1.0}), []) == []
    # Charge 0 and peaks beyond the envelope are ignored
    assert envelope_overlay.theoretical_peaks(MASS, [[0, 1000.0, 5.0, 0.0]], ABUNDANCES) == []
    assert envelope_overlay.theoretical_peaks(MASS, _signal(5, {9: 1.0}), ABUNDANCES) == []


DEMO = Path('example-data/workspaces/default/cache')


@pytest.mark.skipif(not DEMO.exists(), reason='demo workspace not present')
def test_demo_workspace(tmp_path):
    from src.workflow.FileManager import FileManager

    # A copy of the cache index, so the demo workspace is not touched
    cache = tmp_path / 'cache'
    cache.mkdir()
    shutil.copy(DEMO / 'cache.db', cache / 'cache.db')
    (cache / 'files').symlink_to(DEMO.resolve() / 'files')
    file_manager = FileManager(cache, cache)
    datasets = file_manager.get_results_list(['threedim_SN_plot', 'mass_table'])
    assert datasets
    peaks = envelope_overlay.selected_mass_envelope(file_manager, datasets[0], 0, 0)
    assert peaks
    assert all(p['intensity'] > 0 and p['charge'] > 0 for p in peaks)
    assert envelope_overlay.selected_mass_envelope(file_manager, datasets[0], 0, 10 ** 6) == []
    assert envelope_overlay.selected_mass_envelope(file_manager, datasets[0], None, 0) == []
