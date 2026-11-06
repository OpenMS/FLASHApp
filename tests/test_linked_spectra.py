"""Precursor spectrum of the selected spectrum (src/render/linked_spectra.py)."""

import shutil
from pathlib import Path

import polars as pl
import pytest

from src.render import linked_spectra


def _scans():
    return pl.DataFrame({
        'index': [0, 1, 2, 3, 4, 5],
        'Scan': [100, 101, 102, 103, 104, 105],
        'MSLevel': [1, 2, 2, 1, 2, 2],
        'PrecursorMass': [0.0, 1000.0, 2000.0, 0.0, 3000.0, 4000.0],
    })


def _links(**kwargs):
    links = linked_spectra.link_surveys(_scans(), **kwargs)
    return {row['scan_idx']: (row['survey_idx'], row['survey_mass_idx']) for row in links.to_dicts()}


def test_previous_ms1_without_precursor_scans():
    links = _links()

    assert links[0] == (0, None)
    assert links[1] == (0, None)
    assert links[2] == (0, None)
    assert links[3] == (3, None)
    assert links[5] == (3, None)


def test_named_precursor_scan_wins():
    precursor_scans = pl.DataFrame({'index': [1, 2, 4, 5], 'PrecursorScan': [103.0, 0.0, 100.0, 999.0]})
    links = _links(precursor_scans=precursor_scans)

    assert links[1] == (3, None)
    # 0: unknown, 999: no such spectrum; both fall back to the previous MS1
    assert links[2] == (0, None)
    assert links[4] == (0, None)
    assert links[5] == (3, None)


def test_precursor_mass_is_the_closest_within_the_tolerance():
    survey_masses = pl.DataFrame({
        'scan_idx': [0, 0, 0, 3],
        'mass_idx': [0, 1, 2, 0],
        'MonoMass': [1000.008, 1000.002, 2000.5, 3000.001],
    })
    links = _links(survey_masses=survey_masses)

    assert links[1] == (0, 1)
    # 0.5 Da off: not the precursor
    assert links[2] == (0, None)
    # The precursor spectrum of scan 4 is scan 3
    assert links[4] == (3, 0)
    # An MS1 spectrum has no precursor mass
    assert links[0] == (0, None)


def test_spectra_before_the_first_ms1_have_none():
    scans = _scans().with_columns(pl.Series('MSLevel', [2, 2, 2, 1, 2, 2]))
    links = linked_spectra.link_surveys(scans)

    assert links.filter(pl.col('scan_idx') < 3)['survey_idx'].is_null().all()
    assert linked_spectra.survey_of(links, 4) == (3, None)
    assert linked_spectra.survey_of(links, 0) == (None, None)
    assert linked_spectra.survey_of(links, None) == (None, None)
    assert linked_spectra.survey_of(links, 42) == (None, None)


def test_no_spectra():
    links = linked_spectra.link_surveys(pl.DataFrame())
    assert links.is_empty()
    assert list(links.schema) == list(linked_spectra.LINK_SCHEMA)


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
    datasets = file_manager.get_results_list(['scan_table'])
    assert datasets
    for dataset_id in datasets:
        links = linked_spectra.survey_links(file_manager, dataset_id)
        scans = pl.read_parquet(file_manager.get_result_path(dataset_id, 'scan_table'))
        assert links.height == scans.height
        levels = dict(zip(scans['index'].cast(pl.Int64), scans['MSLevel']))
        # Every precursor spectrum is an MS1 spectrum
        assert all(levels[s] == 1 for s in links['survey_idx'].drop_nulls())
        # FLASHDeconv found most precursor masses in their precursor spectrum
        # (a run without MS1 spectra has none)
        if 1 in levels.values():
            ms2 = links.filter(pl.col('scan_idx').is_in([i for i, lvl in levels.items() if lvl > 1]))
            assert ms2['survey_mass_idx'].is_not_null().sum() > ms2.height / 2
        else:
            assert links['survey_idx'].is_null().all()


def _precursor(mz, lower, upper, target=None):
    oms = pytest.importorskip('pyopenms')
    precursor = oms.Precursor()
    precursor.setMZ(mz)
    precursor.setIsolationWindowLowerOffset(lower)
    precursor.setIsolationWindowUpperOffset(upper)
    if target is not None:
        precursor.setMetaValue(linked_spectra.TARGET_MZ_KEY, target)
    return precursor


def test_precursor_window_around_the_target():
    assert linked_spectra.precursor_window(_precursor(1301.77, 2.5, 2.5, 1302.49)) == pytest.approx(
        (1299.99, 1304.99, 1302.49)
    )
    # Without a target, around the precursor m/z
    assert linked_spectra.precursor_window(_precursor(800.0, 1.0, 2.0)) == pytest.approx((799.0, 802.0, 800.0))
    # No offsets: no window
    assert linked_spectra.precursor_window(_precursor(800.0, 0.0, 0.0)) is None


def test_isolation_windows_of_a_missing_file(tmp_path):
    assert linked_spectra.isolation_windows(tmp_path / 'missing.mzML') == ()


@pytest.mark.skipif(not DEMO.exists(), reason='demo workspace not present')
def test_demo_isolation_windows(tmp_path):
    pytest.importorskip('pyopenms')
    from src.workflow.FileManager import FileManager

    cache = tmp_path / 'cache'
    cache.mkdir()
    shutil.copy(DEMO / 'cache.db', cache / 'cache.db')
    (cache / 'files').symlink_to(DEMO.resolve() / 'files')
    file_manager = FileManager(cache, cache)
    dataset_id = 'example_deconvolution'
    scans = pl.read_parquet(file_manager.get_result_path(dataset_id, 'scan_table'))
    levels = dict(zip(scans['index'].cast(pl.Int64), scans['MSLevel']))
    for scan_idx, level in levels.items():
        window = linked_spectra.isolation_window_of(file_manager, dataset_id, scan_idx)
        if level == 1:
            assert window is None
        else:
            lower, upper, target = window
            assert lower < target < upper
            assert upper - lower == pytest.approx(5.0)
    assert linked_spectra.isolation_window_of(file_manager, dataset_id, None) is None
