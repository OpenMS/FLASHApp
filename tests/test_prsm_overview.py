"""Overview columns of the protein table (src/render/prsm_overview.py)."""

import shutil
from pathlib import Path

import polars as pl
import pytest

from src.render import prsm_overview


def _hit(index, **values):
    row = {
        'index': index, 'Scan': 7, 'PrecursorMass': 10000.5, 'ProteoformMass': 10000.0,
        'StartPosition': 2, 'EndPosition': 120, 'length': 130,
        'ModMass': '15.994915', 'ModStart': '4', 'ModEnd': '4',
    }
    row.update(values)
    return row


def _row(frame, index):
    return frame.filter(pl.col('index') == index).to_dicts()[0]


def test_columns_of_a_hit():
    precursors = prsm_overview.precursors(pl.DataFrame({
        'ScanNum': ['7', '7'], 'PrecursorMz': ['1001.05', '1001.05'],
        'PrecursorCharge': ['10', '10'], 'PrecursorMonoisotopicMass': ['10000.4', '10000.4'],
    }))
    row = _row(prsm_overview.overview_columns(pl.DataFrame([_hit(3)]), precursors), 3)

    assert row['residues'] == '2–120 of 130'
    assert row['precursor_mz'] == pytest.approx(1001.05)
    assert row['precursor_charge'] == 10
    # FLASHTnT's precursor mass wins over the spectrum table's
    assert row['precursor_mass'] == pytest.approx(10000.5)
    assert row['delta_mass'] == pytest.approx(0.5)
    assert row['modification_sites'] == '+15.9949 Da at 4'


def test_undetermined_values_are_null():
    hit = _hit(0, PrecursorMass=-1.0, ProteoformMass=-1.0, EndPosition=-1, ModMass='', ModStart='', ModEnd='')
    row = _row(prsm_overview.overview_columns(pl.DataFrame([hit])), 0)

    assert row['residues'] == '2–? of 130'
    for column in ('precursor_mz', 'precursor_charge', 'precursor_mass', 'delta_mass', 'modification_sites'):
        assert row[column] is None, column


def test_precursor_mass_from_the_spectrum_table():
    precursors = prsm_overview.precursors(pl.DataFrame({
        'ScanNum': ['7'], 'PrecursorMonoisotopicMass': ['10000.4'],
    }))
    row = _row(prsm_overview.overview_columns(pl.DataFrame([_hit(1, PrecursorMass=-1.0)]), precursors), 1)

    assert row['precursor_mass'] == pytest.approx(10000.4)
    assert row['precursor_mz'] is None
    assert row['precursor_charge'] is None


def test_unlocalized_modifications_name_their_range():
    hit = _hit(1, ModMass='-2.006326;18.0', ModStart='16;30', ModEnd='97;30')
    row = _row(prsm_overview.overview_columns(pl.DataFrame([hit])), 1)

    assert row['modification_sites'] == '−2.0063 Da in 16–97 (not localized); +18.0 Da at 30'


def test_no_hits():
    frame = prsm_overview.overview_columns(pl.DataFrame())
    assert frame.is_empty()
    assert list(frame.schema) == list(prsm_overview.OVERVIEW_SCHEMA)


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
    datasets = file_manager.get_results_list(['protein_dfs'])
    assert datasets
    for dataset_id in datasets:
        frame = prsm_overview.protein_overview(file_manager, dataset_id).collect()
        assert frame['index'].n_unique() == frame.height > 0
        # The demo spectra carry an isolation m/z but no precursor charge
        assert frame['precursor_mz'].is_not_null().all()
        assert frame['precursor_charge'].is_null().all()
