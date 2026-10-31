"""Tests for the frames of the PrSM validation views (src/render/validation_data.py).

Unit tests on hand-made entries, plus runs on the demo workspace shipped in
example-data/.
"""

import os
import shutil
import sys
from pathlib import Path

import numpy as np
import polars as pl
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.render import validation_data  # noqa: E402
from src.workflow.FileManager import FileManager  # noqa: E402

WORKSPACE = Path(__file__).resolve().parent.parent / 'example-data' / 'workspaces' / 'default'
FD_DATASET = 'example_fd'
TNT_DATASET = 'example_spectrum_aqpz'


def _file_manager(tmp_path, tool):
    cache = tmp_path / tool
    cache.mkdir()
    shutil.copy(WORKSPACE / tool / 'cache' / 'cache.db', cache / 'cache.db')
    (cache / 'files').symlink_to(WORKSPACE / tool / 'cache' / 'files')
    return FileManager(cache, cache)


@pytest.fixture
def fd(tmp_path):
    return _file_manager(tmp_path, 'flashdeconv')


@pytest.fixture
def tnt(tmp_path):
    return _file_manager(tmp_path, 'flashtnt')


# --- match_fragments ---------------------------------------------------------

def _entry(start=10, end=19):
    return {
        'sequence': list('A' * 30),
        'proteoform_start': start,
        'proteoform_end': end,
        # Fragment lengths 1, 2, 3
        'fragment_masses_b': [[100.0], [200.0], [300.0]],
        'fragment_masses_y': [[150.0], [250.0, 260.0], [350.0]],
    }


def test_match_fragments_within_tolerance():
    observed = [200.0 * (1 + 4e-6), 250.0 * (1 - 2e-6), 999.0]
    rows = validation_data.match_fragments(_entry(), observed, 5.0, ['b', 'y'])
    by_name = {r['fragment']: r for r in rows}
    assert set(by_name) == {'b2', 'y2'}
    assert by_name['b2']['ppm_error'] == pytest.approx(4.0, abs=1e-6)
    assert by_name['y2']['ppm_error'] == pytest.approx(-2.0, abs=1e-6)
    # Prefix ions: residue after the bond counted from the start; suffix ions
    # from the end of the proteoform
    assert by_name['b2']['residue'] == 11
    assert by_name['y2']['residue'] == 18


def test_match_fragments_skips_undetermined_terminus():
    observed = [200.0, 250.0]
    rows = validation_data.match_fragments(_entry(end=-1), observed, 5.0, ['b', 'y'])
    assert {r['ion_type'] for r in rows} == {'b'}
    rows = validation_data.match_fragments(_entry(start=-1), observed, 5.0, ['b', 'y'])
    assert {r['ion_type'] for r in rows} == {'y'}


def test_match_fragments_without_masses():
    assert validation_data.match_fragments(_entry(), [], 5.0, ['b']) == []
    assert validation_data.match_fragments(_entry(), None, 5.0, ['b']) == []


def test_annotated_fragments():
    fragments = [
        {'name': 'b5', 'ion_type': 'b', 'kind': 'N', 'start': 10, 'end': 14,
         'charge': 2, 'mono_mass': 500.0, 'observed_mass': 500.001, 'ppm_error': 2.0},
        {'name': 'y[3-5]', 'ion_type': 'b/y', 'kind': 'I', 'start': 13, 'end': 15,
         'charge': 1, 'mono_mass': 300.0, 'observed_mass': 300.0, 'ppm_error': -1.0},
    ]
    rows = validation_data.annotated_fragment_errors(fragments)
    assert [r['ion_type'] for r in rows] == ['b', 'internal']
    assert rows[0]['fragment'] == 'b5 (2+)'
    assert [r['residue'] for r in rows] == [14, 13]


# --- Averagine ---------------------------------------------------------------

def test_averagine_envelope_shape():
    averagine = validation_data.Averagine()
    small = averagine.envelope(1000.0)
    large = averagine.envelope(20000.0)
    assert small.max() == pytest.approx(1.0)
    # The monoisotopic peak is the highest for small, not for large proteins
    assert int(np.argmax(small)) == 0
    assert int(np.argmax(large)) > 5
    # The tail below 0.1 % is cut
    assert large[-1] >= validation_data.MIN_RELATIVE_ABUNDANCE


# --- On the demo results ------------------------------------------------------

def test_fragment_errors_on_demo(tnt):
    errors = validation_data.fragment_errors(tnt, TNT_DATASET).collect()
    assert errors.height > 0
    assert set(errors.columns) == set(validation_data.FRAGMENT_ERROR_SCHEMA)
    # Matched within the FLASHTnT tolerance of the demo (5 ppm)
    assert errors['ppm_error'].abs().max() <= 5.0 + 1e-6
    assert set(errors['ion_type'].unique()) <= {'b', 'c', 'y', 'z'}


def test_isotope_envelopes_on_demo(fd):
    envelopes = validation_data.isotope_envelopes(fd, FD_DATASET).collect()
    assert set(envelopes['side'].unique()) == {'observed', 'theoretical'}
    # Every side of every mass is scaled to a highest isotope of 1
    tops = envelopes.group_by(['scan_idx', 'mass_idx', 'side']).agg(pl.col('abundance').max())
    assert tops['abundance'].min() == pytest.approx(1.0)
    assert tops['abundance'].max() == pytest.approx(1.0)
    # Fragmentation spectra get their precursor envelope as mass -1
    precursors = envelopes.filter(pl.col('mass_idx') == -1)
    assert precursors.height > 0
    scans = pl.read_parquet(WORKSPACE / 'flashdeconv' / 'cache' / 'files' / 'example_fd' / 'scan_table.pq')
    ms2 = set(scans.filter(pl.col('MSLevel') == 2)['index'].cast(pl.Int64).to_list())
    assert set(precursors['scan_idx'].unique().to_list()) <= ms2


def test_precursor_xics_on_demo(fd):
    xics = validation_data.precursor_xics(fd, FD_DATASET).collect()
    assert xics.height > 0
    one = xics.filter(pl.col('scan_idx') == xics['scan_idx'][0])
    # The sum over the charges is the 'All charges' trace
    total = one.filter(pl.col('charge') == 0).sort('rt')
    per_charge = (
        one.filter(pl.col('charge') > 0).group_by('rt').agg(pl.col('intensity').sum()).sort('rt')
    )
    assert total['rt'].to_list() == per_charge['rt'].to_list()
    assert np.allclose(total['intensity'].to_numpy(), per_charge['intensity'].to_numpy())
    assert set(one.filter(pl.col('charge') == 0)['series']) == {'All charges'}


def test_precursor_xics_without_ms1(tnt):
    # The FLASHTnT demo is a single fragmentation spectrum
    assert validation_data.precursor_xics(tnt, TNT_DATASET).collect().is_empty()


# --- Characterization ----------------------------------------------------------

def _fragments(*pairs):
    return [{'ion_type': ion, 'residue': residue} for ion, residue in pairs]


def test_explained_bonds():
    # b ends at its residue, y starts at it: both explain bond 4 (residues 4|5)
    bonds = validation_data.explained_bonds(_fragments(('b', 4), ('y', 5), ('c', 7), ('internal', 3)))
    assert bonds == {4, 7}


def test_characterize_levels():
    entry = {'sequence': list('A' * 11), 'proteoform_start': 0, 'proteoform_end': 10,
             'modifications': [{'start': 3, 'end': 3}]}
    # 3 of 10 bonds, termini known, the one modification localized
    full = validation_data.characterize(entry, _fragments(('b', 0), ('b', 1), ('y', 9)))
    assert full == {'Cleavages': 30.0, 'LocalizedMods': '1/1', 'Characterization': 'Characterized'}
    # An unlocalized mass shift
    entry_shift = dict(entry, modifications=[{'start': 3, 'end': 6}])
    partial = validation_data.characterize(entry_shift, _fragments(('b', 0), ('b', 1), ('y', 9)))
    assert partial['Characterization'] == 'Partial'
    assert partial['LocalizedMods'] == '0/1'
    # Undetermined C-terminus
    entry_open = dict(entry, proteoform_end=-1)
    assert validation_data.characterize(entry_open, _fragments(('b', 0), ('b', 1), ('b', 2)))['Characterization'] == 'Partial'
    # No fragments
    assert validation_data.characterize(entry, [])['Characterization'] == 'Ambiguous'
    assert validation_data.characterize(dict(entry, modifications=[]), [])['LocalizedMods'] == '-'


def test_characterization_on_demo(tnt):
    table = validation_data.characterization(tnt, TNT_DATASET).collect()
    assert table.height == 17
    assert set(table['Characterization']) <= {'Characterized', 'Partial', 'Ambiguous'}
    # The best hit of the demo (AQPZ, score 505) is identified but carries an
    # unlocalized mass shift
    best = table.filter(pl.col('proteoform_index') == 9).row(0, named=True)
    assert best['Characterization'] == 'Partial'
    assert best['LocalizedMods'] == '0/1'
    assert best['Cleavages'] > 25
