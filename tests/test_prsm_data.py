"""Tests for the PrSM mass table frame (src/render/prsm_data.py)."""

import os
import shutil
import sys
from pathlib import Path

import polars as pl
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.render import prsm_data, validation_data  # noqa: E402
from src.workflow.FileManager import FileManager  # noqa: E402

WORKSPACE = Path(__file__).resolve().parent.parent / 'example-data' / 'workspaces' / 'default'
TNT_DATASET = 'example_native_identified'


@pytest.fixture
def tnt(tmp_path):
    cache = tmp_path / 'flashtnt'
    cache.mkdir()
    shutil.copy(WORKSPACE / 'cache' / 'cache.db', cache / 'cache.db')
    (cache / 'files').symlink_to(WORKSPACE / 'cache' / 'files')
    return FileManager(cache, cache)


def _entry():
    return {
        'sequence': list('A' * 30),
        'proteoform_start': 10,
        'proteoform_end': 19,
        'fragment_masses_b': [[100.0], [200.0], [300.0]],
        'fragment_masses_y': [[150.0], [200.0 * (1 + 3e-6)], [350.0]],
    }


def test_mass_matches_keyed_by_mass_index():
    masses = [999.0, 200.0 * (1 + 1e-6), 350.0, None]
    matches = prsm_data.mass_matches(_entry(), masses, 5.0, ['b', 'y'])
    assert set(matches) == {1, 2}
    # One mass explaining two ions keeps both
    assert {m['fragment'] for m in matches[1]} == {'b2', 'y2'}
    assert [m['fragment'] for m in matches[2]] == ['y3']


def test_mass_matches_annotated_fragments_by_nearest_mass():
    entry = _entry()
    entry['annotated_fragments'] = [
        {'kind': 'N', 'ion_type': 'b', 'name': 'b2', 'charge': 3, 'end': 11,
         'mono_mass': 200.0, 'observed_mass': 200.0 * (1 + 2e-6), 'ppm_error': 2.0},
        # Too far from every deconvolved mass
        {'kind': 'C', 'ion_type': 'y', 'name': 'y9', 'charge': 2, 'start': 12,
         'mono_mass': 900.0, 'observed_mass': 900.0, 'ppm_error': 0.0},
    ]
    matches = prsm_data.mass_matches(entry, [150.0, 200.0, 800.0], 5.0, ['b', 'y'])
    assert set(matches) == {1}
    assert matches[1][0]['fragment'] == 'b2 (3+)'


def test_row_takes_numbers_of_closest_match():
    matches = [
        {'fragment': 'y2', 'residue': 18, 'theoretical_mass': 200.0006, 'ppm_error': 3.0},
        {'fragment': 'b2', 'residue': 11, 'theoretical_mass': 200.0, 'ppm_error': -1.0},
    ]
    row = prsm_data._row(4, 7, 199.9998, {'charges': '2–3'}, matches)
    assert row['mass_id'] == 8
    assert row['status'] == prsm_data.MATCHED
    assert row['ions'] == 'b2, y2'
    assert row['position'] == 12
    assert row['ppm_error'] == -1.0
    assert row['mass_error'] == pytest.approx(-0.0002)
    unmatched = prsm_data._row(4, 0, 500.0, {}, None)
    assert unmatched['status'] == prsm_data.NOT_MATCHED
    assert unmatched['ions'] is None and unmatched['ppm_error'] is None


def test_prsm_masses_on_demo(tnt):
    masses = prsm_data.prsm_masses(tnt, TNT_DATASET).collect()
    assert list(masses.columns) == list(prsm_data.PRSM_MASS_SCHEMA)
    # Every mass of the spectrum, for every proteoform
    counts = masses.group_by('proteoform_index').len()
    assert counts['len'].min() == counts['len'].max() > 0
    assert masses['mass_id'].equals(masses['mass_idx'] + 1)
    matched = masses.filter(pl.col('matched'))
    assert matched.height > 0
    assert (masses['status'] == prsm_data.MATCHED).sum() == matched.height
    assert matched['ppm_error'].abs().max() <= 5.0 + 1e-6
    assert masses.filter(~pl.col('matched'))['ions'].null_count() == masses.height - matched.height
    # The mass table supplies intensity and charges
    assert masses['SumIntensity'].null_count() == 0
    assert masses['charges'].null_count() == 0

    # The fragment error plot shows the same proteoforms. Its fragments come
    # from FLASHAnnotator's envelopes in the demo, so they are not the
    # deconvolved masses of this table.
    errors = validation_data.fragment_errors(tnt, TNT_DATASET).collect()
    assert set(errors['proteoform_index']) == set(matched['proteoform_index'])
