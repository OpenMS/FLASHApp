"""Mass shifts of identified proteoforms (src/render/mass_shifts.py).

Hand-made PrSMs, plus a run on the identified lysate example of the demo
workspace shipped in example-data/.
"""
import os
import shutil
import sys
from pathlib import Path

import polars as pl
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.render import hierarchy_data, mass_shifts  # noqa: E402
from src.workflow.FileManager import FileManager  # noqa: E402

WORKSPACE = Path(__file__).resolve().parent.parent / 'example-data' / 'workspaces' / 'default'
LYSATE = 'example_lysate_identified'


def _hits(rows):
    defaults = {
        'index': 0, 'scan_idx': 0, 'Scan': 1, 'proteoform_id': 1, 'accession': 'P1',
        'RetentionTime': 100.0, 'Score': 10.0, 'Proforma': 'PEPTIDE',
        'ModMass': None, 'ModStart': None, 'ModEnd': None,
    }
    return pl.DataFrame(
        [{**defaults, 'index': i, 'proteoform_id': i, **row} for i, row in enumerate(rows)],
        schema_overrides={'ModMass': pl.String, 'ModStart': pl.String, 'ModEnd': pl.String},
    )


def test_one_row_per_shift_and_none_for_unmodified():
    hits = _hits([
        {'ModMass': '42.0106', 'ModStart': '2', 'ModEnd': '2'},
        {'ModMass': '15.9949;79.9663', 'ModStart': '5;9', 'ModEnd': '5;12'},
        {},
    ])
    shifts = mass_shifts.prsm_shifts(hits)
    assert shifts['shift'].to_list() == pytest.approx([42.0106, 15.9949, 79.9663])
    assert shifts['site'].to_list() == ['2', '5', '9–12']
    assert shifts['index'].to_list() == [0, 1, 1]


def test_numeric_mod_mass_column():
    hits = _hits([{}]).with_columns(pl.lit(15.9949).alias('ModMass'), pl.lit(3.0).alias('ModStart'),
                                     pl.lit(3.0).alias('ModEnd'))
    shifts = mass_shifts.prsm_shifts(hits)
    assert shifts['shift'].to_list() == pytest.approx([15.9949])
    assert shifts['site'].to_list() == ['3']


def test_groups_matched_to_unimod():
    # Intact masses: oxidations measured up to 30 mDa off still group and match
    hits = _hits([
        {'ModMass': '15.9949'}, {'ModMass': '16.0185'}, {'ModMass': '16.0230'},
        {'ModMass': '42.0106'}, {'ModMass': '41.9914'},
        {'ModMass': '83.826'},
    ])
    groups, histogram, prsms = mass_shifts.summary(mass_shifts.prsm_shifts(hits))
    assert groups['prsms'].to_list() == [3, 2, 1]
    assert groups['interpretation'].to_list()[:2] == ['Oxidation', 'Acetyl']
    # Nothing in Unimod at 83.83 Da, also after an isotope offset
    assert groups['interpretation'].to_list()[2] == ''
    assert prsms['shift_id'].to_list() == [0, 0, 0, 1, 1, 2]
    labelled = histogram.filter(pl.col('peak'))
    assert labelled.height == 3
    assert any(label.endswith('Oxidation') for label in labelled['label'])
    assert histogram['prsms'].sum() == 6


def test_isotope_offset_named():
    # Oxidation (15.9949) measured one 13C peak too heavy: preferred over the
    # amino acid substitution Unimod lists at that mass
    groups, _, _ = mass_shifts.summary(mass_shifts.prsm_shifts(_hits([{'ModMass': '16.9983'}])))
    assert groups['interpretation'][0] == 'Oxidation (+1 isotope peak)'
    assert groups['unimod'][0] == pytest.approx(16.9983, abs=1e-3)


def test_no_shifts():
    groups, histogram, prsms = mass_shifts.summary(mass_shifts.prsm_shifts(_hits([{}])))
    assert groups.height == histogram.height == prsms.height == 0


def test_lysate_demo(tmp_path):
    cache = tmp_path / 'flashtnt'
    cache.mkdir()
    shutil.copy(WORKSPACE / 'cache' / 'cache.db', cache / 'cache.db')
    (cache / 'files').symlink_to(WORKSPACE / 'cache' / 'files')
    hits = hierarchy_data.prsms(FileManager(cache, cache), LYSATE)
    groups, histogram, prsms = mass_shifts.summary(mass_shifts.prsm_shifts(hits))
    assert groups.height > 0
    assert groups['prsms'].sum() == prsms.height == histogram['prsms'].sum()
    top = groups.row(0, named=True)
    assert top['interpretation'] == 'Oxidation'
    assert set(prsms['index']) <= set(hits['index'])
