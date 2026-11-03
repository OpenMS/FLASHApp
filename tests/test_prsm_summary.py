"""The PrSM summary rows (src/render/prsm_summary.py)."""

from pathlib import Path

import polars as pl
import pytest

from src.render import prsm_summary


def _hit(index, **values):
    row = {
        'index': index, 'Scan': 7, 'RetentionTime': 1234.567, 'accession': 'P1',
        'description': 'Protein one', 'PrecursorMass': 10000.5, 'ProteoformMass': 10000.0,
        'Proforma': 'PEP[+15.9949]TIDE', 'MatchingFragments': 12, 'Coverage(%)': 25.04,
        'StartPosition': 2, 'EndPosition': 120, 'length': 130, 'TagCount': 3,
        'ModCount': 1, 'ModMass': '15.994915', 'ModStart': '4', 'ModEnd': '4',
        'Score': 42, 'ProteoformLevelQvalue': 0.00123,
    }
    row.update(values)
    return row


def _values(frame, index):
    return dict(frame.filter(pl.col('index') == index).select('field', 'value').iter_rows())


def test_fields_of_a_hit_in_order():
    precursors = prsm_summary.precursors(pl.DataFrame({
        'ScanNum': ['7', '7'], 'PrecursorMz': ['1001.05', '1001.05'],
        'PrecursorCharge': ['10', '10'], 'PrecursorMonoisotopicMass': ['10000.4', '10000.4'],
    }))
    rows = prsm_summary.summary_rows(pl.DataFrame([_hit(3)]), precursors)

    assert rows['field'].to_list()[:3] == ['Protein', 'Proteoform', 'Residues']
    assert rows.height == prsm_summary.FIELD_COUNT
    assert rows['row'].n_unique() == rows.height
    values = _values(rows, 3)
    assert values['Protein'] == 'P1 (Protein one)'
    assert values['Proteoform'] == 'PEP[+15.9949]TIDE'
    assert values['Residues'] == '2–120 of 130'
    assert values['Precursor m/z'] == '1001.05'
    assert values['Precursor charge'] == '10'
    # FLASHTnT's precursor mass wins over the spectrum table's
    assert values['Precursor mass (Da)'] == '10000.5'
    assert values['Δ mass (Da)'] == '0.5'
    assert values['Sequence coverage (%)'] == '25.0'
    assert values['Modifications'] == '1 (1 localized)'
    assert values['Modification sites'] == '+15.9949 Da at 4'
    assert values['q-value (proteoform level)'] == '0.0012'


def test_undetermined_values_show_a_dash():
    hit = _hit(
        0, PrecursorMass=-1.0, ProteoformMass=-1.0, EndPosition=-1, ProteoformLevelQvalue=-1.0,
        ModCount=0, ModMass='', ModStart='', ModEnd='',
    )
    values = _values(prsm_summary.summary_rows(pl.DataFrame([hit])), 0)

    assert values['Residues'] == '2–? of 130'
    for field in ('Precursor m/z', 'Precursor charge', 'Precursor mass (Da)',
                  'Proteoform mass (Da)', 'Δ mass (Da)', 'q-value (proteoform level)',
                  'Modification sites'):
        assert values[field] == '-', field
    assert values['Modifications'] == 'none'


def test_unlocalized_modifications_name_their_range():
    hit = _hit(1, ModCount=2, ModMass='-2.006326;18.0', ModStart='16;30', ModEnd='97;30')
    values = _values(prsm_summary.summary_rows(pl.DataFrame([hit])), 1)

    assert values['Modifications'] == '2 (1 localized)'
    assert values['Modification sites'] == '−2.0063 Da in 16–97 (not localized); +18.0 Da at 30'


def test_spectrum_table_without_precursor_columns():
    precursors = prsm_summary.precursors(pl.DataFrame({'ScanNum': ['7'], 'PrecursorMz': ['900.5']}))
    values = _values(prsm_summary.summary_rows(pl.DataFrame([_hit(0)]), precursors), 0)

    assert values['Precursor m/z'] == '900.5'
    assert values['Precursor charge'] == '-'


def test_no_hits():
    assert prsm_summary.summary_rows(pl.DataFrame()).is_empty()


DEMO = Path('example-data/workspaces/default/flashtnt/cache')


@pytest.mark.skipif(not DEMO.exists(), reason='demo workspace not present')
def test_demo_workspace():
    from src.workflow.FileManager import FileManager

    file_manager = FileManager(DEMO.parent, DEMO)
    datasets = file_manager.get_results_list(['protein_dfs'])
    assert datasets
    for dataset_id in datasets:
        rows = prsm_summary.prsm_summary(file_manager, dataset_id).collect()
        hits = rows['index'].n_unique()
        assert rows.height == hits * prsm_summary.FIELD_COUNT
        # The demo spectra carry an isolation m/z but no precursor charge
        assert set(rows.filter(pl.col('field') == 'Precursor charge')['value']) == {'-'}
