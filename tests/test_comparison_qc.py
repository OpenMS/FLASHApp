"""What the identification and quality-control views of a comparison count."""

import polars as pl
import pytest

from src.quant import qc

RUNS = ['A1', 'A2', 'B1']


def _long():
    rows = [
        # id, sample, intensity, status, score, qvalue
        ('p1', 'A1', 100.0, 'detected', None, None),
        ('p1', 'A2', 120.0, 'detected', None, None),
        ('p1', 'B1', 300.0, 'detected', None, None),
        ('p2', 'A1', 50.0, 'detected', None, None),
        ('p2', 'A2', 40.0, 'transferred', 0.9, 0.001),
        ('p2', 'B1', None, 'missing', 0.1, 0.8),
        ('p3', 'A1', None, 'missing', None, None),
        ('p3', 'A2', None, 'missing', None, None),
        ('p3', 'B1', 80.0, 'detected', None, None),
        ('p4', 'A1', None, 'missing', None, None),
        ('p4', 'A2', None, 'missing', None, None),
        ('p4', 'B1', None, 'missing', None, None),
    ]
    return pl.DataFrame(rows, schema={'id': pl.Utf8, 'sample': pl.Utf8, 'intensity': pl.Float64,
                                      'status': pl.Utf8, 'score': pl.Float64, 'qvalue': pl.Float64}, orient='row')


def test_identification_level_prefers_proteoforms():
    assert qc.identification_level(['feature', 'proteoform', 'protein'], 'feature') == 'proteoform'
    assert qc.identification_level(['feature', 'proteoform', 'protein'], 'protein') == 'protein'
    assert qc.identification_level(['feature'], 'feature') == 'feature'


def test_counts_per_run_by_how_they_were_found():
    counts = qc.counts_per_run(_long(), RUNS)
    got = {(r['sample'], r['status']): r['count'] for r in counts.iter_rows(named=True)}
    assert got[('A1', 'detected')] == 2
    assert got[('A2', 'detected')] == 1 and got[('A2', 'matched')] == 1
    assert got[('B1', 'detected')] == 2 and got[('B1', 'matched, identified')] == 0
    assert counts['sample'].to_list()[:3] == ['A1', 'A1', 'A1']


def test_frequency_counts_entities_by_runs():
    table = qc.frequency(_long(), RUNS)
    assert table['detected'].to_list() == [2, 0, 1]   # p2, p3 in one run; p1 in three
    assert table['found'].to_list() == [1, 1, 1]      # p3 in one, p2 in two, p1 in three


def test_overlap_is_a_jaccard_index():
    frame = qc.overlap(_long(), RUNS)
    by = {r['Sample']: r for r in frame.iter_rows(named=True)}
    assert by['A1']['A1'] == 100.0
    assert by['A1']['A2'] == pytest.approx(50.0)       # {p1, p2} vs {p1}
    assert by['A2']['B1'] == pytest.approx(50.0)       # {p1} vs {p1, p3}
    assert by['A1']['B1'] == by['B1']['A1']


def test_identification_table_lists_found_entities_most_often_found_first():
    annotations = pl.DataFrame({'id': ['p1', 'p2', 'p3', 'p4'], 'label': ['one', 'two', 'three', 'four'],
                                'accession': ['sp|X|ONE', 'sp|Y|TWO', None, None]})
    table = qc.identification_table(_long(), annotations, RUNS)
    assert table['label'].to_list() == ['one', 'two', 'three']
    assert table.row(1, named=True)['A2'] == 'matched'
    assert table.row(1, named=True)['B1'] == 'missing'
    assert table['protein_label'].to_list()[0] == 'ONE'


def test_totals_and_runs_loaded_unlike_their_replicates():
    frame = qc.totals(_long(), RUNS)
    assert frame['total'].to_list() == [150.0, 160.0, 380.0]
    assert frame['quantified'].to_list() == [2, 2, 2]
    assert qc.unlike_replicates(frame, {'A1': 'A', 'A2': 'A', 'B1': 'B'}) == []
    low = pl.DataFrame({'sample': ['A1', 'A2', 'A3'], 'total': [10.0, 100.0, 110.0]})
    assert qc.unlike_replicates(low, {'A1': 'A', 'A2': 'A', 'A3': 'A'}) == ['A1']


def test_correlation_needs_shared_values():
    frame = qc.correlation(_long(), RUNS, min_shared=2)
    by = {r['Sample']: r for r in frame.iter_rows(named=True)}
    assert by['A1']['A1'] == 1.0
    assert by['A1']['A2'] == pytest.approx(1.0)          # two shared values always correlate
    assert by['A1']['B1'] is None                         # one shared value
    assert by['A2']['B1'] is None


def test_coefficients_of_variation_within_conditions():
    wide = pl.DataFrame({'A1': [100.0, 50.0, None], 'A2': [100.0, 150.0, 10.0], 'B1': [1.0, 2.0, 3.0]})
    frame = qc.coefficients_of_variation(wide, {'A1': 'A', 'A2': 'A', 'B1': 'B'})
    assert frame['condition'].to_list() == ['A', 'A']    # B has one replicate, the third row one value
    assert frame['cv'].to_list()[0] == pytest.approx(0.0)
    assert frame['cv'].to_list()[1] == pytest.approx(70.71, abs=0.01)


def test_match_quality_keeps_scored_matches_only():
    frame = qc.match_quality(_long(), RUNS)
    assert frame.select('sample', 'outcome').rows() == [('A2', 'matched'), ('B1', 'not matched')]
