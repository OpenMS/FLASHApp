"""Protein → proteoform → PrSM levels (src/render/hierarchy_data.py)."""

import polars as pl

from src.render import hierarchy_data


def _hits(qvalues=None):
    rows = [
        # index, Scan, accession, Proforma, Start, End, Score
        (0, 10, 'P1', 'ACDEK', 1, 5, 30.0),
        (1, 11, 'P1', 'ACDEK', 1, 5, 40.0),
        (2, 12, 'P1', 'CDEK', 2, 5, 20.0),
        (3, 13, 'P2', 'MKLV', 0, -1, 50.0),
        (4, 14, 'P2', '', 3, 9, 10.0),
    ]
    frame = pl.DataFrame(
        rows, orient='row',
        schema=['index', 'Scan', 'accession', 'Proforma', 'StartPosition', 'EndPosition', 'Score'],
    ).with_columns(
        pl.col('Scan').alias('scan_idx'),
        pl.lit('desc').alias('description'),
        pl.lit(1000.0).alias('ProteoformMass'),
        pl.lit(1001.0).alias('PrecursorMass'),
        pl.lit(1.5).alias('RetentionTime'),
        pl.lit(12).alias('MatchingFragments'),
        pl.lit(0).alias('ModCount'),
        pl.lit(0.0).alias('ModMass'),
        pl.Series('ProteoformLevelQvalue', qvalues if qvalues is not None else [-1.0] * len(rows)),
    )
    return frame


def test_groups_are_numbered_by_best_score():
    hits = hierarchy_data.assign_groups(_hits())
    groups = hierarchy_data.protein_groups_from(hits)
    # P2 has the best hit (50), so it is protein 1
    assert groups['accession'].to_list() == ['P2', 'P1']
    assert groups['protein_group'].to_list() == [1, 2]
    p1 = groups.filter(pl.col('accession') == 'P1').row(0, named=True)
    assert p1['proteoforms'] == 2
    assert p1['prsms'] == 3
    assert p1['best_score'] == 40.0
    assert p1['best_index'] == 1
    assert p1['best_qvalue'] is None


def test_proteoforms_by_proforma_with_fallback():
    hits = hierarchy_data.assign_groups(_hits())
    forms = hierarchy_data.proteoforms_from(hits)
    assert forms['proteoform_id'].to_list() == [1, 2, 3, 4]
    assert forms['protein_group'].to_list() == [1, 1, 2, 2]
    acdek = forms.filter(pl.col('proteoform') == 'ACDEK').row(0, named=True)
    assert acdek['prsms'] == 2 and acdek['best_index'] == 1
    assert (acdek['start'], acdek['end']) == (1, 5)
    # Undetermined termini become nulls, not 0 / -1
    mklv = forms.filter(pl.col('proteoform') == 'MKLV').row(0, named=True)
    assert mklv['start'] is None and mklv['end'] is None
    # A hit without ProForma falls back to accession and residue range
    assert forms.filter(pl.col('proteoform').str.starts_with('P2|3|9')).height == 1


def test_prsm_rows_belong_to_their_proteoform():
    hits = hierarchy_data.assign_groups(_hits())
    rows = hierarchy_data.prsm_rows_from(hits)
    forms = hierarchy_data.proteoforms_from(hits)
    acdek = forms.filter(pl.col('proteoform') == 'ACDEK')['proteoform_id'][0]
    assert sorted(rows.filter(pl.col('proteoform_id') == acdek)['Scan'].to_list()) == [10, 11]
    assert rows['qvalue'].null_count() == rows.height


def test_qvalue_cutoff_and_key_figures():
    hits = hierarchy_data.assign_groups(_hits([0.0, 0.005, 0.2, 0.01, 0.5]))
    assert sorted(hits['Scan'].to_list()) == [10, 11, 13]
    figures = hierarchy_data.key_figures(hits)
    assert figures == {'proteins': 2, 'proteoforms': 2, 'prsms': 3, 'cutoff': 'q ≤ 0.01'}
    assert hierarchy_data.level_title('Proteins', figures['proteins']) == 'Proteins (2)'


def test_without_qvalues_nothing_is_filtered():
    hits = hierarchy_data.assign_groups(_hits())
    assert hits.height == 5
    figures = hierarchy_data.key_figures(hits)
    assert figures['cutoff'] is None


def test_protein_level_qvalues_beside_proteoform_level():
    hits = _hits([0.0, 0.005, 0.002, 0.01, 0.008]).with_columns(
        pl.Series('ProteinLevelQvalue', [0.001, 0.001, 0.001, -1.0, 0.004]),
    )
    hits = hierarchy_data.assign_groups(hits)
    groups = hierarchy_data.protein_groups_from(hits)
    p1 = groups.filter(pl.col('accession') == 'P1').row(0, named=True)
    p2 = groups.filter(pl.col('accession') == 'P2').row(0, named=True)
    assert (p1['best_qvalue'], p1['best_protein_qvalue']) == (0.0, 0.001)
    # FLASHTnT's -1 (not determined) is left out of the lowest value
    assert p2['best_protein_qvalue'] == 0.004
    forms = hierarchy_data.proteoforms_from(hits)
    assert forms.filter(pl.col('proteoform') == 'MKLV')['best_protein_qvalue'][0] is None
    rows = hierarchy_data.prsm_rows_from(hits)
    assert rows.filter(pl.col('Scan') == 14)['protein_qvalue'][0] == 0.004


def test_protein_level_qvalues_missing_from_older_output():
    hits = hierarchy_data.assign_groups(_hits([0.0] * 5))
    assert hierarchy_data.protein_groups_from(hits)['best_protein_qvalue'].null_count() == 2
    assert hierarchy_data.prsm_rows_from(hits)['protein_qvalue'].null_count() == 5


def test_most_rows_per_table():
    hits = hierarchy_data.assign_groups(_hits())
    assert hierarchy_data.most_rows(hits) == {
        'protein_group_table': 2, 'proteoform_table': 2, 'proteoform_prsm_table': 2,
    }
    assert hierarchy_data.most_rows(hits.clear())['proteoform_table'] == 0


def test_proteoform_map_one_point_per_prsm():
    hits = hierarchy_data.assign_groups(_hits())
    points = hierarchy_data.proteoform_map_from(hits)
    assert points.height == hits.height
    assert set(points.columns) >= {
        'rt', 'mass', 'score', 'index', 'scan_idx', 'proteoform_id', 'protein_group', 'label',
    }
    best = points.filter(pl.col('index') == 3).row(0, named=True)
    # The point carries the ids the hierarchy tables select by
    assert (best['protein_group'], best['score'], best['scan_idx']) == (1, 50.0, 13)
    assert best['label'] == 'P2 · MKLV<br>Scan 13'


def test_proteoform_map_mass_falls_back_to_precursor_mass():
    hits = hierarchy_data.assign_groups(_hits().with_columns(
        pl.Series('ProteoformMass', [1000.0, -1.0, -1.0, 2000.0, -1.0]),
        pl.Series('PrecursorMass', [1001.0, 1101.0, -1.0, 2001.0, 0.0]),
    ))
    points = hierarchy_data.proteoform_map_from(hits)
    masses = dict(zip(points['index'].to_list(), points['mass'].to_list()))
    # FLASHTnT's proteoform mass first, else the precursor mass; neither: no point
    assert masses == {0: 1000.0, 1: 1101.0, 3: 2000.0}


def test_proteoform_map_shortens_long_proteoforms():
    long_form = 'A' * 100
    hits = hierarchy_data.assign_groups(_hits().with_columns(
        pl.Series('Proforma', [long_form] * 5)
    ))
    label = hierarchy_data.proteoform_map_from(hits)['label'][0]
    assert label.startswith('P2 · ' + 'A' * 39 + '…<br>')


def test_proteoform_map_empty():
    hits = hierarchy_data.assign_groups(_hits()).clear()
    assert hierarchy_data.proteoform_map_from(hits).height == 0
