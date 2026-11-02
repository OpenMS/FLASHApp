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
    assert hierarchy_data.protein_list_title(figures) == (
        'Proteins · 2 proteins, 2 proteoforms, 3 PrSMs at q ≤ 0.01'
    )


def test_without_qvalues_nothing_is_filtered():
    hits = hierarchy_data.assign_groups(_hits())
    assert hits.height == 5
    figures = hierarchy_data.key_figures(hits)
    assert figures['cutoff'] is None
    assert hierarchy_data.protein_list_title(
        {'proteins': 1, 'proteoforms': 1, 'prsms': 1, 'cutoff': None}
    ) == 'Proteins · 1 protein, 1 proteoform, 1 PrSM'
