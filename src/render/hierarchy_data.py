"""Input frames of the protein → proteoform → PrSM tables.

Identifications are browsed in three levels: a protein list (each protein
with its number of proteoforms and its best hit), the proteoforms of a protein
(each with its number of PrSMs), and the PrSMs of a proteoform. This module
builds those levels from FLASHTnT's protein table.

A PrSM is the best hit of a spectrum (`is_best`), and when FLASHTnT estimated
q-values (a target/decoy search) only hits within QVALUE_CUTOFF count. A
proteoform is identified by its ProForma string (the sequence with its
modifications); hits without one fall back to accession, residue range and
modification masses. Like insight_data, the module is free of Streamlit and
openms-insight.
"""

import polars as pl

from src.render import insight_data

QVALUE_CUTOFF = 0.01


def _with_missing(column):
    """FLASHTnT writes -1 for a value it did not determine."""
    return pl.when(pl.col(column) < 0).then(None).otherwise(pl.col(column))


def _undetermined_position(column):
    return pl.when(pl.col(column) <= 0).then(None).otherwise(pl.col(column))


def prsms(file_manager, dataset_id):
    """The PrSMs (best hit of each spectrum, within the q-value cutoff when
    q-values exist) with `protein_group` and `proteoform_id` columns."""
    hits = insight_data.proteins(file_manager, dataset_id).filter(pl.col('is_best')).collect()
    return assign_groups(hits)


def has_qvalues(hits):
    return (
        'ProteoformLevelQvalue' in hits.columns
        and hits.filter(pl.col('ProteoformLevelQvalue') >= 0).height > 0
    )


def assign_groups(hits):
    """Filter the best hits by q-value (when there are any) and number their
    proteins and proteoforms by their best score."""
    if has_qvalues(hits):
        hits = hits.filter(
            pl.col('ProteoformLevelQvalue').is_between(0, QVALUE_CUTOFF)
        )
    fallback = pl.concat_str(
        [
            pl.col('accession').fill_null(''),
            pl.col('StartPosition').cast(pl.String).fill_null(''),
            pl.col('EndPosition').cast(pl.String).fill_null(''),
            pl.col('ModMass').cast(pl.String).fill_null(''),
        ],
        separator='|',
    )
    hits = hits.with_columns(
        pl.when(pl.col('Proforma').fill_null('') != '')
        .then(pl.col('Proforma'))
        .otherwise(fallback)
        .alias('proteoform_key'),
        pl.col('accession').fill_null('').alias('accession'),
    )
    # Proteins and proteoforms ordered by their best score, numbered from 1
    protein_order = (
        hits.group_by('accession')
        .agg(pl.col('Score').max().alias('_best'))
        .sort(['_best', 'accession'], descending=[True, False], nulls_last=True)
        .with_row_index('protein_group', offset=1)
        .select(pl.col('protein_group').cast(pl.Int64), 'accession')
    )
    proteoform_order = (
        hits.group_by(['accession', 'proteoform_key'])
        .agg(pl.col('Score').max().alias('_best'))
        .join(protein_order, on='accession')
        .sort(['protein_group', '_best', 'proteoform_key'],
              descending=[False, True, False], nulls_last=True)
        .with_row_index('proteoform_id', offset=1)
        .select(pl.col('proteoform_id').cast(pl.Int64), 'accession', 'proteoform_key')
    )
    return (
        hits.join(protein_order, on='accession', how='left')
        .join(proteoform_order, on=['accession', 'proteoform_key'], how='left')
        .sort(['protein_group', 'proteoform_id', 'Score'],
              descending=[False, False, True], nulls_last=True)
    )


def _best_of():
    """Columns of the highest-scoring PrSM of each group."""
    return [
        pl.col('index').sort_by('Score', descending=True, nulls_last=True).first().alias('best_index'),
        pl.col('scan_idx').sort_by('Score', descending=True, nulls_last=True).first().alias('scan_idx'),
    ]


def protein_groups_from(hits):
    return (
        hits.group_by(['protein_group', 'accession'])
        .agg(
            pl.col('description').first(),
            pl.col('proteoform_id').n_unique().alias('proteoforms'),
            pl.len().cast(pl.Int64).alias('prsms'),
            pl.col('Score').max().alias('best_score'),
            _with_missing('ProteoformLevelQvalue').min().alias('best_qvalue'),
            *_best_of(),
        )
        .sort('protein_group')
    )


def proteoforms_from(hits):
    return (
        hits.group_by(['proteoform_id', 'protein_group', 'accession', 'proteoform_key'])
        .agg(
            # 1-based residue numbers; 0 or -1 when FLASHTnT did not determine the terminus
            _undetermined_position('StartPosition').first().alias('start'),
            _undetermined_position('EndPosition').first().alias('end'),
            _with_missing('ProteoformMass').sort_by('Score', descending=True, nulls_last=True)
            .first().alias('mass'),
            pl.col('ModCount').first().alias('modifications'),
            pl.len().cast(pl.Int64).alias('prsms'),
            pl.col('Score').max().alias('best_score'),
            _with_missing('ProteoformLevelQvalue').min().alias('best_qvalue'),
            *_best_of(),
        )
        .rename({'proteoform_key': 'proteoform'})
        .sort('proteoform_id')
    )


def prsm_rows_from(hits):
    return hits.select(
        'proteoform_id', 'index', 'scan_idx', 'Scan', 'accession',
        pl.col('RetentionTime'), 'PrecursorMass',
        _with_missing('ProteoformMass').alias('ProteoformMass'),
        'MatchingFragments', 'Score',
        _with_missing('ProteoformLevelQvalue').alias('qvalue'),
    )


# Longest ProForma string shown in a proteoform map point's hover text
_HOVER_PROTEOFORM = 40


def proteoform_map_from(hits):
    """One point per PrSM for the proteoform heatmap: retention time, mass and
    score, the ids the hierarchy tables select by, and a hover label.

    The mass is the proteoform mass FLASHTnT determined, else the precursor
    mass; PrSMs with neither (single spectra without precursor information)
    have no place on the map and are left out.
    """
    if not hits.height:
        return pl.DataFrame(schema={
            'rt': pl.Float64, 'mass': pl.Float64, 'score': pl.Float64, 'index': pl.Int64,
            'scan_idx': pl.Int64, 'proteoform_id': pl.Int64, 'protein_group': pl.Int64,
            'accession': pl.String, 'label': pl.String,
        })
    proteoform = pl.col('proteoform_key')
    shortened = (
        pl.when(proteoform.str.len_chars() > _HOVER_PROTEOFORM)
        .then(proteoform.str.slice(0, _HOVER_PROTEOFORM - 1) + '…')
        .otherwise(proteoform)
    )
    return (
        hits.with_columns(
            pl.coalesce(
                _with_missing('ProteoformMass'),
                pl.when(pl.col('PrecursorMass') > 0).then(pl.col('PrecursorMass')),
            ).cast(pl.Float64).alias('mass'),
        )
        .filter(pl.col('mass').is_not_null() & pl.col('RetentionTime').is_not_null())
        .select(
            pl.col('RetentionTime').cast(pl.Float64).alias('rt'),
            'mass',
            pl.col('Score').cast(pl.Float64).alias('score'),
            pl.col('index').cast(pl.Int64),
            pl.col('scan_idx').cast(pl.Int64),
            'proteoform_id', 'protein_group', 'accession',
            pl.concat_str(
                [
                    pl.col('accession'), pl.lit(' · '), shortened,
                    pl.lit('<br>Scan '), pl.col('Scan').cast(pl.String).fill_null('-'),
                ],
            ).alias('label'),
        )
    )


def key_figures(hits):
    """Counts of the protein list's title: proteins, proteoforms, PrSMs and
    the cutoff they were counted at."""
    return {
        'proteins': hits['accession'].n_unique() if hits.height else 0,
        'proteoforms': hits['proteoform_id'].n_unique() if hits.height else 0,
        'prsms': hits.height,
        'cutoff': f'q ≤ {QVALUE_CUTOFF:g}' if 'ProteoformLevelQvalue' in hits.columns
        and hits.filter(pl.col('ProteoformLevelQvalue') >= 0).height else None,
    }


def most_rows(hits):
    """The most rows each table can show at once: all proteins, the
    proteoforms of the protein with the most, the PrSMs of the proteoform
    with the most."""
    if not hits.height:
        return {'protein_group_table': 0, 'proteoform_table': 0, 'proteoform_prsm_table': 0}
    return {
        'protein_group_table': hits['protein_group'].n_unique(),
        'proteoform_table': hits.group_by('protein_group')
        .agg(pl.col('proteoform_id').n_unique())['proteoform_id'].max(),
        'proteoform_prsm_table': hits.group_by('proteoform_id').len()['len'].max(),
    }


def level_title(name, count):
    """Title of one level's table with its total, e.g. "Proteins (12)"."""
    return f'{name} ({count})'
