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
        'proteoform_id', 'index', 'scan_idx', 'Scan',
        pl.col('RetentionTime'), 'PrecursorMass',
        _with_missing('ProteoformMass').alias('ProteoformMass'),
        'MatchingFragments', 'Score',
        _with_missing('ProteoformLevelQvalue').alias('qvalue'),
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


def protein_list_title(figures):
    def count(n, word):
        return f"{n} {word}{'' if n == 1 else 's'}"
    title = (
        f"Proteins · {count(figures['proteins'], 'protein')}, "
        f"{count(figures['proteoforms'], 'proteoform')}, {count(figures['prsms'], 'PrSM')}"
    )
    if figures['cutoff']:
        title += f" at {figures['cutoff']}"
    return title
