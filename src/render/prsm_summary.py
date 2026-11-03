"""Input frame of the PrSM summary: the key facts of one identification.

For every hit of FLASHTnT's protein table, one row per field (protein,
proteoform in ProForma notation, residue range, scan, retention time,
precursor m/z, charge and mass, proteoform mass and its difference to the
precursor, matched fragments, coverage, tags, modifications with their
localization, score, q-value), in a fixed order. The summary cell filters
these rows by the selected hit.

Precursor m/z and charge come from FLASHDeconv's MS2 spectrum table
(spec2.tsv); without it those fields show "-". Like insight_data, the module
is free of Streamlit and openms-insight.
"""

import polars as pl

from src.render import insight_data

MISSING = '-'

SCHEMA = {'row': pl.Int64, 'index': pl.Int64, 'order': pl.Int64, 'field': pl.String, 'value': pl.String}

# Columns of FLASHDeconv's spectrum table used here
_PRECURSOR_COLUMNS = {
    'ScanNum': pl.Int64,
    'PrecursorMz': pl.Float64,
    'PrecursorCharge': pl.Int64,
    'PrecursorMonoisotopicMass': pl.Float64,
}


def _known(column):
    """FLASHTnT and FLASHDeconv write -1, 0 or NaN for a value they did not determine."""
    value = pl.col(column)
    return pl.when(value.is_not_nan() & (value > 0)).then(value).otherwise(None)


def _fixed(expr, decimals):
    return pl.when(expr.is_null()).then(pl.lit(MISSING)).otherwise(expr.round(decimals).cast(pl.String))


def _text(expr):
    return pl.when(expr.is_null() | (expr.cast(pl.String) == '')).then(pl.lit(MISSING)).otherwise(expr.cast(pl.String))


def precursors(spectrum_table):
    """One row per MS2 scan: precursor m/z, charge and monoisotopic mass,
    from a frame with FLASHDeconv's spec2.tsv columns."""
    present = [c for c in _PRECURSOR_COLUMNS if c in spectrum_table.columns]
    if 'ScanNum' not in present:
        return pl.DataFrame(schema={'Scan': pl.Int64, **{c: t for c, t in _PRECURSOR_COLUMNS.items() if c != 'ScanNum'}})
    frame = spectrum_table.select(
        [pl.col(c).cast(_PRECURSOR_COLUMNS[c], strict=False) for c in present]
    )
    for column, dtype in _PRECURSOR_COLUMNS.items():
        if column not in frame.columns:
            frame = frame.with_columns(pl.lit(None, dtype=dtype).alias(column))
    return (
        frame.rename({'ScanNum': 'Scan'})
        .group_by('Scan', maintain_order=True)
        .agg(pl.all().drop_nulls().first())
    )


def read_precursors(path):
    """precursors() of a spec2.tsv file, or None when there is none."""
    if path is None:
        return None
    table = pl.read_csv(
        path, separator='\t', columns=list(_PRECURSOR_COLUMNS), infer_schema_length=0,
        null_values=['nan', 'NaN', ''],
    )
    return precursors(table)


def _modifications(hits):
    """Per hit: 'n (k localized)' and one line per modification:
    '<mass> Da at <residue>' or '<mass> Da in <start>–<end> (not localized)'."""
    split = (
        hits.select(
            'index',
            pl.col('ModMass').cast(pl.String).fill_null('').str.split(';').alias('mass'),
            pl.col('ModStart').cast(pl.String).fill_null('').str.split(';').alias('start'),
            pl.col('ModEnd').cast(pl.String).fill_null('').str.split(';').alias('end'),
        )
        .explode(['mass', 'start', 'end'])
        .filter(pl.col('mass') != '')
        .with_columns(
            pl.col('mass').cast(pl.Float64, strict=False),
            # Residue positions in the protein, 1-based as FLASHTnT writes them
            pl.col('start').cast(pl.Int64, strict=False),
            pl.col('end').cast(pl.Int64, strict=False),
        )
        .with_columns((pl.col('start') == pl.col('end')).alias('localized'))
        .with_columns(
            pl.format(
                '{} Da {}',
                pl.when(pl.col('mass') >= 0).then(pl.lit('+')).otherwise(pl.lit('−'))
                + pl.col('mass').abs().round(4).cast(pl.String),
                pl.when(pl.col('localized'))
                .then(pl.format('at {}', pl.col('start')))
                .otherwise(pl.format('in {}–{} (not localized)', pl.col('start'), pl.col('end'))),
            ).alias('line')
        )
    )
    return split.group_by('index', maintain_order=True).agg(
        pl.len().alias('mods'),
        pl.col('localized').sum().alias('localized'),
        pl.col('line'),
    )


def summary_rows(hits, precursor_table=None):
    """Long frame (row, index, order, field, value) of the hits of the
    protein table; `index` is the hit, `row` a unique row id."""
    if hits.is_empty():
        return pl.DataFrame(schema=SCHEMA)
    hits = hits.with_columns(pl.col('index').cast(pl.Int64), pl.col('Scan').cast(pl.Int64))
    if precursor_table is not None and not precursor_table.is_empty():
        hits = hits.join(precursor_table, on='Scan', how='left')
    else:
        hits = hits.with_columns(
            pl.lit(None, dtype=pl.Float64).alias('PrecursorMz'),
            pl.lit(None, dtype=pl.Int64).alias('PrecursorCharge'),
            pl.lit(None, dtype=pl.Float64).alias('PrecursorMonoisotopicMass'),
        )
    hits = hits.join(_modifications(hits), on='index', how='left').with_columns(
        pl.col('mods').fill_null(0), pl.col('localized').fill_null(0)
    )
    return pl.concat([
        hits.select(
            'index',
            pl.lit(order, dtype=pl.Int64).alias('order'),
            pl.lit(name).alias('field'),
            expr.cast(pl.String).alias('value'),
        )
        for order, (name, expr) in enumerate(_fields())
    ]).sort(['index', 'order']).select(
        (pl.col('index') * FIELD_SLOTS + pl.col('order')).alias('row'), 'index', 'order', 'field', 'value',
    )


def _fields():
    """(field name, value expression) of the summary, in display order."""
    precursor_mass = pl.coalesce(_known('PrecursorMass'), _known('PrecursorMonoisotopicMass'))
    proteoform_mass = _known('ProteoformMass')
    start = pl.when(pl.col('StartPosition') > 0).then(pl.col('StartPosition').cast(pl.String)).otherwise(pl.lit('?'))
    end = pl.when(pl.col('EndPosition') > 0).then(pl.col('EndPosition').cast(pl.String)).otherwise(pl.lit('?'))
    qvalue = pl.when(pl.col('ProteoformLevelQvalue') >= 0).then(pl.col('ProteoformLevelQvalue')).otherwise(None)
    description = pl.col('description').cast(pl.String).fill_null('')
    return [
        ('Protein', pl.when(description == '').then(pl.col('accession').cast(pl.String))
            .otherwise(pl.format('{} ({})', pl.col('accession'), description))),
        ('Proteoform', _text(pl.col('Proforma'))),
        # "?" marks a terminus FLASHTnT did not determine
        ('Residues', pl.format('{}–{} of {}', start, end, pl.col('length'))),
        ('Scan', _text(pl.col('Scan'))),
        ('Retention time (s)', _fixed(_known('RetentionTime'), 2)),
        ('Precursor m/z', _fixed(_known('PrecursorMz'), 4)),
        ('Precursor charge', _text(pl.when(pl.col('PrecursorCharge') > 0).then(pl.col('PrecursorCharge')))),
        ('Precursor mass (Da)', _fixed(precursor_mass, 4)),
        ('Proteoform mass (Da)', _fixed(proteoform_mass, 4)),
        ('Δ mass (Da)', _fixed(precursor_mass - proteoform_mass, 4)),
        ('Matched fragments', _text(pl.col('MatchingFragments'))),
        ('Sequence coverage (%)', _fixed(pl.col('Coverage(%)').cast(pl.Float64), 1)),
        ('Tags', _text(pl.col('TagCount'))),
        ('Modifications', pl.when(pl.col('mods') == 0).then(pl.lit('none'))
            .otherwise(pl.format('{} ({} localized)', pl.col('mods'), pl.col('localized')))),
        ('Modification sites', pl.when(pl.col('mods') == 0).then(pl.lit(MISSING))
            .otherwise(pl.col('line').list.join('; '))),
        ('Score', _text(pl.col('Score'))),
        ('q-value (proteoform level)', _fixed(qvalue, 4)),
    ]


FIELD_COUNT = len(_fields())
# Room for the field numbers of one hit in the row id
FIELD_SLOTS = 100


def prsm_summary(file_manager, dataset_id):
    """summary_rows() of a dataset, as a lazy frame."""
    hits = insight_data.proteins(file_manager, dataset_id).collect()
    path = insight_data.result_path(file_manager, dataset_id, 'spec2_tsv')
    return summary_rows(hits, read_precursors(path)).lazy()
