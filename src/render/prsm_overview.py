"""Overview columns of the protein table.

The protein table shows a few essential columns; a button in each row opens
the whole row as an overview (openms-insight `Table(row_overview=True)`).
This module adds the columns that only the overview shows: the residue range,
precursor m/z and charge, the precursor mass, the difference between precursor
and proteoform mass, and where each modification sits.

Precursor m/z and charge come from FLASHDeconv's MS2 spectrum table
(spec2.tsv); without it they stay empty. Like insight_data, the module is free
of Streamlit and openms-insight.
"""

import polars as pl

from src.render import insight_data

# Columns of FLASHDeconv's spectrum table used here
_PRECURSOR_COLUMNS = {
    'ScanNum': pl.Int64,
    'PrecursorMz': pl.Float64,
    'PrecursorCharge': pl.Int64,
    'PrecursorMonoisotopicMass': pl.Float64,
}

OVERVIEW_SCHEMA = {
    'index': pl.Int64,
    'residues': pl.String,
    'precursor_mz': pl.Float64,
    'precursor_charge': pl.Int64,
    'precursor_mass': pl.Float64,
    'delta_mass': pl.Float64,
    'modification_sites': pl.String,
}


def _known(column):
    """FLASHTnT and FLASHDeconv write -1, 0 or NaN for a value they did not determine."""
    value = pl.col(column).cast(pl.Float64)
    return pl.when(value.is_not_nan() & (value > 0)).then(value).otherwise(None)


def precursors(spectrum_table):
    """One row per MS2 scan (`Scan`): precursor m/z, charge and monoisotopic
    mass, from a frame with FLASHDeconv's spec2.tsv columns."""
    frame = spectrum_table.select(
        [pl.col(c).cast(t, strict=False) for c, t in _PRECURSOR_COLUMNS.items() if c in spectrum_table.columns]
    )
    for column, dtype in _PRECURSOR_COLUMNS.items():
        if column not in frame.columns:
            frame = frame.with_columns(pl.lit(None, dtype=dtype).alias(column))
    return (
        frame.rename({'ScanNum': 'Scan'})
        .drop_nulls('Scan')
        .group_by('Scan', maintain_order=True)
        .agg(pl.all().drop_nulls().first())
    )


def read_precursors(path):
    """precursors() of a spec2.tsv file, or None without one."""
    if path is None:
        return None
    header = pl.read_csv(path, separator='\t', n_rows=0, infer_schema_length=0).columns
    table = pl.read_csv(
        path, separator='\t', columns=[c for c in _PRECURSOR_COLUMNS if c in header],
        infer_schema_length=0, null_values=['nan', 'NaN', ''],
    )
    return precursors(table)


def modification_sites(hits):
    """Per hit (`index`): one entry per modification, '<mass> Da at <residue>'
    or '<mass> Da in <first>–<last> (not localized)', joined by '; '.
    Positions are FLASHTnT's 1-based residue numbers in the protein."""
    sites = (
        hits.select(
            pl.col('index').cast(pl.Int64),
            pl.col('ModMass').cast(pl.String).fill_null('').str.split(';').alias('mass'),
            pl.col('ModStart').cast(pl.String).fill_null('').str.split(';').alias('start'),
            pl.col('ModEnd').cast(pl.String).fill_null('').str.split(';').alias('end'),
        )
        .explode(['mass', 'start', 'end'])
        .filter(pl.col('mass') != '')
        .with_columns(
            pl.col('mass').cast(pl.Float64, strict=False),
            pl.col('start').cast(pl.Int64, strict=False),
            pl.col('end').cast(pl.Int64, strict=False),
        )
        .with_columns(
            pl.format(
                '{} Da {}',
                pl.when(pl.col('mass') >= 0).then(pl.lit('+')).otherwise(pl.lit('−'))
                + pl.col('mass').abs().round(4).cast(pl.String),
                pl.when(pl.col('start') == pl.col('end'))
                .then(pl.format('at {}', pl.col('start')))
                .otherwise(pl.format('in {}–{} (not localized)', pl.col('start'), pl.col('end'))),
            ).alias('site')
        )
    )
    return sites.group_by('index', maintain_order=True).agg(
        pl.col('site').str.join('; ').alias('modification_sites')
    )


def overview_columns(hits, precursor_table=None):
    """Frame keyed by `index` with the columns of OVERVIEW_SCHEMA."""
    if hits.is_empty():
        return pl.DataFrame(schema=OVERVIEW_SCHEMA)
    hits = hits.with_columns(pl.col('index').cast(pl.Int64), pl.col('Scan').cast(pl.Int64))
    if precursor_table is not None and not precursor_table.is_empty():
        hits = hits.join(precursor_table, on='Scan', how='left')
    else:
        hits = hits.with_columns(
            pl.lit(None, dtype=pl.Float64).alias('PrecursorMz'),
            pl.lit(None, dtype=pl.Int64).alias('PrecursorCharge'),
            pl.lit(None, dtype=pl.Float64).alias('PrecursorMonoisotopicMass'),
        )
    # "?" marks a terminus FLASHTnT did not determine
    start = pl.when(pl.col('StartPosition') > 0).then(pl.col('StartPosition').cast(pl.String)).otherwise(pl.lit('?'))
    end = pl.when(pl.col('EndPosition') > 0).then(pl.col('EndPosition').cast(pl.String)).otherwise(pl.lit('?'))
    precursor_mass = pl.coalesce(_known('PrecursorMass'), _known('PrecursorMonoisotopicMass'))
    frame = hits.select(
        'index',
        pl.format('{}–{} of {}', start, end, pl.col('length')).alias('residues'),
        _known('PrecursorMz').alias('precursor_mz'),
        pl.when(pl.col('PrecursorCharge') > 0).then(pl.col('PrecursorCharge')).cast(pl.Int64).alias('precursor_charge'),
        precursor_mass.alias('precursor_mass'),
        (precursor_mass - _known('ProteoformMass')).alias('delta_mass'),
    )
    return frame.join(modification_sites(hits), on='index', how='left').select(list(OVERVIEW_SCHEMA))


def protein_overview(file_manager, dataset_id):
    """overview_columns() of a dataset, as a lazy frame."""
    hits = insight_data.proteins(file_manager, dataset_id).collect()
    path = insight_data.result_path(file_manager, dataset_id, 'spec2_tsv')
    return overview_columns(hits, read_precursors(path)).lazy()
