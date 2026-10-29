"""Input frames for the openms-insight viewer components.

Every function turns the per-dataset caches the workflows write (see
src/parse/) into the shape one Insight component expects. They return lazy
frames: the Insight component collects them only when it builds its own cache,
which happens once per dataset (see src/render/render.py).

Kept free of Streamlit and openms-insight so it can be tested on its own.
"""

import polars as pl

from src.render.sequence_data_store import build_table


def result_path(file_manager, dataset_id, name_tag):
    """Path of a stored result, or None when it was not written."""
    path = file_manager.get_result_path(dataset_id, name_tag)
    if path is None or not path.exists():
        return None
    return path


def scan(file_manager, dataset_id, name_tag):
    """Lazy frame of a stored parquet result."""
    return pl.scan_parquet(result_path(file_manager, dataset_id, name_tag))


def scan_table(file_manager, dataset_id):
    """One row per spectrum; `index` is the spectrum index all selections use."""
    return scan(file_manager, dataset_id, 'scan_table').with_columns(
        pl.col('index').cast(pl.Int64)
    )


def mass_table(file_manager, dataset_id):
    """One row per deconvolved mass: scan_idx, mass_idx and the mass columns."""
    lf = scan(file_manager, dataset_id, 'mass_table')
    list_columns = [
        c for c, t in lf.collect_schema().items() if isinstance(t, pl.List)
    ]
    return (
        lf.rename({'index': 'scan_idx'})
        .with_columns(
            pl.col('scan_idx').cast(pl.Int64),
            pl.int_ranges(0, pl.col('MonoMass').list.len()).alias('mass_idx'),
        )
        .explode(list_columns + ['mass_idx'])
        .drop_nulls('mass_idx')
        .with_columns(pl.col('mass_idx').cast(pl.Int64))
        .select(['scan_idx', 'mass_idx'] + list_columns)
    )


def spectra(file_manager, dataset_id, name_tag):
    """One row per spectrum with list columns (deconv/combined spectrum caches)."""
    return scan(file_manager, dataset_id, name_tag).with_columns(
        pl.col('index').cast(pl.Int64)
    )


def _peaks_per_mass(lf, column, kind):
    """Rows of (scan_idx, mass_idx, kind, mass, charge, intensity lists) for the
    signal or noise peaks of every mass. A peak is [index, m/z, intensity,
    charge]; its mass is m/z times charge, as FLASHViewer draws it."""
    return (
        lf.select(pl.col('index').cast(pl.Int64).alias('scan_idx'), pl.col(column).alias('peaks'))
        .with_columns(pl.int_ranges(0, pl.col('peaks').list.len()).alias('mass_idx'))
        .explode(['peaks', 'mass_idx'])
        .drop_nulls('mass_idx')
        .select(
            'scan_idx',
            pl.col('mass_idx').cast(pl.Int64),
            pl.lit(kind).alias('kind'),
            pl.col('peaks').list.eval(
                pl.element().list.get(1, null_on_oob=True)
                * pl.element().list.get(3, null_on_oob=True)
            ).alias('mass'),
            pl.col('peaks').list.eval(pl.element().list.get(3, null_on_oob=True)).alias('charge'),
            pl.col('peaks').list.eval(pl.element().list.get(2, null_on_oob=True)).alias('intensity'),
        )
    )


def signal_peaks(file_manager, dataset_id):
    """Signal and noise peaks per mass for the 3D plot.

    Rows with mass_idx -1 hold, for a fragmentation spectrum, the peaks of its
    precursor: the mass of the precursor spectrum within 0.01 Da of the
    precursor mass. They are shown while no mass is selected.
    """
    threedim = scan(file_manager, dataset_id, 'threedim_SN_plot')
    peaks = pl.concat([
        _peaks_per_mass(threedim, 'SignalPeaks', 'Signal'),
        _peaks_per_mass(threedim, 'NoisyPeaks', 'Noise'),
    ])

    scans = scan_table(file_manager, dataset_id)
    masses = (
        spectra(file_manager, dataset_id, 'deconv_spectrum')
        .select(
            pl.col('index').alias('precursor_idx'),
            pl.col('MonoMass'),
            pl.int_ranges(0, pl.col('MonoMass').list.len()).alias('precursor_mass_idx'),
        )
        .explode(['MonoMass', 'precursor_mass_idx'])
        .drop_nulls('precursor_mass_idx')
    )
    precursors = (
        threedim.select(
            pl.col('index').cast(pl.Int64).alias('scan_idx'),
            pl.col('PrecursorScan').cast(pl.Int64),
        )
        .filter(pl.col('PrecursorScan') > 0)
        .join(scans.select(pl.col('index').alias('scan_idx'), 'PrecursorMass'), on='scan_idx')
        .join(
            scans.select(pl.col('Scan').alias('PrecursorScan'), pl.col('index').alias('precursor_idx'))
            .unique('PrecursorScan', keep='first'),
            on='PrecursorScan',
        )
        .join(masses, on='precursor_idx')
        .filter((pl.col('MonoMass') - pl.col('PrecursorMass')).abs() < 1e-2)
        .group_by('scan_idx')
        .agg(pl.col('precursor_idx').first(), pl.col('precursor_mass_idx').min())
    )
    precursor_peaks = (
        precursors.join(
            peaks.rename({'scan_idx': 'precursor_idx', 'mass_idx': 'precursor_mass_idx'}),
            on=['precursor_idx', 'precursor_mass_idx'],
        )
        .select(
            'scan_idx', pl.lit(-1, dtype=pl.Int64).alias('mass_idx'),
            'kind', 'mass', 'charge', 'intensity',
        )
    )
    return pl.concat([peaks, precursor_peaks])


def fdr_densities(file_manager, dataset_id, target_tag, decoy_tag):
    """Target and decoy QScore densities in one frame, `kind` telling them apart."""
    frames = []
    for name_tag, kind in ((target_tag, 'target'), (decoy_tag, 'decoy')):
        if result_path(file_manager, dataset_id, name_tag) is None:
            continue
        frames.append(
            scan(file_manager, dataset_id, name_tag).select(
                pl.col('x').cast(pl.Float64),
                pl.col('y').cast(pl.Float64),
                pl.lit(kind).alias('kind'),
            )
        )
    if not frames:
        return pl.LazyFrame(
            schema={'x': pl.Float64, 'y': pl.Float64, 'kind': pl.String}
        )
    return pl.concat(frames)


def _scan_index_of(file_manager, dataset_id, lf):
    """Add scan_idx (the spectrum index) to a frame with a Scan number column."""
    scans = (
        scan_table(file_manager, dataset_id)
        .select(pl.col('Scan').cast(pl.Int64), pl.col('index').alias('scan_idx'))
        .unique('Scan', keep='first')
    )
    return lf.with_columns(pl.col('Scan').cast(pl.Int64)).join(scans, on='Scan', how='left', maintain_order='left')


def proteins(file_manager, dataset_id):
    """Proteoform table with scan_idx and is_best (the highest-scoring hit of
    its spectrum; the first one on ties, rows without a scan always count)."""
    lf = scan(file_manager, dataset_id, 'protein_dfs')
    lf = lf.with_columns(pl.col('index').cast(pl.Int64)).with_row_index('_order')
    best = (
        lf.filter(pl.col('Scan').is_not_null())
        .sort(['Scan', 'Score', '_order'], descending=[False, True, False], nulls_last=True)
        .group_by('Scan', maintain_order=True)
        .agg(pl.col('_order').first().alias('_best'))
    )
    lf = (
        lf.join(best, on='Scan', how='left', maintain_order='left')
        .with_columns(
            (pl.col('_best').is_null() | (pl.col('_order') == pl.col('_best'))).alias('is_best')
        )
        .drop(['_order', '_best'])
    )
    return _scan_index_of(file_manager, dataset_id, lf)


def tags(file_manager, dataset_id):
    """Sequence tags with scan_idx (tags belong to spectra, not proteoforms)."""
    lf = scan(file_manager, dataset_id, 'tag_dfs')
    return _scan_index_of(file_manager, dataset_id, lf)


def proteoform_sequences(file_manager, dataset_id):
    """One row per proteoform: its sequence data (see sequence_data_store) and
    the deconvolved masses and precursor mass of its spectrum."""
    path = result_path(file_manager, dataset_id, 'sequence_data')
    if path.suffix == '.pq':
        sequences = pl.scan_parquet(path)
    else:
        # Results written before the parquet store hold a pickled dict
        sequences = pl.from_arrow(
            build_table(file_manager.get_results(dataset_id, ['sequence_data'])['sequence_data'])
        ).lazy()
    proteoform_scans = proteins(file_manager, dataset_id).select(
        pl.col('index').alias('proteoform_index'), 'scan_idx'
    )
    spectrum = spectra(file_manager, dataset_id, 'sequence_view').select(
        pl.col('index').alias('scan_idx'), 'MonoMass', 'PrecursorMass'
    )
    return (
        sequences.with_columns(pl.col('proteoform_index').cast(pl.Int64))
        .join(proteoform_scans, on='proteoform_index', how='left', maintain_order='left')
        .join(spectrum, on='scan_idx', how='left', maintain_order='left')
    )


def _split_numbers(column):
    return pl.col(column).str.split(',').list.eval(
        pl.element().str.strip_chars().cast(pl.Float64, strict=False)
    )


def quant_traces(file_manager, dataset_id):
    """Mass traces of every feature group: one row per trace with list columns
    of m/z, retention time and intensity."""
    lf = scan(file_manager, dataset_id, 'quant_dfs').select(
        'FeatureGroupIndex', 'Charges', 'MZs', 'RTs', 'Intensities'
    )
    return (
        lf.with_columns(pl.int_ranges(0, pl.col('Charges').list.len()).alias('trace'))
        .explode(['Charges', 'MZs', 'RTs', 'Intensities', 'trace'])
        .drop_nulls('trace')
        .select(
            pl.col('FeatureGroupIndex').cast(pl.Int64),
            pl.col('trace').cast(pl.Int64),
            pl.col('Charges').cast(pl.Int64).alias('charge'),
            _split_numbers('MZs').alias('mz'),
            _split_numbers('RTs').alias('rt'),
            _split_numbers('Intensities').alias('intensity'),
        )
        # Legend in charge order, as FLASHViewer lists it
        .sort(['FeatureGroupIndex', 'charge', 'trace'])
    )


def quant_table(file_manager, dataset_id):
    """Feature groups (list columns dropped: the table shows scalars)."""
    lf = scan(file_manager, dataset_id, 'quant_dfs')
    return lf.select(
        [c for c, t in lf.collect_schema().items() if not isinstance(t, pl.List)]
    )


def sequence_view_spectra(file_manager, dataset_id):
    """Deconvolved masses and precursor mass per spectrum, for matching a
    user-given sequence (FLASHDeconv)."""
    return spectra(file_manager, dataset_id, 'sequence_view')
