"""Input frames of the run-level QC panels.

Distributions over a whole run, to catch a failed acquisition, calibration
drift or poor deconvolution before reading identifications (MetaDraw's QC
histograms, TDAuditor's run metrics). Everything comes from results
FLASHDeconv and FLASHTnT already wrote; like insight_data, the module is free
of Streamlit and openms-insight.

Each function returns a lazy frame of (series, x, y) points for a SeriesPlot.
"""

import numpy as np
import polars as pl

from src.render import insight_data, validation_data

POINT_SCHEMA = {'series': pl.String, 'x': pl.Float64, 'y': pl.Float64}
HISTOGRAM_BINS = 50


def _points(rows):
    return pl.DataFrame(rows, schema=POINT_SCHEMA, orient='row').lazy()


def histogram(values_by_series, bins=HISTOGRAM_BINS, value_range=None, integer=False):
    """Counts per bin (bin centres as x) of every series, on shared bins.

    `integer` counts every whole value on its own (charges)."""
    present = {k: np.asarray(v, dtype=float) for k, v in values_by_series.items()}
    present = {k: v[np.isfinite(v)] for k, v in present.items()}
    present = {k: v for k, v in present.items() if v.size}
    if not present:
        return _points([])
    every = np.concatenate(list(present.values()))
    lo, hi = value_range if value_range else (every.min(), every.max())
    if integer:
        edges = np.arange(np.floor(lo) - 0.5, np.ceil(hi) + 1.5)
    else:
        if hi <= lo:
            hi = lo + 1
        edges = np.linspace(lo, hi, bins + 1)
    centres = (edges[:-1] + edges[1:]) / 2
    rows = []
    for series, values in present.items():
        counts, _ = np.histogram(values, bins=edges)
        rows.extend((series, float(c), float(n)) for c, n in zip(centres, counts))
    return _points(rows)


def _masses(file_manager, dataset_id):
    """Every deconvolved mass with its spectrum's MS level and retention time."""
    scans = insight_data.scan_table(file_manager, dataset_id).select(
        pl.col('index').alias('scan_idx'), 'MSLevel', 'RT'
    )
    return (
        insight_data.mass_table(file_manager, dataset_id)
        .join(scans, on='scan_idx')
        .with_columns(pl.format('MS{}', pl.col('MSLevel')).alias('series'))
        .collect()
    )


def _by_level(masses, column):
    return {
        series: group[column].to_numpy()
        for (series,), group in masses.group_by(['series'], maintain_order=True)
    }


def qscore_distribution(file_manager, dataset_id):
    masses = _masses(file_manager, dataset_id).sort('MSLevel')
    return histogram(_by_level(masses, 'QScore'), value_range=(0.0, 1.0))


def mass_distribution(file_manager, dataset_id):
    masses = _masses(file_manager, dataset_id).sort('MSLevel')
    return histogram(_by_level(masses, 'MonoMass'))


def charge_distribution(file_manager, dataset_id):
    """Masses per charge state: a mass counts for every charge from its
    lowest to its highest."""
    masses = _masses(file_manager, dataset_id).sort('MSLevel')
    charges = masses.select(
        'series',
        pl.int_ranges(pl.col('MinCharges'), pl.col('MaxCharges') + 1).alias('charge'),
    ).explode('charge').drop_nulls('charge')
    return histogram(_by_level(charges, 'charge'), integer=True)


def masses_per_spectrum(file_manager, dataset_id):
    """Number of deconvolved masses of every spectrum against its retention
    time, per MS level."""
    scans = insight_data.scan_table(file_manager, dataset_id).collect().sort('RT')
    return (
        scans.select(
            pl.format('MS{}', pl.col('MSLevel')).alias('series'),
            pl.col('RT').cast(pl.Float64).alias('x'),
            pl.col('#Masses').cast(pl.Float64).alias('y'),
        )
        .sort(['series', 'x'])
        .lazy()
    )


def identifications_over_time(file_manager, dataset_id, q_threshold=0.01, window=60.0):
    """Proteoform-spectrum matches per `window` seconds of retention time: the
    best hit of every spectrum, and those within the q-value threshold (when
    FLASHTnT estimated q-values)."""
    proteins = insight_data.proteins(file_manager, dataset_id).filter(pl.col('is_best'))
    scans = insight_data.scan_table(file_manager, dataset_id).select(
        pl.col('index').alias('scan_idx'), 'RT'
    )
    hits = proteins.join(scans, on='scan_idx').select('RT', 'ProteoformLevelQvalue').collect()
    if hits.is_empty():
        return _points([])
    rows = []
    bins = (hits['RT'] // window).cast(pl.Int64)
    for series, mask in (
        ('Best hit per spectrum', np.ones(hits.height, dtype=bool)),
        (f'q ≤ {q_threshold:g}', (
            (hits['ProteoformLevelQvalue'] >= 0) & (hits['ProteoformLevelQvalue'] <= q_threshold)
        ).fill_null(False).to_numpy()),
    ):
        if series != 'Best hit per spectrum' and not mask.any():
            continue
        counts = pl.DataFrame({'bin': bins.to_numpy()[mask]}).group_by('bin').len().sort('bin')
        rows.extend(
            (series, (b + 0.5) * window, float(n))
            for b, n in counts.iter_rows()
        )
    return _points(rows)


def fragment_error_over_time(file_manager, dataset_id):
    """Median fragment mass error (ppm) of the best hit of every spectrum
    against its retention time: a trend shows calibration drift."""
    errors = validation_data.fragment_errors(file_manager, dataset_id)
    proteins = insight_data.proteins(file_manager, dataset_id).filter(pl.col('is_best'))
    scans = insight_data.scan_table(file_manager, dataset_id).select(
        pl.col('index').alias('scan_idx'), 'RT'
    )
    return (
        errors.group_by('proteoform_index')
        .agg(pl.col('ppm_error').median().alias('y'))
        .join(proteins.select(pl.col('index').alias('proteoform_index'), 'scan_idx'),
              on='proteoform_index')
        .join(scans, on='scan_idx')
        .select(
            pl.lit('Median fragment error').alias('series'),
            pl.col('RT').cast(pl.Float64).alias('x'),
            pl.col('y').cast(pl.Float64),
        )
        .sort('x')
    )


CHROMATOGRAM_SERIES = {
    'tic': 'TIC (raw)',
    'bpc': 'Base peak (raw)',
    'deconv': 'Deconvolved TIC',
}


def chromatogram_level(file_manager, dataset_id):
    """Lowest MS level of a run's spectra (1 unless it has MS2 only)."""
    level = insight_data.scan_table(file_manager, dataset_id).select(
        pl.col('MSLevel').min()
    ).collect().item()
    return 1 if level is None else int(level)


def chromatograms(file_manager, dataset_id, ms_level=None):
    """Total ion, base peak and deconvolved total intensity of every spectrum
    of one MS level against its retention time. Without
    `ms_level`, the lowest level the run has: MS1, or MS2 for runs of
    fragmentation spectra only.

    One row per (series, spectrum) with the spectrum `index`, so a click on a
    point selects that spectrum. The deconvolved series is the summed
    intensity of the spectrum's deconvolved masses."""
    if ms_level is None:
        ms_level = chromatogram_level(file_manager, dataset_id)
    scans = insight_data.scan_table(file_manager, dataset_id).filter(pl.col('MSLevel') == ms_level).select('index', pl.col('RT').cast(pl.Float64).alias('x'))
    raw = insight_data.spectra(file_manager, dataset_id, 'anno_spectrum').select(
        'index',
        pl.col('SumIntensity_Anno').list.sum().cast(pl.Float64).alias('tic'),
        pl.col('SumIntensity_Anno').list.max().cast(pl.Float64).alias('bpc'),
    )
    frame = scans.join(raw, on='index', how='left')
    value_columns = ['tic', 'bpc']
    if insight_data.result_path(file_manager, dataset_id, 'deconv_spectrum') is not None:
        deconv = insight_data.spectra(file_manager, dataset_id, 'deconv_spectrum').select(
            'index', pl.col('SumIntensity').list.sum().cast(pl.Float64).alias('deconv'),
        )
        frame = frame.join(deconv, on='index', how='left')
        value_columns.append('deconv')
    order = {name: rank for rank, name in enumerate(CHROMATOGRAM_SERIES)}
    return (
        frame.unpivot(index=['index', 'x'], on=value_columns, variable_name='key', value_name='y')
        .with_columns(
            pl.col('y').fill_null(0.0),
            pl.col('key').replace_strict(order, return_dtype=pl.Int64).alias('rank'),
            pl.col('key').replace_strict(CHROMATOGRAM_SERIES, return_dtype=pl.String).alias('series'),
        )
        .sort(['rank', 'x'])
        .select('series', 'x', 'y', 'index')
    )
