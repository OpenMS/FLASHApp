"""What a comparison of runs was made from: the stored runs behind it, the
PrSMs of every run and the reporter ions of every labelled spectrum.

The identification and quality-control views of a comparison
(src/render/quant_cells.py) follow a proteoform to its PrSMs in each run and
into the run's own sequence view, and show the reporter ions and the
consensus features themselves. Streamlit-free.

A comparison names its runs by the base name of their mzML (``runs`` of the
manifest; for isobaric labels the labelled samples, and ``fractions`` the
runs). ``sources.runs`` holds the dataset ids the runs were stored under.
A run whose dataset is gone, or that was stored under another id (the
example results), is found by the mzML its analysis settings name.
"""

import re
from pathlib import Path

import polars as pl

from src.quant import contract

SETTINGS_KEY = 'analysis_settings'
# The time stamp a run's dataset id ends with (<file>_<YYYYmmdd-HHMMSS>)
_STAMP = re.compile(r'_\d{8}-\d{6}$')


def _stem(path):
    return Path(str(path)).stem if path else None


def _run_of_dataset(file_manager, dataset_id):
    """The run name (mzML base name) a stored run was made from, or None."""
    try:
        if not file_manager.result_exists(dataset_id, SETTINGS_KEY):
            return None
        record = file_manager.get_results(dataset_id, [SETTINGS_KEY])[SETTINGS_KEY]
    except Exception:  # a damaged or unknown dataset
        return None
    return _stem(record.get('mzml')) if isinstance(record, dict) else None


def run_datasets(file_manager, manifest):
    """Run name -> dataset id of the stored runs a comparison was made from
    (the fractions of an isobaric comparison). Runs without a stored dataset
    are left out."""
    sources = list((manifest.get('sources') or {}).get('runs') or [])
    names = list(manifest.get('fractions') or []) if manifest.get('kind') == 'isobaric' else list(manifest.get('runs') or [])
    if manifest.get('kind') == 'isobaric':
        # The fractions are dataset ids: their run names come from them
        names = [_STAMP.sub('', f) for f in names]
    found = {}
    for dataset_id in sources:
        if file_manager.result_exists(dataset_id, SETTINGS_KEY):
            name = _run_of_dataset(file_manager, dataset_id) or _STAMP.sub('', dataset_id)
            found.setdefault(name, dataset_id)
    missing = [n for n in names if n not in found]
    if missing:
        try:
            stored = file_manager.get_results_list([SETTINGS_KEY])
        except Exception:
            stored = []
        for dataset_id in sorted(stored, reverse=True):
            name = _run_of_dataset(file_manager, dataset_id)
            if name in missing and name not in found:
                found[name] = dataset_id
    return {name: found[name] for name in names if name in found}


# -- PrSMs ---------------------------------------------------------------------

PRSM_SCHEMA = {
    'prsm_row': pl.Int64, 'key': pl.Utf8, 'Run': pl.Utf8, 'Scan': pl.Int64, 'RT': pl.Float64,
    'accession': pl.Utf8, 'Proforma': pl.Utf8, 'ProteoformMass': pl.Float64,
    'PrecursorMass': pl.Float64, 'qvalue': pl.Float64, 'MassErrorPPM': pl.Float64,
}


def _number(column):
    value = pl.col(column).cast(pl.Float64, strict=False)
    return pl.when(value < 0).then(None).otherwise(value)


def prsms(identifications, level):
    """The PrSMs of FLASHQuantConsensus's -out_id (contract.read_identifications),
    one row each, with ``key``: the id of the entity of ``level`` they
    identify in the comparison (the proteoform key, the protein accession or
    the consensus feature). ``prsm_row`` numbers them."""
    if identifications is None or identifications.is_empty():
        return pl.DataFrame(schema=PRSM_SCHEMA)
    key = {
        'proteoform': pl.col('proteoform'),
        'protein': pl.col('ProteinAccession'),
        'feature': pl.col('ConsensusIndex').cast(pl.Utf8),
    }[level]
    columns = identifications.columns
    frame = identifications.select(
        key.alias('key'), pl.col('Run'),
        pl.col('Scan').cast(pl.Int64, strict=False),
        (_number('RetentionTime') / 60 if 'RetentionTime' in columns else pl.lit(None, pl.Float64)).alias('RT'),
        pl.col('ProteinAccession').alias('accession'), pl.col('Proforma'),
        _number('ProteoformMass').alias('ProteoformMass'),
        (_number('PrecursorMass') if 'PrecursorMass' in columns else pl.lit(None, pl.Float64)).alias('PrecursorMass'),
        _number('ProteoformLevelQvalue').alias('qvalue'),
        (pl.col('MassErrorPPM').cast(pl.Float64, strict=False) if 'MassErrorPPM' in columns
         else pl.lit(None, pl.Float64)).alias('MassErrorPPM'),
    ).unique(['key', 'Run', 'Scan'], maintain_order=True).sort(['key', 'Run', 'Scan'])
    return frame.with_row_index('prsm_row').with_columns(pl.col('prsm_row').cast(pl.Int64)).select(
        list(PRSM_SCHEMA)).cast(PRSM_SCHEMA)


def comparison_prsms(file_manager, dataset_id, level):
    """The PrSMs a comparison was made with (empty without identifications)."""
    if not file_manager.result_exists(dataset_id, 'identifications_tsv'):
        return pl.DataFrame(schema=PRSM_SCHEMA)
    path = file_manager.get_results(dataset_id, ['identifications_tsv'])['identifications_tsv']
    return prsms(contract.read_identifications(path), level)


def prsm_in_run(hits, scan):
    """(index, scan_idx) of the best hit of spectrum ``scan`` in a run's
    protein table (``insight_data.proteins``), or None."""
    if scan is None or hits is None or 'Scan' not in hits.columns:
        return None
    row = hits.filter(pl.col('Scan') == int(scan))
    if 'is_best' in row.columns:
        best = row.filter(pl.col('is_best'))
        row = best if best.height else row
    if row.height == 0:
        return None
    first = row.row(0, named=True)
    return int(first['index']), (int(first['scan_idx']) if first.get('scan_idx') is not None else None)


# -- reporter ions --------------------------------------------------------------

def reporter_ions(frames, channels, samples):
    """Every reporter ion of every labelled spectrum, one row each:
    ``fraction``, ``Scan``, ``channel``, ``sample`` and ``intensity`` (null
    where the spectrum has no ion of the label). ``frames`` maps a fraction
    to its stored reporter table (``isobaric_dfs``); ``channels`` and
    ``samples`` are parallel (the comparison's manifest)."""
    parts = []
    names = dict(zip(channels, samples))
    for fraction, frame in frames.items():
        present = [c for c in channels if c in frame.columns]
        if not present or 'Scan' not in frame.columns:
            continue
        long = frame.select('Scan', *present).unpivot(
            index='Scan', on=present, variable_name='channel', value_name='intensity',
        )
        parts.append(long.with_columns(
            pl.lit(fraction).alias('fraction'),
            pl.col('intensity').cast(pl.Float64),
        ))
    schema = {'fraction': pl.Utf8, 'Scan': pl.Int64, 'channel': pl.Utf8, 'sample': pl.Utf8, 'intensity': pl.Float64}
    if not parts:
        return pl.DataFrame(schema=schema)
    out = pl.concat(parts, how='vertical_relaxed').with_columns(
        pl.when(pl.col('intensity') > 0).then(pl.col('intensity')).otherwise(None).alias('intensity'),
        pl.col('channel').replace_strict(names, default=None).alias('sample'),
    )
    return out.select(list(schema)).cast(schema)


def labels_found(ions, channels):
    """How many spectra carry 0, 1, ... all labels: ``labels``, ``spectra``."""
    per_spectrum = ions.group_by('fraction', 'Scan').agg(
        pl.col('intensity').is_not_null().sum().cast(pl.Int64).alias('labels'),
    )
    counts = per_spectrum.group_by('labels').agg(pl.len().cast(pl.Int64).alias('spectra'))
    grid = pl.DataFrame({'labels': list(range(len(channels) + 1))}, schema={'labels': pl.Int64})
    return grid.join(counts, on='labels', how='left').with_columns(pl.col('spectra').fill_null(0)).sort('labels')


def channel_presence(ions, channels):
    """Per label, the share (%) of the labelled spectra that carry its
    reporter ion and the median log2 intensity where it does: ``channel``,
    ``sample``, ``present``, ``median``. A spectrum counts as labelled when it
    carries any reporter ion."""
    labelled = ions.group_by('fraction', 'Scan').agg(pl.col('intensity').is_not_null().any().alias('any')).filter('any')
    spectra = max(labelled.height, 1)
    kept = ions.join(labelled.select('fraction', 'Scan'), on=['fraction', 'Scan'], how='inner')
    stats = kept.group_by('channel').agg(
        pl.col('sample').first(),
        (100.0 * pl.col('intensity').is_not_null().sum() / spectra).alias('present'),
        pl.col('intensity').log(2).median().alias('median'),
    )
    order = pl.DataFrame({'channel': list(channels)})
    return order.join(stats, on='channel', how='left').with_columns(pl.col('present').fill_null(0.0))


def reporter_ratios(ions, reference):
    """log2 ratio of every reporter ion to the ion of the ``reference``
    channel in the same spectrum, for the spectra carrying both:
    ``channel``, ``sample``, ``log2_ratio``."""
    ref = ions.filter((pl.col('channel') == reference) & pl.col('intensity').is_not_null()).select(
        'fraction', 'Scan', pl.col('intensity').alias('_ref'),
    )
    return (
        ions.filter((pl.col('channel') != reference) & pl.col('intensity').is_not_null())
        .join(ref, on=['fraction', 'Scan'], how='inner')
        .select('channel', 'sample', (pl.col('intensity') / pl.col('_ref')).log(2).alias('log2_ratio'))
    )


# -- consensus features ----------------------------------------------------------

def feature_overview(long, annotations, samples):
    """One row per consensus feature: ``RowIndex``, ``id``, ``label``,
    ``mass``, ``rt`` (min), ``MinCharge``, ``MaxCharge``, ``Detected`` (runs
    that detected it), ``Matched`` (runs it was matched into), ``Missing`` and
    ``proteoform`` (the proteoform it was assigned to, when the runs were
    identified), most often detected first."""
    runs = long.filter(pl.col('sample').is_in(list(samples)))
    counts = runs.group_by('id').agg(
        (pl.col('status') == 'detected').sum().cast(pl.Int64).alias('Detected'),
        pl.col('status').is_in(['transferred', 'transferred_identified']).sum().cast(pl.Int64).alias('Matched'),
    ).with_columns((len(samples) - pl.col('Detected') - pl.col('Matched')).alias('Missing'))
    keep = [c for c in ('label', 'mass', 'rt', 'MinCharge', 'MaxCharge', 'parent') if c in annotations.columns]
    table = annotations.select('id', *keep).join(counts, on='id', how='inner')
    if 'parent' in table.columns:
        table = table.rename({'parent': 'proteoform'})
    else:
        table = table.with_columns(pl.lit(None, pl.Utf8).alias('proteoform'))
    sort = ['Detected', 'Matched'] + (['mass'] if 'mass' in table.columns else [])
    return (
        table.sort(sort, descending=[True, True] + [False] * (len(sort) - 2))
        .with_row_index('RowIndex').with_columns(pl.col('RowIndex').cast(pl.Int64))
    )


def feature_values(long, overview, samples):
    """Every consensus feature in every run: ``RowIndex`` (of
    ``feature_overview``), ``Run``, ``log2_intensity`` and ``Found``
    (detected, matched or missing)."""
    found = (
        pl.when(pl.col('status') == 'detected').then(pl.lit('detected'))
        .when(pl.col('status').is_in(['transferred', 'transferred_identified'])).then(pl.lit('matched'))
        .otherwise(pl.lit('missing'))
    )
    return (
        long.filter(pl.col('sample').is_in(list(samples)))
        .join(overview.select('id', 'RowIndex'), on='id', how='inner')
        .select(
            'RowIndex', pl.col('sample').alias('Run'),
            pl.when(pl.col('intensity') > 0).then(pl.col('intensity').log(2)).otherwise(None).alias('log2_intensity'),
            found.alias('Found'),
        )
        .sort('RowIndex')
    )
