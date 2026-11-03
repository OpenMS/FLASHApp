"""Isobaric labelling (TMT, iTRAQ): the reporter ions FLASHDeconv measures in
every fragmentation spectrum (-iq:type, -out_quant), and the comparison of
the labelled samples in the shared quantification format
(src/quant/contract.py).

In an isobaric experiment every run holds all samples, one per channel. The
runs of one workflow start are treated as fractions of the same labelled
samples, so a sample's value sums over the runs. Levels:

- ``feature``: a precursor mass, its spectra linked across runs within
  ``MASS_TOLERANCE_PPM``;
- ``proteoform`` and ``protein`` (identified runs): the proteoform the
  spectrum was identified as (FLASHTnT PrSMs at ``QVALUE`` or below), summed.

Streamlit-free.
"""

from pathlib import Path

import polars as pl

from src.quant import contract

# FLASHDeconv's -iq:type: label, then (channel name, reporter m/z) in channel order
PLEXES = {
    'tmt6plex': ('TMT 6-plex', [
        ('126', 126.127726), ('127', 127.124761), ('128', 128.134436),
        ('129', 129.131471), ('130', 130.141145), ('131', 131.138180),
    ]),
    'tmt10plex': ('TMT 10-plex', [
        ('126', 126.127726), ('127N', 127.124761), ('127C', 127.131081), ('128N', 128.128116),
        ('128C', 128.134436), ('129N', 129.131471), ('129C', 129.137790), ('130N', 130.134825),
        ('130C', 130.141145), ('131', 131.138180),
    ]),
    'tmt11plex': ('TMT 11-plex', [
        ('126', 126.127726), ('127N', 127.124761), ('127C', 127.131081), ('128N', 128.128116),
        ('128C', 128.134436), ('129N', 129.131471), ('129C', 129.137790), ('130N', 130.134825),
        ('130C', 130.141145), ('131N', 131.138180), ('131C', 131.144500),
    ]),
    'tmt16plex': ('TMTpro 16-plex', [
        ('126', 126.127726), ('127N', 127.124761), ('127C', 127.131081), ('128N', 128.128116),
        ('128C', 128.134436), ('129N', 129.131471), ('129C', 129.137790), ('130N', 130.134825),
        ('130C', 130.141145), ('131N', 131.138180), ('131C', 131.144500), ('132N', 132.141535),
        ('132C', 132.147855), ('133N', 133.144890), ('133C', 133.151210), ('134N', 134.148245),
    ]),
    'tmt18plex': ('TMTpro 18-plex', [
        ('126', 126.127726), ('127N', 127.124761), ('127C', 127.131081), ('128N', 128.128116),
        ('128C', 128.134436), ('129N', 129.131471), ('129C', 129.137790), ('130N', 130.134825),
        ('130C', 130.141145), ('131N', 131.138180), ('131C', 131.144500), ('132N', 132.141535),
        ('132C', 132.147855), ('133N', 133.144890), ('133C', 133.151210), ('134N', 134.148245),
        ('134C', 134.154565), ('135N', 135.151600),
    ]),
    'itraq4plex': ('iTRAQ 4-plex', [
        ('114', 114.1112), ('115', 115.1082), ('116', 116.1116), ('117', 117.1149),
    ]),
    'itraq8plex': ('iTRAQ 8-plex', [
        ('113', 113.1078), ('114', 114.1112), ('115', 115.1082), ('116', 116.1116),
        ('117', 117.1149), ('118', 118.1120), ('119', 119.1153), ('121', 121.1220),
    ]),
}
DEFAULT_PLEX = 'tmt10plex'

MASS_TOLERANCE_PPM = 10.0
QVALUE = 0.01
CHANNELS_FILE = 'isobaric-channels.tsv'

_FIXED = {
    'ScanNum': ('Scan', pl.Int64), 'PrecursorScanNum': ('PrecursorScan', pl.Int64),
    'PrecursorFeatureIndex': ('Feature', pl.Int64), 'PrecursorMZ': ('PrecursorMZ', pl.Float64),
    'PrecursorCharge': ('Charge', pl.Int64), 'PrecursorMonoisotopicMass': ('PrecursorMass', pl.Float64),
    'PrecursorQscore': ('Qscore', pl.Float64), 'RetentionTime': ('RT', pl.Float64),
    'QuantifiedChannelCount': ('QuantifiedChannels', pl.Int64),
}


def label(plex):
    return PLEXES.get(plex, (plex, []))[0]


def channel_names(plex):
    return [name for name, _ in PLEXES.get(plex, ('', []))[1]]


def reporter_mz(plex):
    return [mz for _, mz in PLEXES.get(plex, ('', []))[1]]


def read_quant(path, plex):
    """FLASHDeconv -out_quant: one row per fragmentation spectrum with its
    precursor (Scan, PrecursorScan, Feature, PrecursorMZ, Charge,
    PrecursorMass, Qscore, RT in minutes, QuantifiedChannels) and one column
    per channel (the corrected reporter intensity, named as the channel).
    An empty file (no reporter ions found) gives an empty table."""
    names = channel_names(plex)
    schema = {new: dtype for new, dtype in _FIXED.values()} | {name: pl.Float64 for name in names}
    if path is None or not Path(path).exists() or Path(path).stat().st_size == 0:
        return pl.DataFrame(schema=schema)
    table = pl.read_csv(path, separator='\t', infer_schema_length=0, null_values=['', 'NaN', 'nan'])
    if table.is_empty():
        return pl.DataFrame(schema=schema)
    columns = [
        pl.col(old).cast(pl.Float64, strict=False).cast(dtype, strict=False).alias(new)
        if dtype == pl.Int64 else pl.col(old).cast(dtype, strict=False).alias(new)
        for old, (new, dtype) in _FIXED.items() if old in table.columns
    ]
    for i, name in enumerate(names, start=1):
        column = f'QuantForCh{i}'
        columns.append(
            (pl.col(column).cast(pl.Float64, strict=False) if column in table.columns
             else pl.lit(None, pl.Float64)).alias(name)
        )
    frame = table.select(columns).with_columns((pl.col('RT') / 60).alias('RT'))
    return frame.select(list(schema)).cast(schema)


def default_channels(plex):
    """channel, sample, condition: one sample per channel, no conditions yet."""
    names = channel_names(plex)
    return pl.DataFrame({'channel': names, 'sample': names, 'condition': [''] * len(names)})


def align_channels(table, plex):
    """The channel table for `plex`: rows of `table` for its channels, the
    rest from the default; sample names made unique and never empty."""
    default = default_channels(plex)
    if table is None or table.is_empty():
        return default
    table = table.select(
        pl.col('channel').cast(pl.Utf8),
        pl.col('sample').cast(pl.Utf8).fill_null(''),
        pl.col('condition').cast(pl.Utf8).fill_null(''),
    ).unique('channel', keep='first')
    merged = default.select('channel').join(table, on='channel', how='left').with_columns(
        pl.when(pl.col('sample').is_null() | (pl.col('sample').str.strip_chars() == ''))
        .then(pl.col('channel')).otherwise(pl.col('sample').str.strip_chars()).alias('sample'),
        pl.col('condition').fill_null('').str.strip_chars(),
    )
    return merged.with_columns(
        pl.when(pl.col('sample').is_duplicated())
        .then(pl.col('sample') + ' (' + pl.col('channel') + ')')
        .otherwise(pl.col('sample')).alias('sample'),
    )


def channels_path(workflow_dir):
    return Path(workflow_dir, CHANNELS_FILE)


def load_channels(workflow_dir, plex):
    path = channels_path(workflow_dir)
    table = pl.read_csv(path, separator='\t', infer_schema_length=0) if path.exists() else None
    return align_channels(table, plex)


def save_channels(workflow_dir, table, plex):
    table = align_channels(table, plex)
    table.write_csv(channels_path(workflow_dir), separator='\t')
    return table


def read_prsms(path):
    """The best target PrSM of every scan at QVALUE or below: Scan, accession,
    description, proforma, mass, qvalue and the proteoform key."""
    schema = {
        'Scan': pl.Int64, 'accession': pl.Utf8, 'description': pl.Utf8, 'proforma': pl.Utf8,
        'mass': pl.Float64, 'qvalue': pl.Float64, 'proteoform': pl.Utf8,
    }
    if path is None or not Path(path).exists():
        return pl.DataFrame(schema=schema)
    table = pl.read_csv(path, separator='\t', infer_schema_length=0, null_values=['', 'nan', 'NaN'])
    if table.is_empty() or 'Scan' not in table.columns:
        return pl.DataFrame(schema=schema)
    for column in ('Proforma', 'ProteinDescription', 'Score', 'ProteoformLevelQvalue'):
        if column not in table.columns:
            table = table.with_columns(pl.lit(None, pl.Utf8).alias(column))
    table = table.select(
        pl.col('Scan').cast(pl.Float64, strict=False).cast(pl.Int64, strict=False),
        pl.col('ProteinAccession').alias('accession'),
        pl.col('ProteinDescription').alias('description'),
        pl.col('Proforma').alias('proforma'),
        # FLASHTnT reports -1 for a proteoform without a precursor mass (native mode)
        pl.when(pl.col('ProteoformMass').cast(pl.Float64, strict=False) > 0)
        .then(pl.col('ProteoformMass').cast(pl.Float64, strict=False)).alias('mass'),
        pl.col('ProteoformLevelQvalue').cast(pl.Float64, strict=False).alias('qvalue'),
        pl.col('Score').cast(pl.Float64, strict=False).alias('_score'),
    ).filter(
        pl.col('Scan').is_not_null() & pl.col('accession').is_not_null()
        & ~pl.col('accession').str.contains('DECOY_', literal=True)
        & (pl.col('qvalue').is_null() | (pl.col('qvalue') <= QVALUE))
    )
    table = table.sort(['Scan', '_score'], descending=[False, True], nulls_last=True).unique('Scan', keep='first')
    return table.with_columns(
        pl.when(pl.col('proforma').is_not_null())
        .then(pl.col('accession') + '|' + pl.col('proforma'))
        .otherwise(pl.col('accession') + '|' + pl.col('mass').round(1).cast(pl.Utf8) + ' Da')
        .alias('proteoform'),
    ).select(list(schema)).sort('Scan')


def _cluster_masses(masses, ppm=MASS_TOLERANCE_PPM):
    """Cluster index of every mass: sorted masses within `ppm` of the
    previous one share a cluster."""
    order = sorted(range(len(masses)), key=lambda i: masses[i])
    clusters = [0] * len(masses)
    current, previous = -1, None
    for i in order:
        if previous is None or masses[i] - previous > previous * ppm * 1e-6:
            current += 1
        clusters[i] = current
        previous = masses[i]
    return clusters


def _long(scans, key, channels):
    """Sum the channel intensities of `scans` per `key`: the contract's long
    table (id, sample, intensity, status, score, qvalue)."""
    parts = []
    for row in channels.iter_rows(named=True):
        parts.append(scans.group_by(key).agg(
            pl.col(row['channel']).fill_nan(None).sum().alias('intensity'),
        ).select(
            pl.col(key).cast(pl.Utf8).alias('id'), pl.lit(row['sample']).alias('sample'), 'intensity',
        ))
    long = pl.concat(parts)
    return long.with_columns(
        pl.when(pl.col('intensity') > 0).then(pl.col('intensity')).otherwise(None).alias('intensity'),
    ).with_columns(
        pl.when(pl.col('intensity').is_not_null()).then(pl.lit('detected')).otherwise(pl.lit('missing')).alias('status'),
        pl.lit(None, pl.Float64).alias('score'),
        pl.lit(None, pl.Float64).alias('qvalue'),
    ).sort('id', 'sample')


def build(runs, channels):
    """The levels of an isobaric comparison, {level: (long, annotations)}.

    `runs` lists (run name, read_quant frame, read_prsms frame or None);
    `channels` is the channel table (channel, sample, condition). Spectra
    without a precursor mass (and no identification) or without any
    reporter ion are left out."""
    names = list(channels['channel'])
    frames = []
    for run, quant, prsms in runs:
        frame = quant.with_columns(pl.lit(run).alias('Run'))
        if prsms is not None and not prsms.is_empty():
            frame = frame.join(prsms, on='Scan', how='left')
        else:
            frame = frame.with_columns(
                pl.lit(None, pl.Utf8).alias('proteoform'), pl.lit(None, pl.Utf8).alias('accession'),
                pl.lit(None, pl.Utf8).alias('description'), pl.lit(None, pl.Utf8).alias('proforma'),
                pl.lit(None, pl.Float64).alias('mass'), pl.lit(None, pl.Float64).alias('qvalue'),
            )
        frames.append(frame)
    if not frames:
        return {}
    # Without an MS1 survey (native runs) FLASHDeconv has no precursor mass:
    # the identified proteoform's mass stands in, and where FLASHTnT reports
    # none either (-1 in native mode), the proteoform itself is the feature
    scans = pl.concat(frames, how='diagonal_relaxed').with_columns(
        pl.when(pl.col('PrecursorMass').fill_nan(None) > 0).then(pl.col('PrecursorMass'))
        .when(pl.col('mass').fill_nan(None) > 0).then(pl.col('mass'))
        .otherwise(None).alias('PrecursorMass'),
    ).filter(
        (pl.col('PrecursorMass').is_not_null() | pl.col('proteoform').is_not_null())
        & pl.any_horizontal([pl.col(n).fill_nan(None).fill_null(0) > 0 for n in names])
    )
    if scans.is_empty():
        return {}
    with_mass = scans.filter(pl.col('PrecursorMass').is_not_null())
    with_mass = with_mass.with_columns(
        pl.Series('feature', _cluster_masses(with_mass['PrecursorMass'].to_list()), dtype=pl.Int64),
    )
    offset = (with_mass['feature'].max() + 1) if not with_mass.is_empty() else 0
    by_form = scans.filter(pl.col('PrecursorMass').is_null()).with_columns(
        (pl.col('proteoform').rank('dense').cast(pl.Int64) - 1 + offset).alias('feature'),
    )
    scans = pl.concat([with_mass, by_form], how='diagonal_relaxed')

    feature_info = scans.group_by('feature').agg(
        pl.col('PrecursorMass').median().alias('mass'),
        pl.col('RT').median().alias('rt'),
        pl.len().alias('Spectra'),
        pl.col('Run').n_unique().alias('Runs'),
        pl.col('proteoform').drop_nulls().unique().alias('_forms'),
    ).with_columns(
        pl.col('feature').cast(pl.Utf8).alias('id'),
        pl.coalesce(
            pl.col('mass').round(2).cast(pl.Utf8) + ' Da @ ' + pl.col('rt').round(2).cast(pl.Utf8) + ' min',
            pl.col('_forms').list.first().str.split('|').list.last() + ' @ ' + pl.col('rt').round(2).cast(pl.Utf8) + ' min',
        ).alias('label'),
        pl.when(pl.col('_forms').list.len() == 1).then(pl.col('_forms').list.first()).otherwise(None).alias('parent'),
        (pl.col('_forms').list.len() > 1).alias('ambiguous'),
    ).drop('_forms', 'feature').sort('mass')
    feature_long = _long(scans, 'feature', channels).join(
        feature_info.select('id', 'parent'), on='id', how='left',
    )
    levels = {'feature': (feature_long.select(list(contract.LONG_SCHEMA)).cast(contract.LONG_SCHEMA), feature_info)}

    identified = scans.filter(pl.col('proteoform').is_not_null())
    if identified.is_empty():
        return levels
    forms = identified.group_by('proteoform').agg(
        pl.col('accession').first(),
        pl.col('description').drop_nulls().first(),
        pl.col('proforma').drop_nulls().first(),
        pl.col('mass').filter(pl.col('mass') > 0).median(),
        pl.col('RT').median().alias('rt'),
        pl.col('qvalue').min(),
        pl.len().alias('PrSMs'),
        pl.col('Run').n_unique().alias('IdentifiedRuns'),
    ).with_columns(
        pl.col('proteoform').alias('id'),
        # Short enough for a plot label: the protein's name and the start of the sequence
        (pl.col('accession').str.split('|').list.last().str.replace('_.*$', '') + ' ' + pl.coalesce(
            pl.when(pl.col('proforma').str.len_chars() > 12)
            .then(pl.col('proforma').str.slice(0, 12) + '…').otherwise(pl.col('proforma')),
            pl.col('mass').round(1).cast(pl.Utf8) + ' Da',
        )).alias('label'),
    )
    form_long = _long(identified, 'proteoform', channels).join(
        forms.select('id', pl.col('accession').alias('parent')), on='id', how='left',
    )
    levels['proteoform'] = (form_long.select(list(contract.LONG_SCHEMA)).cast(contract.LONG_SCHEMA), forms)

    protein_long = _long(identified, 'accession', channels).with_columns(pl.lit(None, pl.Utf8).alias('parent'))
    proteins = forms.group_by('accession').agg(
        pl.col('description').drop_nulls().first(),
        pl.len().alias('proteoforms'),
        pl.col('PrSMs').sum(),
    ).with_columns(
        pl.col('accession').alias('id'),
        pl.col('accession').str.split('|').list.get(-1, null_on_oob=True).alias('label'),
    )
    levels['protein'] = (protein_long.select(list(contract.LONG_SCHEMA)).cast(contract.LONG_SCHEMA), proteins)
    return levels


def design(channels):
    """sample -> (condition, replicate) for session.store: replicates
    numbered per condition in channel order."""
    counts, result = {}, {}
    for row in channels.iter_rows(named=True):
        counts[row['condition']] = counts.get(row['condition'], 0) + 1
        result[row['sample']] = (row['condition'], counts[row['condition']])
    return result
