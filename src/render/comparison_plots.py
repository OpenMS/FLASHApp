"""Plots of the identification and quality-control views of a comparison.

Each function takes a frame from src/quant/qc.py and returns the settings of an
openms-insight ``SeriesPlot`` (bars or boxes) with its tidy table as ``data``;
the cells (src/render/quant_cells.py) draw it with ``series.show_series``.
"""

import polars as pl

from src.render.series import thin

# Detected by the run, matched with and without an identification of the run
HOW_COLORS = {'detected': '#3157e9', 'matched, identified': '#7f9cf5', 'matched': '#ff7501'}
OUTCOME_COLORS = {'matched, identified': '#7f9cf5', 'matched': '#ff7501', 'not matched': '#9E9E9E'}


def counts_per_run(counts, noun, sample_noun='Run'):
    """Bars: entities each run detected and matched, side by side."""
    shown = [how for how in HOW_COLORS
             if how == 'detected' or counts.filter(counts['status'] == how)['count'].sum()]
    data = (
        counts.filter(pl.col('status').is_in(shown)).select('sample', 'status', 'count')
        .with_columns(pl.format('{}: {} {}', 'status', 'count', pl.lit(noun)).alias('hover'))
    )
    return dict(
        data=data, x_column='sample', y_column='count', series_column='status', series_colors=HOW_COLORS,
        mode='bars', hover_column='hover', title=f'{noun.capitalize()} per {sample_noun.lower()}',
        x_label=sample_noun, y_label=noun.capitalize(),
    )


def frequency(table, noun, sample_noun='Run'):
    """Grouped bars: entities found in 1 ... n runs."""
    data = (
        table.select(pl.col('runs').cast(pl.Utf8), 'detected', 'found')
        .unpivot(index='runs', variable_name='how', value_name='count')
    )
    return dict(
        data=data, x_column='runs', y_column='count', series_column='how',
        series_names={'detected': 'detected', 'found': 'detected or matched'},
        series_colors={'detected': HOW_COLORS['detected'], 'found': HOW_COLORS['matched']},
        mode='bars', title=f'{noun.capitalize()} found in how many {sample_noun.lower()}s',
        x_label=f'{sample_noun}s', y_label=noun.capitalize(),
    )


def _median(values):
    values = sorted(values)
    n = len(values)
    if not n:
        return None
    return values[n // 2] if n % 2 else (values[n // 2 - 1] + values[n // 2]) / 2


def totals(frame, colors, conditions, sample_noun='Run'):
    """Bars: summed intensity per run (sample), coloured by condition."""
    median = _median([t for t in frame['total'].to_list() if t > 0])
    data = frame.with_columns(
        pl.col('sample').replace_strict(conditions, default='no condition').alias('condition'),
        pl.format(
            '{}: total {}, {} quantified, {}% of the median', 'sample',
            pl.col('total').round_sig_figs(3), 'quantified', (pl.col('relative').fill_null(0) * 100).round(0),
        ).alias('hover'),
    ).select('sample', 'total', 'condition', 'hover')
    title = f'Total intensity per {sample_noun.lower()}'
    if median is not None:
        title += f' (median {median:.3g})'
    return dict(
        data=data, x_column='sample', y_column='total', series_column='condition',
        series_colors={**colors, 'no condition': '#9E9E9E'}, mode='bars', hover_column='hover',
        title=title, x_label=sample_noun, y_label='Summed intensity',
    )


def coefficients_of_variation(frame, colors):
    """Box plots: CV of the entities within each condition."""
    medians = [
        f"{condition} {_median(frame.filter(frame['condition'] == condition)['cv'].to_list()):.0f}%"
        for condition in sorted(set(frame['condition'].to_list()))
    ]
    return dict(
        data=thin(frame.select('condition', 'cv').sort('condition'), ['condition'], 'cv'),
        x_column='condition', y_column='cv', series_column='condition', series_colors=colors, mode='box',
        show_legend=False, title=f'Variation between replicates (CV), median {", ".join(medians)}',
        x_label='Condition', y_label='CV (%)', y_range_mode='nonnegative',
    )


def match_quality(frame, samples, sample_noun='Run'):
    """Box plots: match score of the values matched between runs and of the
    rejected matches, per run."""
    order = {s: i for i, s in enumerate(samples)}
    data = frame.select('sample', 'score', 'outcome').sort(
        pl.col('sample').replace_strict(order, default=len(order)))
    return dict(
        data=thin(data, ['sample', 'outcome'], 'score'), x_column='sample', y_column='score',
        series_column='outcome', series_colors=OUTCOME_COLORS, mode='box', title='Match between runs: match score',
        x_label=sample_noun, y_label='Match score', y_range_mode='normal',
    )


# -- reporter ions (isobaric labels) -------------------------------------------

def reporter_intensities(presence, ions, colors, conditions):
    """Box plots: log2 intensity of every label's reporter ions over the
    labelled spectra, coloured by condition; the share of spectra carrying the
    ion is part of the label."""
    parts = []
    for row in presence.iter_rows(named=True):
        channel, sample = row['channel'], row['sample'] or row['channel']
        values = ions.filter(ions['channel'] == channel)['intensity'].drop_nulls()
        values = values.filter(values > 0)
        parts.append(pl.DataFrame({
            'label': f'{channel} ({sample}) · {row["present"]:.0f}%',
            'log2_intensity': values.log(2),
            'condition': conditions.get(sample, 'no condition'),
        }, schema={'label': pl.Utf8, 'log2_intensity': pl.Float64, 'condition': pl.Utf8}))
    data = pl.concat(parts) if parts else pl.DataFrame(
        schema={'label': pl.Utf8, 'log2_intensity': pl.Float64, 'condition': pl.Utf8})
    return dict(
        data=thin(data, ['label'], 'log2_intensity'), x_column='label', y_column='log2_intensity',
        series_column='condition', series_colors={**colors, 'no condition': '#9E9E9E'}, mode='box',
        title='Reporter ion intensity per label (share of spectra with the ion)',
        x_label='Label (sample)', y_label='log2 intensity', y_range_mode='normal',
    )


def reporter_ratios(ratios, reference, samples):
    """Box plots: log2 ratio of every label to the reference label in the same
    spectrum (0 is equal intensity)."""
    parts = [
        pl.DataFrame({'label': f'{channel} ({sample})', 'log2_ratio': values},
                     schema={'label': pl.Utf8, 'log2_ratio': pl.Float64})
        for channel, sample in samples
        if len(values := ratios.filter(ratios['channel'] == channel)['log2_ratio'].drop_nulls())
    ]
    data = pl.concat(parts) if parts else pl.DataFrame(schema={'label': pl.Utf8, 'log2_ratio': pl.Float64})
    return dict(
        data=thin(data, ['label'], 'log2_ratio'), x_column='label', y_column='log2_ratio', mode='box',
        series_colors=None, show_legend=False, title=f'Ratio to label {reference} in the same spectrum',
        x_label='Label (sample)', y_label=f'log2 (label / {reference})', y_range_mode='normal',
    )


def labels_found(table, plex):
    """Bars: spectra carrying 0, 1, ... all reporter ions."""
    total = int(table['spectra'].sum()) or 1
    data = table.select(
        pl.col('labels').cast(pl.Utf8), 'spectra',
        pl.format('{} labels: {} spectra ({}%)', 'labels', 'spectra',
                  (pl.col('spectra') * 100.0 / total).round(0)).alias('hover'),
    )
    return dict(
        data=data, x_column='labels', y_column='spectra', mode='bars', hover_column='hover', show_legend=False,
        title=f'Reporter ions found per spectrum ({plex})', x_label='Labels with a reporter ion',
        y_label='Spectra',
    )


# -- consensus features (label-free) -------------------------------------------

FOUND_COLORS = {'detected': '#3157e9', 'matched': '#ff7501', 'missing': '#9E9E9E'}


def feature_intensity(values, samples, sample_noun='Run'):
    """Box plots per run: log2 intensity of the features the run detected and
    of those matched into it from other runs."""
    order = {s: i for i, s in enumerate(samples)}
    data = (
        values.filter(pl.col('Found').is_in(['detected', 'matched']))
        .select('Run', 'log2_intensity', 'Found')
        .sort(pl.col('Run').replace_strict(order, default=len(order)))
    )
    return dict(
        data=thin(data, ['Run', 'Found'], 'log2_intensity'), x_column='Run', y_column='log2_intensity',
        series_column='Found', series_colors=FOUND_COLORS, mode='box',
        title=f'Feature intensity per {sample_noun.lower()}', x_label=sample_noun, y_label='log2 intensity',
        y_range_mode='normal',
    )
