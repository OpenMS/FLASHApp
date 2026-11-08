"""Plots of the identification and quality-control views of a comparison.

Each function takes a frame from src/quant/qc.py and returns a Plotly figure;
the cells (src/render/quant_cells.py) show it with ``show_fig``.
"""

import plotly.graph_objects as go

# Detected by the run, matched with and without an identification of the run
HOW_COLORS = {'detected': '#3157e9', 'matched, identified': '#7f9cf5', 'matched': '#ff7501'}
OUTCOME_COLORS = {'matched, identified': '#7f9cf5', 'matched': '#ff7501', 'not matched': '#9E9E9E'}


def counts_per_run(counts, noun, sample_noun='Run'):
    """Stacked bars: entities each run detected and matched."""
    fig = go.Figure()
    for how, color in HOW_COLORS.items():
        part = counts.filter(counts['status'] == how)
        if part['count'].sum() == 0 and how != 'detected':
            continue
        fig.add_bar(x=part['sample'].to_list(), y=part['count'].to_list(), name=how, marker_color=color,
                    hovertemplate=f'%{{x}}<br>{how}: %{{y}}<extra></extra>')
    fig.update_layout(
        barmode='stack', showlegend=True, title=f'{noun.capitalize()} per {sample_noun.lower()}',
        xaxis_title=sample_noun, yaxis_title=noun.capitalize(), legend_title_text='Found',
    )
    return fig


def frequency(table, noun, sample_noun='Run'):
    """Grouped bars: entities found in 1 ... n runs."""
    fig = go.Figure()
    x = [str(k) for k in table['runs'].to_list()]
    fig.add_bar(x=x, y=table['detected'].to_list(), name='detected', marker_color=HOW_COLORS['detected'],
                hovertemplate=f'Detected in %{{x}} {sample_noun.lower()}s: %{{y}}<extra></extra>')
    fig.add_bar(x=x, y=table['found'].to_list(), name='detected or matched', marker_color=HOW_COLORS['matched'],
                hovertemplate=f'Found in %{{x}} {sample_noun.lower()}s: %{{y}}<extra></extra>')
    fig.update_layout(
        barmode='group', title=f'{noun.capitalize()} found in how many {sample_noun.lower()}s',
        xaxis_title=f'{sample_noun}s', yaxis_title=noun.capitalize(), legend_title_text='',
    )
    return fig


def totals(frame, colors, conditions, sample_noun='Run'):
    """Bars: summed intensity per run (sample), coloured by condition, with
    the median as a line."""
    fig = go.Figure()
    samples = frame['sample'].to_list()
    groups = {}
    for sample in samples:
        groups.setdefault(conditions.get(sample, 'no condition'), []).append(sample)
    by_sample = {r['sample']: r for r in frame.iter_rows(named=True)}
    for condition, members in groups.items():
        fig.add_bar(
            x=members, y=[by_sample[s]['total'] for s in members], name=condition,
            marker_color=colors.get(condition, '#9E9E9E'),
            customdata=[[by_sample[s]['quantified'], (by_sample[s]['relative'] or 0) * 100] for s in members],
            hovertemplate='%{x}<br>Total %{y:.3s}<br>%{customdata[0]} quantified<br>'
                          '%{customdata[1]:.0f}% of the median<extra></extra>',
        )
    positive = sorted(t for t in frame['total'].to_list() if t > 0)
    if positive:
        n = len(positive)
        median = positive[n // 2] if n % 2 else (positive[n // 2 - 1] + positive[n // 2]) / 2
        fig.add_hline(y=median, line_dash='dot', line_color='grey', annotation_text='median')
    fig.update_layout(
        title=f'Total intensity per {sample_noun.lower()}', xaxis_title=sample_noun,
        yaxis_title='Summed intensity', xaxis={'categoryorder': 'array', 'categoryarray': samples},
        legend_title_text='Condition',
    )
    return fig


def coefficients_of_variation(frame, colors):
    """Box plots: CV of the entities within each condition."""
    fig = go.Figure()
    for condition in sorted(set(frame['condition'].to_list())):
        values = frame.filter(frame['condition'] == condition)['cv'].to_list()
        fig.add_box(y=values, name=condition, marker_color=colors.get(condition, '#9E9E9E'), boxmean=True,
                    hovertemplate=f'{condition}<br>CV %{{y:.1f}}%<extra></extra>')
        if values:
            median = sorted(values)[len(values) // 2]
            fig.add_annotation(x=condition, y=median, text=f'median {median:.0f}%', showarrow=False,
                               yshift=12, font={'size': 11})
    fig.update_layout(title='Variation between replicates (CV)', yaxis_title='CV (%)', xaxis_title='Condition',
                      showlegend=False, yaxis_rangemode='tozero')
    return fig


def match_quality(frame, samples, sample_noun='Run'):
    """Box plots: match score of the values matched between runs and of the
    rejected matches, per run."""
    fig = go.Figure()
    for outcome, color in OUTCOME_COLORS.items():
        part = frame.filter(frame['outcome'] == outcome)
        if part.height == 0:
            continue
        fig.add_box(x=part['sample'].to_list(), y=part['score'].to_list(), name=outcome, marker_color=color,
                    boxpoints='outliers')
    fig.update_layout(
        boxmode='group', title='Match between runs: match score',
        xaxis_title=sample_noun, yaxis_title='Match score', legend_title_text='',
        xaxis={'categoryorder': 'array', 'categoryarray': list(samples)},
    )
    return fig


# -- reporter ions (isobaric labels) -------------------------------------------

def reporter_intensities(presence, ions, colors, conditions):
    """Box plots: log2 intensity of every label's reporter ions over the
    labelled spectra, coloured by condition, with the share of spectra
    carrying the ion above each box."""
    fig = go.Figure()
    for row in presence.iter_rows(named=True):
        channel, sample = row['channel'], row['sample'] or row['channel']
        values = ions.filter(ions['channel'] == channel)['intensity'].drop_nulls()
        condition = conditions.get(sample, 'no condition')
        name = f'{channel}<br>{sample}'
        fig.add_box(
            y=values.log(2).to_list(), name=name, marker_color=colors.get(condition, '#9E9E9E'), boxpoints=False,
            hovertemplate=f'{channel} ({sample}, {condition})<br>log2 %{{y:.2f}}<extra></extra>',
        )
        fig.add_annotation(x=name, y=1, yref='paper', yanchor='bottom', text=f'{row["present"]:.0f}%',
                           showarrow=False, font={'size': 11})
    fig.update_layout(
        title='Reporter ion intensity per label (share of spectra with the ion on top)',
        xaxis_title='Label (sample)', yaxis_title='log2 intensity', showlegend=False, margin={'t': 70},
    )
    return fig


def reporter_ratios(ratios, reference, samples):
    """Box plots: log2 ratio of every label to the reference label in the same
    spectrum, with the line of equal intensity."""
    fig = go.Figure()
    for channel, sample in samples:
        values = ratios.filter(ratios['channel'] == channel)['log2_ratio'].to_list()
        if not values:
            continue
        fig.add_box(y=values, name=f'{channel}<br>{sample}', marker_color='#3157e9', boxpoints=False,
                    hovertemplate=f'{channel} ({sample})<br>log2 ratio %{{y:.2f}}<extra></extra>')
    fig.add_hline(y=0, line_dash='dot', line_color='grey')
    fig.update_layout(
        title=f'Ratio to label {reference} in the same spectrum', xaxis_title='Label (sample)',
        yaxis_title=f'log2 (label / {reference})', showlegend=False,
    )
    return fig


def labels_found(table, plex):
    """Bars: spectra carrying 0, 1, ... all reporter ions."""
    fig = go.Figure()
    total = int(table['spectra'].sum()) or 1
    fig.add_bar(
        x=[str(n) for n in table['labels'].to_list()], y=table['spectra'].to_list(), marker_color='#3157e9',
        customdata=[100.0 * s / total for s in table['spectra'].to_list()],
        hovertemplate='%{x} labels: %{y} spectra (%{customdata:.0f}%)<extra></extra>',
    )
    fig.update_layout(
        title=f'Reporter ions found per spectrum ({plex})', xaxis_title='Labels with a reporter ion',
        yaxis_title='Spectra',
    )
    return fig


# -- consensus features (label-free) -------------------------------------------

FOUND_COLORS = {'detected': '#3157e9', 'matched': '#ff7501', 'missing': '#9E9E9E'}


def feature_intensity(values, samples, sample_noun='Run'):
    """Box plots per run: log2 intensity of the features the run detected and
    of those matched into it from other runs."""
    fig = go.Figure()
    for found in ('detected', 'matched'):
        part = values.filter((values['Found'] == found) & values['log2_intensity'].is_not_null())
        if part.height == 0:
            continue
        fig.add_box(x=part['Run'].to_list(), y=part['log2_intensity'].to_list(), name=found,
                    marker_color=FOUND_COLORS[found], boxpoints=False)
    fig.update_layout(
        boxmode='group', title=f'Feature intensity per {sample_noun.lower()}', xaxis_title=sample_noun,
        yaxis_title='log2 intensity', legend_title_text='Found',
        xaxis={'categoryorder': 'array', 'categoryarray': list(samples)},
    )
    return fig
