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
