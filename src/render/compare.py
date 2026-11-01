"""Views of a comparison of FLASHQuant runs (see src/parse/quant_compare.py).

The views and their wording follow OpenDDA (OpenMS/quantms-web): Compare
groups, Volcano plot, Sample similarity (PCA), Clustered heatmap and the
diagnostic plots of its Fine-tune pages, on proteoforms instead of proteins.

The comparison is stored as one dataset of its own cache: the linked
proteoforms (`consensus`), their feature groups per run (`members`), the
test results (`statistics`), normalized abundances per run (`normalized`,
long), the tables after each step (`filtered`, `imputed`, `normalized_wide`),
the z-scores of the most significant proteoforms (`zscores`), the
experimental design (`design`) and the settings (`settings`). The views
share the identifier `consensus` (a linked proteoform's ConsensusIndex):
clicking a point in the volcano plot or a row in the table shows its
abundances and its feature groups in every run.
"""

import polars as pl
import streamlit as st

from openms_insight import ClusteredHeatmap, PCAPlot, SeriesPlot, Table, VolcanoPlot, defer_reruns

from src.common.common import show_fig
from src.parse.differential import describe_settings
from src.render import insight_data
from src.render.postprocessing_plots import (
    imputed_value_preview, missingness_vs_intensity, pvalue_histogram, sample_distributions,
)
from src.render.render import RECLICK, TABLE_HEIGHT, MissingResult, _component, _fixed

DATASET = 'comparison'
RESULTS = [
    'consensus', 'members', 'statistics', 'normalized', 'filtered', 'imputed', 'normalized_wide',
    'zscores', 'design', 'settings',
]

CONDITION_COLORS = ['#3366CC', '#DC3912', '#FF9900', '#109618', '#990099', '#0099C6']

# OpenDDA's cut-offs for calling a change
FC_CUTOFF = 1.0
P_CUTOFF = 0.05

_P = {'formatter': 'scientific', 'formatterParams': {'precision': 3}}
_FC = {'formatter': 'signed', 'formatterParams': {'precision': 2}}

CONSENSUS_COLUMNS = [
    {'title': 'Proteoform', 'field': 'Proteoform', 'sorter': 'string',
     'headerTooltip': 'Monoisotopic mass and apex retention time of the linked feature groups.'},
    {'title': 'Change', 'field': 'Change', 'sorter': 'string',
     'headerTooltip': f'Adjusted p-value ≤ {P_CUTOFF} and |log2 fold change| ≥ {FC_CUTOFF:g}.'},
    {'title': 'log2 fold change', 'field': 'log2FC', 'sorter': 'number', **_FC},
    {'title': 'Adjusted p-value', 'field': 'p-adj', 'sorter': 'number', **_P,
     'headerTooltip': 'Empty: removed by the filter, so not tested.'},
    {'title': 'p-value', 'field': 'p-value', 'sorter': 'number', **_P},
    {'title': 'Monoisotopic Mass', 'field': 'MonoisotopicMass', 'sorter': 'number', **_fixed()},
    {'title': 'Retention Time (min)', 'field': 'RTMinutes', 'sorter': 'number', **_fixed(2)},
    {'title': 'Runs', 'field': 'Runs', 'sorter': 'number',
     'headerTooltip': 'In how many runs the proteoform was quantified.'},
]

MEMBER_COLUMNS = [
    {'title': 'Run', 'field': 'Run', 'sorter': 'string'},
    {'title': 'Feature Group', 'field': 'FeatureGroupIndex', 'sorter': 'number',
     'headerTooltip': 'Index of the feature group in the run (FLASHQuant Viewer).'},
    {'title': 'Monoisotopic Mass', 'field': 'MonoisotopicMass', 'sorter': 'number', **_fixed()},
    {'title': 'Retention Time (min)', 'field': 'RTMinutes', 'sorter': 'number', **_fixed(2)},
    {'title': 'Quantity', 'field': 'Quantity', 'sorter': 'number', **_fixed(0)},
]


def _scan(file_manager, name):
    return insight_data.scan(file_manager, DATASET, name)


def design(file_manager):
    """Run -> condition of the stored comparison."""
    frame = _scan(file_manager, 'design').collect()
    return dict(zip(frame['Run'], frame['Condition']))


def settings(file_manager):
    return file_manager.get_results(DATASET, ['settings'])['settings']


def comparison_groups(conditions):
    """(reference, compared): a positive log2FC means higher in the compared
    condition (openms_insight subtracts the first condition's mean from the
    last one's, in sorted order)."""
    names = sorted(set(conditions))
    return (names[0], names[-1]) if len(names) >= 2 else None


def condition_colors(names):
    return {name: CONDITION_COLORS[i % len(CONDITION_COLORS)] for i, name in enumerate(sorted(names))}


def change_labels(compared, fc=FC_CUTOFF, padj=P_CUTOFF):
    """OpenDDA's Change column: higher/lower in the compared condition, or –."""
    significant = (pl.col('p-adj') <= padj) & (pl.col('log2FC').abs() >= fc)
    return (
        pl.when(significant & (pl.col('log2FC') > 0)).then(pl.lit(f'higher in {compared}'))
        .when(significant & (pl.col('log2FC') < 0)).then(pl.lit(f'lower in {compared}'))
        .otherwise(pl.lit('–'))
    )


def consensus_table(file_manager, compared):
    """Linked proteoforms with their test results (null where not tested),
    most significant first."""
    consensus = _scan(file_manager, 'consensus').select(
        'ConsensusIndex', 'Proteoform', 'MonoisotopicMass', 'RT', 'Runs',
    )
    tests = _scan(file_manager, 'statistics').select('ConsensusIndex', 'log2FC', 'p-value', 'p-adj')
    return (
        consensus.join(tests, on='ConsensusIndex', how='left')
        .with_columns(
            (pl.col('RT') / 60).alias('RTMinutes'),
            change_labels(compared).alias('Change'),
        )
        .drop('RT')
        .sort('p-adj', nulls_last=True)
    )


def members(file_manager):
    return (
        _scan(file_manager, 'members')
        .with_columns((pl.col('RT') / 60).alias('RTMinutes'))
        .drop('RT')
        .with_row_index('MemberIndex')
        .with_columns(pl.col('MemberIndex').cast(pl.Int64))
    )


def count_significant(statistics, fc=FC_CUTOFF, padj=P_CUTOFF):
    """Proteoforms up and down at |log2FC| >= fc and adjusted p <= padj."""
    tested = statistics.drop_nulls(['log2FC', 'p-adj'])
    significant = tested.filter(pl.col('p-adj') <= padj)
    return (
        significant.filter(pl.col('log2FC') >= fc).height,
        significant.filter(pl.col('log2FC') <= -fc).height,
    )


def _metadata(run_conditions):
    return pl.DataFrame(
        {'sample_id': list(run_conditions), 'group': list(run_conditions.values())},
        schema={'sample_id': pl.String, 'group': pl.String},
    )


def render_comparison(file_manager, state_manager, key='flashquant_compare'):
    for name in RESULTS:
        if insight_data.result_path(file_manager, DATASET, name) is None:
            raise MissingResult(name)
    run_conditions = design(file_manager)
    runs = list(run_conditions)
    direction = comparison_groups(run_conditions.values())
    colors = condition_colors(set(run_conditions.values()))
    statistics = _scan(file_manager, 'statistics').collect()
    used = settings(file_manager)

    # One view at a time: Insight components drawn inside a hidden st.tabs panel get no height
    views = ['Result', 'Sample similarity (PCA)', 'Clustered heatmap', 'Diagnostics']
    view = st.segmented_control(
        'View', views, default=views[0], key=f'{key}_view', label_visibility='collapsed',
    ) or views[0]
    if view == views[0]:
        _render_result(file_manager, state_manager, key, statistics, direction, colors, used)
    elif view == views[1]:
        _render_pca(file_manager, state_manager, key, run_conditions, colors)
    elif view == views[2]:
        _component(
            ClusteredHeatmap, 'compare_heatmap', file_manager, DATASET, ['zscores', 'design'],
            lambda: dict(
                data=_scan(file_manager, 'zscores'), id_col='Proteoform',
                metadata=_metadata(run_conditions),
                title='Most significant proteoforms (z-scores, clustered)',
                x_label='Runs', y_label='Proteoforms',
                colorscale=[[0, '#6699E0'], [0.5, '#FFFFFF'], [1, '#E06666']],
                reversescale=False, intensity_label='Z-score', group_colors=colors,
            ),
        )(key=f'{key}_heatmap', state_manager=state_manager, height=700)
        st.caption(
            'The 50 proteoforms with the smallest adjusted p-values, each scaled to mean 0 and '
            'standard deviation 1 across the runs. Runs of one condition should cluster together.'
        )
    else:
        _render_diagnostics(file_manager, runs, run_conditions, statistics, used)


def _render_result(file_manager, state_manager, key, statistics, direction, colors, used):
    reference, compared = direction
    n_conditions = len(colors)
    up, down = count_significant(statistics)
    metrics = st.columns(3)
    metrics[0].metric('Proteoforms compared', f"{statistics['p-adj'].is_not_null().sum():,}")
    metrics[1].metric(
        f'Higher in {compared}', f'{up:,}',
        help=f'Adjusted p-value ≤ {P_CUTOFF} and at least twice as abundant.',
    )
    metrics[2].metric(
        f'Lower in {compared}', f'{down:,}',
        help=f'Adjusted p-value ≤ {P_CUTOFF} and at most half as abundant.',
    )
    note = f'A positive log2 fold change means higher in **{compared}** than in **{reference}**.'
    if n_conditions > 2:
        note += (
            ' With more than two conditions the p-value tests whether any condition differs; the '
            f'fold change still compares **{compared}** with **{reference}**.'
        )
    st.caption(
        note + f' Counts use |log2 fold change| ≥ {FC_CUTOFF:g} and adjusted p-value ≤ {P_CUTOFF}; '
        'explore other cut-offs on the volcano plot.'
    )
    with st.expander('Settings used'):
        st.markdown('\n'.join(f'- {line}' for line in describe_settings(used)))

    table = _component(
        Table, 'compare_table', file_manager, DATASET, ['consensus', 'statistics', 'design'],
        lambda: dict(
            data=consensus_table(file_manager, compared), column_definitions=CONSENSUS_COLUMNS,
            title='All linked proteoforms, most significant first', index_field='ConsensusIndex',
            interactivity={'consensus': 'ConsensusIndex'}, go_to_fields=['Proteoform'],
            reclick=RECLICK,
        ),
    )
    profile = _component(
        SeriesPlot, 'compare_profile', file_manager, DATASET, ['normalized', 'settings'],
        lambda: dict(
            data=_scan(file_manager, 'normalized'),
            filters={'consensus': 'ConsensusIndex'},
            x_column='Run', y_column='Abundance', series_column='Condition',
            series_colors=colors, mode='markers', y_range_mode='normal',
            title='Abundance per run', x_label='Run',
            y_label=f"Abundance ({used['postproc-transform']}, normalized)",
        ),
    )
    member_table = _component(
        Table, 'compare_members', file_manager, DATASET, ['members'],
        lambda: dict(
            data=members(file_manager), column_definitions=MEMBER_COLUMNS,
            title='Feature groups per run', index_field='MemberIndex',
            filters={'consensus': 'ConsensusIndex'}, go_to_fields=[],
            reclick=RECLICK, pagination=False,
        ),
    )

    with defer_reruns():
        left, right = st.columns(2)
        with left:
            st.markdown(f'Right: higher in **{compared}**. Left: higher in **{reference}**.')
            fc = st.slider(
                'Minimum change (|log2 fold change|)', 0.5, 3.0, FC_CUTOFF, 0.1,
                help='1 = at least twice (or half) as abundant.', key=f'{key}_fc',
            )
            padj = st.slider(
                'Maximum adjusted p-value', 0.001, 0.1, P_CUTOFF, 0.001,
                help='0.05 means about 5 % of the proteoforms called significant are expected to be false positives.',
                key=f'{key}_p',
            )
            volcano = _component(
                VolcanoPlot, 'compare_volcano', file_manager, DATASET, ['statistics'],
                lambda: dict(
                    data=_scan(file_manager, 'statistics').drop_nulls(['log2FC', 'p-adj']),
                    log2fc_column='log2FC', pvalue_column='p-adj', label_column='Proteoform',
                    interactivity={'consensus': 'ConsensusIndex'},
                    title=f'{compared} vs. {reference}',
                    up_color='#E74C3C', down_color='#3498DB', ns_color='#95A5A6',
                    show_threshold_lines=True, threshold_line_style='dash',
                    x_label='log2 fold change', y_label='-log10(adjusted p-value)',
                ),
            )
            volcano(
                key=f'{key}_volcano', state_manager=state_manager, height=500,
                fc_threshold=fc, p_threshold=padj, max_labels=10,
            )
            up, down = count_significant(statistics, fc, padj)
            st.markdown(f'**Higher in {compared}:** {up} | **Lower in {compared}:** {down}')
        with right:
            profile(key=f'{key}_profile', state_manager=state_manager, height=500)
            st.caption('Click a proteoform in the volcano plot or the table to show it here.')
        table(key=f'{key}_table', state_manager=state_manager, height=TABLE_HEIGHT)
        st.caption("The table's download button saves it as CSV, in the current sort and filter.")
        member_table(key=f'{key}_members', state_manager=state_manager, height=260)


def _render_pca(file_manager, state_manager, key, run_conditions, colors):
    runs = list(run_conditions)
    st.markdown(
        'Do runs of the same condition look alike, and do the conditions differ? Each dot is one '
        'run, placed so that runs with similar proteoform profiles sit close together.'
    )
    if len(runs) < 3:
        st.info('PCA needs at least three runs.')
        return
    wide = _scan(file_manager, 'normalized_wide').select('Proteoform', *runs).collect()
    available = wide.height
    if available <= 20:
        top_n = available
    else:
        top_n = st.slider(
            'Proteoforms used (most variable first)', 20, min(5000, available),
            min(500, available), 10, key=f'{key}_pca_top_n',
            help='PCA is computed only on the proteoforms with the highest variance across runs.',
        )
    variance = wide.select(pl.concat_list(runs).list.var().alias('variance'))['variance']
    top = wide.with_columns(variance.alias('_variance')).sort('_variance', descending=True).head(top_n).drop('_variance')
    try:
        pca = _component(
            PCAPlot, 'compare_pca', file_manager, DATASET, ['normalized_wide', 'design'],
            lambda: dict(
                data=top.lazy(), metadata=_metadata(run_conditions), sample_id_field='sample_id',
                group_field='group', n_components=min(5, len(runs)), title='Run PCA',
                group_colors=colors,
            ),
            extra=top_n,
        )
    except ValueError as e:
        st.error(f'PCA computation failed: {e}')
        return
    pc_columns = pca.get_pc_columns()
    x_column, y_column = st.columns(2)
    pc_x = x_column.selectbox('X-axis component', pc_columns, index=0, key=f'{key}_pc_x')
    pc_y = y_column.selectbox(
        'Y-axis component', pc_columns, index=1 if len(pc_columns) > 1 else 0, key=f'{key}_pc_y',
    )
    pca(
        key=f'{key}_pca', state_manager=state_manager, height=550,
        pc_x=int(pc_x.replace('PC', '')), pc_y=int(pc_y.replace('PC', '')),
    )
    st.markdown(
        '**Explained variance:** '
        + ', '.join(f'{c} {r * 100:.1f}%' for c, r in zip(pc_columns, pca.get_variance_ratio()))
    )
    st.caption(
        'Runs of one condition should cluster together and conditions should separate; an '
        'isolated run is a candidate outlier worth checking before trusting the statistics.'
    )


def _render_diagnostics(file_manager, runs, run_conditions, statistics, used):
    st.markdown('The plots behind each step of the comparison, as on OpenDDA\'s Fine-tune pages.')
    filtered = _scan(file_manager, 'filtered').collect().to_pandas()
    imputed = _scan(file_manager, 'imputed').collect().to_pandas()
    normalized = _scan(file_manager, 'normalized_wide').collect().to_pandas()
    show_fig(missingness_vs_intensity(filtered, runs), 'compare-missing-values')
    st.caption(
        'If proteoforms with missing values sit clearly lower than complete ones, values are missing '
        'because they fall below the detection limit (fill them with the smallest value, MNAR).'
    )
    show_fig(imputed_value_preview(filtered, imputed, runs), 'compare-imputed-values')
    show_fig(
        sample_distributions(
            imputed, normalized, runs, run_conditions,
            after_is_log=used['postproc-transform'] != 'None',
        ),
        'compare-run-distributions',
    )
    st.caption('Well-normalized runs have aligned medians; a run that stays offset is worth checking.')
    show_fig(pvalue_histogram(statistics.to_pandas()), 'compare-p-values')
    st.caption(
        'A flat histogram with a peak near 0 is what a sound test looks like; a U-shape or a peak '
        'near 1 suggests the scale or the test does not fit the data.'
    )
