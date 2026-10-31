"""Views of a comparison of FLASHQuant runs (see src/parse/quant_compare.py).

The comparison is stored as one dataset of its own cache: the linked
proteoforms (`consensus`), their feature groups per run (`members`), the
test results (`statistics`), normalized abundances per run (`normalized`),
the z-scores of the most significant proteoforms (`zscores`) and the
experimental design (`design`). All views share the identifier
`consensus` (a linked proteoform's ConsensusIndex): clicking a point in the
volcano plot or a row in the table shows its abundances and its feature
groups in every run.
"""

import polars as pl
import streamlit as st

from openms_insight import ClusteredHeatmap, SeriesPlot, Table, VolcanoPlot, defer_reruns

from src.render import insight_data
from src.render.render import RECLICK, TABLE_HEIGHT, MissingResult, _component, _fixed

DATASET = 'comparison'
RESULTS = ['consensus', 'members', 'statistics', 'normalized', 'zscores', 'design']

CONDITION_COLORS = ['#3366CC', '#DC3912', '#FF9900', '#109618', '#990099', '#0099C6']

_P = {'formatter': 'scientific', 'formatterParams': {'precision': 3}}
_FC = {'formatter': 'signed', 'formatterParams': {'precision': 2}}

CONSENSUS_COLUMNS = [
    {'title': 'Proteoform', 'field': 'Proteoform', 'sorter': 'string',
     'headerTooltip': 'Monoisotopic mass and apex retention time of the linked feature groups.'},
    {'title': 'Monoisotopic Mass', 'field': 'MonoisotopicMass', 'sorter': 'number', **_fixed()},
    {'title': 'Retention Time (min)', 'field': 'RTMinutes', 'sorter': 'number', **_fixed(2)},
    {'title': 'Runs', 'field': 'Runs', 'sorter': 'number',
     'headerTooltip': 'In how many runs the proteoform was quantified.'},
    {'title': 'log2 Fold Change', 'field': 'log2FC', 'sorter': 'number', **_FC},
    {'title': 'p-value', 'field': 'p-value', 'sorter': 'number', **_P},
    {'title': 'Adjusted p-value', 'field': 'p-adj', 'sorter': 'number', **_P,
     'headerTooltip': 'Benjamini-Hochberg. Empty: quantified in too few runs to be tested.'},
]

MEMBER_COLUMNS = [
    {'title': 'Run', 'field': 'Run', 'sorter': 'string'},
    {'title': 'Feature Group', 'field': 'FeatureGroupIndex', 'sorter': 'number',
     'headerTooltip': 'Index of the feature group in the run (FLASHQuant Viewer).'},
    {'title': 'Monoisotopic Mass', 'field': 'MonoisotopicMass', 'sorter': 'number', **_fixed()},
    {'title': 'Retention Time (min)', 'field': 'RTMinutes', 'sorter': 'number', **_fixed(2)},
    {'title': 'Quantity', 'field': 'Quantity', 'sorter': 'number', **_fixed(0)},
]


def conditions(file_manager):
    """Conditions of the stored design, in the order they were tested."""
    design = insight_data.scan(file_manager, DATASET, 'design').collect()
    return sorted(design['Condition'].unique().to_list())


def condition_colors(names):
    return {name: CONDITION_COLORS[i % len(CONDITION_COLORS)] for i, name in enumerate(names)}


def consensus_table(file_manager):
    """Linked proteoforms with their test results (null where not tested)."""
    consensus = insight_data.scan(file_manager, DATASET, 'consensus').select(
        'ConsensusIndex', 'Proteoform', 'MonoisotopicMass', 'RT', 'Runs',
    )
    statistics = insight_data.scan(file_manager, DATASET, 'statistics')
    tests = [c for c in ('log2FC', 'p-value', 'p-adj') if c in statistics.collect_schema().names()]
    return (
        consensus.join(statistics.select('ConsensusIndex', *tests), on='ConsensusIndex', how='left')
        .with_columns((pl.col('RT') / 60).alias('RTMinutes'))
        .drop('RT')
        .sort('p-adj', nulls_last=True)
    )


def members(file_manager):
    return (
        insight_data.scan(file_manager, DATASET, 'members')
        .with_columns((pl.col('RT') / 60).alias('RTMinutes'))
        .drop('RT')
        .with_row_index('MemberIndex')
        .with_columns(pl.col('MemberIndex').cast(pl.Int64))
    )


def render_comparison(file_manager, state_manager, key='flashquant_compare'):
    for name in RESULTS:
        if insight_data.result_path(file_manager, DATASET, name) is None:
            raise MissingResult(name)
    names = conditions(file_manager)
    colors = condition_colors(names)
    two_conditions = len(names) == 2
    columns = [
        c for c in CONSENSUS_COLUMNS if two_conditions or c['field'] != 'log2FC'
    ]

    table = _component(
        Table, 'compare_table', file_manager, DATASET, ['consensus', 'statistics'],
        lambda: dict(
            data=consensus_table(file_manager), column_definitions=columns,
            title='Linked proteoforms', index_field='ConsensusIndex',
            interactivity={'consensus': 'ConsensusIndex'}, go_to_fields=['Proteoform'],
            reclick=RECLICK,
        ),
    )
    profile = _component(
        SeriesPlot, 'compare_profile', file_manager, DATASET, ['normalized'],
        lambda: dict(
            data=insight_data.scan(file_manager, DATASET, 'normalized'),
            filters={'consensus': 'ConsensusIndex'},
            x_column='Run', y_column='Abundance', series_column='Condition',
            series_colors=colors, mode='markers', y_range_mode='normal',
            title='Abundance per run', x_label='Run',
            y_label='log2 abundance (median normalized)',
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
    design = insight_data.scan(file_manager, DATASET, 'design').collect()
    heatmap = _component(
        ClusteredHeatmap, 'compare_heatmap', file_manager, DATASET, ['zscores', 'design'],
        lambda: dict(
            data=insight_data.scan(file_manager, DATASET, 'zscores'), id_col='Proteoform',
            metadata=design.rename({'Run': 'sample_id', 'Condition': 'group'}),
            title='Most significant proteoforms (z-scores, clustered)',
            x_label='Runs', y_label='Proteoforms',
            colorscale=[[0, '#6699E0'], [0.5, '#FFFFFF'], [1, '#E06666']],
            reversescale=False, intensity_label='Z-score', group_colors=colors,
        ),
    )

    with defer_reruns():
        left, right = st.columns(2)
        with left:
            if two_conditions:
                volcano = _component(
                    VolcanoPlot, 'compare_volcano', file_manager, DATASET, ['statistics'],
                    lambda: dict(
                        data=insight_data.scan(file_manager, DATASET, 'statistics')
                        .drop_nulls(['log2FC', 'p-adj']),
                        log2fc_column='log2FC', pvalue_column='p-adj',
                        label_column='Proteoform', interactivity={'consensus': 'ConsensusIndex'},
                        title=f'{names[1]} vs. {names[0]}',
                        y_label='-log10(adjusted p-value)',
                    ),
                )
                fc_column, p_column = st.columns(2)
                fc = fc_column.number_input(
                    'Fold change threshold (|log2|)', 0.0, 10.0, 1.0, 0.1, key=f'{key}_fc',
                )
                p = p_column.number_input(
                    'Adjusted p-value threshold', 0.0001, 1.0, 0.05, 0.01, format='%.4f',
                    key=f'{key}_p',
                )
                volcano(key=f'{key}_volcano', state_manager=state_manager, height=450,
                        fc_threshold=fc, p_threshold=p)
            else:
                heatmap(key=f'{key}_heatmap', state_manager=state_manager, height=450)
        with right:
            profile(key=f'{key}_profile', state_manager=state_manager, height=450)
        table(key=f'{key}_table', state_manager=state_manager, height=TABLE_HEIGHT)
        member_table(key=f'{key}_members', state_manager=state_manager, height=260)
        if two_conditions:
            heatmap(key=f'{key}_heatmap', state_manager=state_manager, height=700)
