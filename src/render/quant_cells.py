"""The Quantification components of the Viewer: one function per Layout
Manager cell, each drawing exactly one component of a comparison of runs
(src/quant/session.py) with no page layout of its own.

Every builder has the signature of the shared quantification views,
``build_<view>(session, state_manager, key, height)``, so swapping them for
the shared ones later only changes imports. ``settings_form(session,
container)`` draws the design and analysis settings into a page or a dialog.

The cells of one comparison share the identifier ``consensus`` (the row of an
entity in the level shown): clicking a proteoform in the table or the volcano
plot shows it in the profile and the per-run table.
"""

import polars as pl
import streamlit as st

from openms_insight import ClusteredHeatmap, PCAPlot, SeriesPlot, Table, VolcanoPlot

from src.common.common import show_fig
from src.parse import differential as diff
from src.quant.session import LEVEL_NAMES, RESULTS, STATUS_LABELS
from src.render import insight_data
from src.render.postprocessing_plots import (
    imputed_value_preview, missingness_vs_intensity, pvalue_histogram, sample_distributions,
)
from src.render.render import RECLICK, MissingResult, _component, _fixed


CONDITION_COLORS = ['#3366CC', '#DC3912', '#FF9900', '#109618', '#990099', '#0099C6']

# OpenDDA's cut-offs for calling a change
FC_CUTOFF = 1.0
P_CUTOFF = 0.05

_P = {'formatter': 'scientific', 'formatterParams': {'precision': 3}}
_FC = {'formatter': 'signed', 'formatterParams': {'precision': 2}}


def comparison_groups(conditions):
    """(reference, compared): a positive log2FC means higher in the compared
    condition (openms_insight subtracts the first condition's mean from the
    last one's, in sorted order)."""
    names = sorted(set(conditions))
    return (names[0], names[-1]) if len(names) >= 2 else None


def condition_colors(names):
    return {name: CONDITION_COLORS[i % len(CONDITION_COLORS)] for i, name in enumerate(sorted(names))}

STATUS_COLORS = {
    'detected': '#2E7D32', 'transferred, identified': '#1565C0',
    'transferred': '#EF6C00', 'missing': '#9E9E9E',
}
NOUN = {'feature': 'feature', 'proteoform': 'proteoform', 'protein': 'protein'}


def _noun(session, plural=False):
    noun = NOUN[session.level]
    return noun + 's' if plural else noun


def _cutoffs(session):
    settings = session.settings
    return float(settings.get('fc-cutoff', FC_CUTOFF)), float(settings.get('padj-cutoff', P_CUTOFF))


def _ready(session):
    """Run the comparison chain if needed; False (with a message) when the
    comparison cannot run yet."""
    problem = session.analyse()
    if problem:
        st.info(problem)
        return False
    return True


def _sources(session, *names):
    return [session.result_name(n) for n in names] + ['quant_design', 'quant_settings']


def _change(compared, fc, padj):
    significant = (pl.col('p-adj') <= padj) & (pl.col('log2FC').abs() >= fc)
    return (
        pl.when(significant & (pl.col('log2FC') > 0)).then(pl.lit(f'higher in {compared}'))
        .when(significant & (pl.col('log2FC') < 0)).then(pl.lit(f'lower in {compared}'))
        .otherwise(pl.lit('–'))
    )


def _sample(session):
    """'Sample' for labelled samples (TMT, iTRAQ), else 'Run'."""
    return session.sample_noun.capitalize()


def _table_columns(session, runs):
    level = session.level
    noun = _noun(session).capitalize()
    columns = [{'title': noun, 'field': 'Proteoform', 'sorter': 'string',
                'headerTooltip': {
                    'feature': 'Monoisotopic mass and apex retention time of the consensus feature.',
                    'proteoform': 'The protein, the residues the proteoform covers and its mass shifts. The ProForma is in its own column.',
                    'protein': 'The protein the proteoforms belong to.',
                }[level]}]
    if level == 'proteoform':
        columns.append({'title': 'Protein', 'field': 'protein_label', 'sorter': 'string'})
        names = session.annotations().collect_schema().names()
        if 'start' in names:
            columns.append({'title': 'Residues', 'field': 'residues', 'sorter': 'string',
                            'headerTooltip': 'First and last residue of the proteoform in the protein sequence.'})
        if 'shifts' in names:
            columns.append({'title': 'Mass shifts', 'field': 'shifts', 'sorter': 'string',
                            'headerTooltip': 'Mass shifts (Da) at a residue, or within a range when not localized.'})
        columns.append({'title': 'ProForma', 'field': 'proforma', 'sorter': 'string'})
    if level == 'protein':
        columns.append({'title': 'Description', 'field': 'description', 'sorter': 'string'})
    columns += [
        {'title': 'Change', 'field': 'Change', 'sorter': 'string'},
        {'title': 'log2 fold change', 'field': 'log2FC', 'sorter': 'number', **_FC},
        {'title': 'Adjusted p-value', 'field': 'p-adj', 'sorter': 'number', **_P,
         'headerTooltip': 'Empty: removed by the filter, so not tested.'},
    ]
    if level != 'protein':
        columns.append({'title': 'Mass (Da)', 'field': 'mass', 'sorter': 'number', **_fixed(2)})
    if level == 'feature':
        columns.append({'title': 'RT (min)', 'field': 'rt', 'sorter': 'number', **_fixed(2)})
        columns.append({'title': 'Identified as', 'field': 'parent', 'sorter': 'string',
                        'headerTooltip': 'The proteoform of the PrSMs assigned to this feature (FLASHTnT).'})
    if level == 'protein':
        columns.append({'title': 'Proteoforms', 'field': 'proteoforms', 'sorter': 'number'})
    if 'PrSMs' in session.annotations().collect_schema().names():
        columns.append({'title': 'PrSMs', 'field': 'PrSMs', 'sorter': 'number',
                        'headerTooltip': 'PrSMs (FLASHTnT, q-value ≤ 1%) assigned to it in all runs.'})
    columns.append({'title': 'Detected in', 'field': 'Detected', 'sorter': 'number',
                    'headerTooltip': 'Runs where it was detected; in the others its value was transferred (match between runs) or is missing.'})
    for run in runs:
        columns.append({'title': run, 'field': run, 'sorter': 'string',
                        'headerTooltip': 'detected; transferred, identified (match between runs, identified in the run); transferred (match between runs); missing'})
    return columns


def build_summary(session, state_manager, key, height=None):
    """Counts of the comparison, as metrics."""
    del state_manager, height
    if not _ready(session):
        return
    reference, compared = comparison_groups(session.conditions().values())
    fc, padj = _cutoffs(session)
    statistics = session.result('statistics').collect()
    tested = statistics.drop_nulls(['log2FC', 'p-adj'])
    significant = tested.filter(pl.col('p-adj') <= padj)
    up = significant.filter(pl.col('log2FC') >= fc).height
    down = significant.filter(pl.col('log2FC') <= -fc).height
    with st.container(key=f'{key}_summary', border=True):
        metrics = st.columns(4)
        metrics[0].metric(f'{_sample(session)}s compared', len(session.conditions()))
        metrics[1].metric(f'{_noun(session, True).capitalize()} tested', f'{tested.height:,}',
                          help=f'Of {statistics.height:,}; the others were removed by the filter.')
        metrics[2].metric(f'Higher in {compared}', f'{up:,}')
        metrics[3].metric(f'Lower in {compared}', f'{down:,}')
        st.caption(
            f'{compared} vs. {reference}: adjusted p-value ≤ {padj:g} and |log2 fold change| ≥ {fc:g}. '
            'Change the conditions and cut-offs with **Edit comparison**.'
        )


def _comparison_table(session):
    _, compared = comparison_groups(session.conditions().values())
    fc, padj = _cutoffs(session)
    table = session.table()
    if session.level == 'proteoform':
        table = table.with_columns(
            pl.col('accession').str.split('|').list.get(-1, null_on_oob=True).alias('protein_label'),
        )
        if 'start' in table.columns:
            table = table.with_columns(
                (pl.col('start').cast(pl.Utf8) + '\u2013' + pl.col('end').cast(pl.Utf8)).alias('residues'),
            )
    return table.with_columns(_change(compared, fc, padj).alias('Change')).sort('p-adj', nulls_last=True).lazy()


def build_table(session, state_manager, key, height):
    """Every entity of the level, most significant first, with its status in
    every run."""
    if not _ready(session):
        return
    runs = session.runs
    _component(
        Table, f'quant_table_{session.level}', session.file_manager, session.dataset_id,
        _sources(session, 'statistics') + [f'quant_{session.level}', f'quant_{session.level}_annotations'],
        lambda: dict(
            data=_comparison_table(session), column_definitions=_table_columns(session, runs),
            title=f'{LEVEL_NAMES[session.level]} of the compared {_sample(session).lower()}s', index_field='ConsensusIndex',
            interactivity={'consensus': 'ConsensusIndex'}, go_to_fields=['Proteoform'],
            reclick=RECLICK,
        ),
    )(key=key, state_manager=state_manager, height=height)


def build_volcano(session, state_manager, key, height):
    if not _ready(session):
        return
    reference, compared = comparison_groups(session.conditions().values())
    fc, padj = _cutoffs(session)
    _component(
        VolcanoPlot, f'quant_volcano_{session.level}', session.file_manager, session.dataset_id,
        _sources(session, 'statistics'),
        lambda: dict(
            data=session.result('statistics').drop_nulls(['log2FC', 'p-adj']),
            log2fc_column='log2FC', pvalue_column='p-adj', label_column='Proteoform',
            interactivity={'consensus': 'ConsensusIndex'},
            title=f'{compared} vs. {reference}',
            up_color='#E74C3C', down_color='#3498DB', ns_color='#95A5A6',
            show_threshold_lines=True, threshold_line_style='dash',
            x_label='log2 fold change', y_label='-log10(adjusted p-value)',
        ),
    )(key=key, state_manager=state_manager, height=height, fc_threshold=fc, p_threshold=padj, max_labels=10)


def build_profile(session, state_manager, key, height):
    """The selected entity's normalized abundance in every run (sample), by condition."""
    if not _ready(session):
        return
    colors = condition_colors(set(session.conditions().values()))
    settings = session.settings
    _component(
        SeriesPlot, f'quant_profile_{session.level}', session.file_manager, session.dataset_id,
        _sources(session, 'normalized'),
        lambda: dict(
            data=session.result('normalized').with_columns(
                (pl.col('Run') + ': ' + pl.col('Status')).alias('Hover'),
            ),
            filters={'consensus': 'ConsensusIndex'},
            x_column='Run', y_column='Abundance', series_column='Condition',
            series_colors=colors, mode='markers', y_range_mode='normal', hover_column='Hover',
            title=f'Abundance per {_sample(session).lower()}', x_label=_sample(session),
            y_label=f"Abundance ({settings['postproc-transform']}, normalized)",
        ),
    )(key=key, state_manager=state_manager, height=height)


def _runs_frame(session):
    """The selected entity in every run: raw abundance, status, match score
    and q-value (the level's own values; summed for proteoforms and
    proteins)."""
    statistics = session.result('statistics').select('ConsensusIndex', 'id')
    conditions = session.conditions()
    return (
        session.long().join(statistics, on='id', how='inner')
        .with_columns(
            pl.col('status').replace_strict(STATUS_LABELS, default='missing').alias('Status'),
            pl.col('sample').replace_strict(conditions, default='', return_dtype=pl.Utf8).alias('Condition'),
        )
        .with_row_index('RowIndex')
        .with_columns(pl.col('RowIndex').cast(pl.Int64))
        .select('RowIndex', 'ConsensusIndex', pl.col('sample').alias('Run'), 'Condition',
                pl.col('intensity').alias('Abundance'), 'Status',
                pl.col('score').alias('MatchScore'), pl.col('qvalue').alias('MatchQvalue'))
    )


def build_runs(session, state_manager, key, height):
    """The selected entity in every run, with how its value was obtained."""
    if not _ready(session):
        return
    _component(
        Table, f'quant_runs_{session.level}', session.file_manager, session.dataset_id,
        _sources(session, 'statistics') + [f'quant_{session.level}'],
        lambda: dict(
            data=_runs_frame(session),
            column_definitions=[
                {'title': _sample(session), 'field': 'Run', 'sorter': 'string'},
                {'title': 'Condition', 'field': 'Condition', 'sorter': 'string'},
                {'title': 'Abundance', 'field': 'Abundance', 'sorter': 'number',
                 'formatter': 'scientific', 'formatterParams': {'precision': 3},
                 'headerTooltip': 'Raw signal (intensity × seconds), before normalization; summed over the features of a proteoform or protein.'},
                {'title': 'Status', 'field': 'Status', 'sorter': 'string',
                 'headerTooltip': 'detected: FLASHDeconv found it in this run. transferred: measured where another run found it (match between runs, FDR 1%); identified: and a PrSM of this run confirms it.'},
                {'title': 'Match score', 'field': 'MatchScore', 'sorter': 'number', **_fixed(3)},
                {'title': 'Match q-value', 'field': 'MatchQvalue', 'sorter': 'number', **_P},
            ],
            title=f'Selected {_noun(session)} per {_sample(session).lower()}', index_field='RowIndex',
            filters={'consensus': 'ConsensusIndex'}, go_to_fields=[], reclick=RECLICK, pagination=False,
        ),
    )(key=key, state_manager=state_manager, height=height)


def _metadata(conditions):
    return pl.DataFrame(
        {'sample_id': list(conditions), 'group': list(conditions.values())},
        schema={'sample_id': pl.String, 'group': pl.String},
    )


def build_pca(session, state_manager, key, height):
    if not _ready(session):
        return
    conditions = session.conditions()
    runs = list(conditions)
    if len(runs) < 3:
        st.info(f'PCA needs at least three {_sample(session).lower()}s with a condition.')
        return
    colors = condition_colors(set(conditions.values()))
    wide = session.result('normalized_wide').select('Proteoform', *runs).collect()
    variance = wide.select(pl.concat_list(runs).list.var().alias('v'))['v']
    top = wide.with_columns(variance.alias('_v')).sort('_v', descending=True).head(500).drop('_v')
    try:
        pca = _component(
            PCAPlot, f'quant_pca_{session.level}', session.file_manager, session.dataset_id,
            _sources(session, 'normalized_wide'),
            lambda: dict(
                data=top.lazy(), metadata=_metadata(conditions), sample_id_field='sample_id',
                group_field='group', n_components=min(5, len(runs)), title='Sample similarity (PCA)',
                group_colors=colors,
            ),
        )
    except ValueError as e:
        st.info(f'PCA could not be computed: {e}')
        return
    pca(key=key, state_manager=state_manager, height=height, pc_x=1, pc_y=2)
    ratios = pca.get_variance_ratio()
    st.caption('Explained variance: ' + ', '.join(
        f'PC{i + 1} {r * 100:.1f}%' for i, r in enumerate(ratios[:3])
    ) + f'. The 500 most variable {_noun(session, True)}.')


def build_heatmap(session, state_manager, key, height):
    if not _ready(session):
        return
    conditions = session.conditions()
    colors = condition_colors(set(conditions.values()))
    _component(
        ClusteredHeatmap, f'quant_heatmap_{session.level}', session.file_manager, session.dataset_id,
        _sources(session, 'zscores'),
        lambda: dict(
            data=session.result('zscores'), id_col='Proteoform', metadata=_metadata(conditions),
            title=f'Most significant {_noun(session, True)} (z-scores, clustered)',
            x_label=f'{_sample(session)}s', y_label=LEVEL_NAMES[session.level],
            colorscale=[[0, '#6699E0'], [0.5, '#FFFFFF'], [1, '#E06666']],
            reversescale=False, intensity_label='Z-score', group_colors=colors,
        ),
    )(key=key, state_manager=state_manager, height=height)


def build_diagnostics(session, state_manager, key, height=None):
    """The plots behind each step of the comparison chain."""
    del state_manager, height
    if not _ready(session):
        return
    conditions = session.conditions()
    runs = list(conditions)
    filtered = session.result('filtered').collect().to_pandas()
    imputed = session.result('imputed').collect().to_pandas()
    normalized = session.result('normalized_wide').collect().to_pandas()
    statistics = session.result('statistics').collect().to_pandas()
    settings = session.settings
    view = st.segmented_control(
        'Step', ['Missing values', 'Imputed values', 'Normalization', 'p-values'],
        default='Missing values', key=f'{key}_step', label_visibility='collapsed',
    ) or 'Missing values'
    if view == 'Missing values':
        show_fig(missingness_vs_intensity(filtered, runs), f'{key}-missing-values')
    elif view == 'Imputed values':
        show_fig(imputed_value_preview(filtered, imputed, runs), f'{key}-imputed-values')
    elif view == 'Normalization':
        show_fig(sample_distributions(
            imputed, normalized, runs, conditions, after_is_log=settings['postproc-transform'] != 'None',
        ), f'{key}-run-distributions')
    else:
        show_fig(pvalue_histogram(statistics), f'{key}-p-values')


BUILDERS = {
    'quant_summary': (build_summary, None),
    'quant_table': (build_table, 376),
    'quant_volcano': (build_volcano, 480),
    'quant_profile': (build_profile, 480),
    'quant_runs': (build_runs, 300),
    'quant_pca': (build_pca, 480),
    'quant_heatmap': (build_heatmap, 600),
    'quant_diagnostics': (build_diagnostics, None),
}


def render_cell(comp_name, session, state_manager, key):
    builder, height = BUILDERS[comp_name]
    try:
        builder(session, state_manager, key, height)
    except MissingResult as e:
        st.info(f'This comparison has no data for this view ({e}).')


# -- settings ---------------------------------------------------------------

def settings_form(session, container=None):
    """The design (condition and replicate per run) and the analysis settings
    of a comparison. Returns True when they were saved."""
    used = session.settings
    with container if container is not None else st.container():
        st.markdown(f'**{_sample(session)}s and conditions**')
        st.caption(f'{_sample(session)}s with the same condition are replicates. A {_sample(session).lower()} without a condition is left out.')
        design = st.data_editor(
            session.design().to_pandas(),
            column_config={
                'sample': st.column_config.TextColumn(_sample(session), disabled=True),
                'condition': st.column_config.TextColumn('Condition'),
                'replicate': st.column_config.NumberColumn('Replicate', min_value=1, step=1),
            },
            hide_index=True, use_container_width=True, key=f'{session.dataset_id}_design_editor',
        )
        level = st.segmented_control(
            'Compare', [l for l in session.levels], default=session.level,
            format_func=lambda l: LEVEL_NAMES[l], key=f'{session.dataset_id}_level_choice',
            help='Proteoforms and proteins need identifications (FLASHTnT).',
        ) or session.level
        st.markdown('**Calling a change**')
        c1, c2 = st.columns(2)
        fc = c1.number_input('Minimum |log2 fold change|', 0.0, 5.0, float(used.get('fc-cutoff', FC_CUTOFF)), 0.1,
                             key=f'{session.dataset_id}_fc', help='1 = at least twice (or half) as abundant.')
        padj = c2.number_input('Maximum adjusted p-value', 0.0001, 0.5, float(used.get('padj-cutoff', P_CUTOFF)), 0.005,
                               format='%.4f', key=f'{session.dataset_id}_padj')
        chosen = {}
        with st.expander('Analysis steps (filter, missing values, normalization, test)'):
            def choice(label, options, setting, column=st, help=None):
                value = used[setting]
                return column.selectbox(label, options, index=options.index(value) if value in options else 0,
                                        key=f'{session.dataset_id}_{setting}', help=help)
            a, b = st.columns(2)
            method = choice('1. Filter by', diff.FILTER_METHODS, 'postproc-filter-method', a)
            chosen['postproc-filter-method'] = method
            setting, label = {
                'Low Repeatability': ('postproc-filter-max-missing', 'Max missing ratio (%)'),
                'Low Abundance': ('postproc-filter-abundance-percentile', 'Threshold percentile (%)'),
                'Low Variance': ('postproc-filter-variance-percentile', 'Threshold percentile (%)'),
            }[method]
            chosen[setting] = b.slider(label, 0.0, 100.0, float(used[setting]), 5.0, key=f'{session.dataset_id}_{setting}')
            a, b = st.columns(2)
            impute = choice('2. Why are values missing?', diff.IMPUTE_CLASSES, 'postproc-impute-class', a)
            chosen['postproc-impute-class'] = impute
            if impute.startswith('MAR'):
                chosen['postproc-impute-mar-strategy'] = choice('Fill with the condition', diff.MAR_STRATEGIES,
                                                                'postproc-impute-mar-strategy', b)
            else:
                chosen['postproc-impute-mnar-scope'] = choice('Smallest value of', diff.MNAR_SCOPES,
                                                              'postproc-impute-mnar-scope', b)
            a, b, c = st.columns(3)
            chosen['postproc-transform'] = choice('3. Transformation', diff.TRANSFORMS, 'postproc-transform', a)
            chosen['postproc-normalization'] = choice('Normalization', diff.NORMALIZATIONS, 'postproc-normalization', b)
            chosen['postproc-scaling'] = choice('Scaling', diff.SCALINGS, 'postproc-scaling', c)
            a, b = st.columns(2)
            n_conditions = max(len({c for c in design['condition'] if str(c).strip()}), 2)
            chosen['postproc-stat-test'] = choice('4. Test', diff.test_options(n_conditions), 'postproc-stat-test', a,
                                                  help='limma_like: moderated t-test (recommended).')
            chosen['postproc-fdr'] = choice('Correction', diff.FDR_METHODS, 'postproc-fdr', b)
        if st.button('Apply', type='primary', key=f'{session.dataset_id}_apply'):
            frame = pl.from_pandas(design).with_columns(
                pl.col('condition').fill_null('').cast(pl.Utf8).str.strip_chars(),
                pl.col('replicate').fill_null(1).cast(pl.Int64),
            )
            session.save_design(frame)
            session.save_settings({**used, **chosen, 'level': level, 'fc-cutoff': fc, 'padj-cutoff': padj})
            return True
    return False
