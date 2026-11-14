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

from openms_insight import ClusteredHeatmap, Heatmap, PCAPlot, SeriesPlot, Table, VolcanoPlot

from src.parse import differential as diff
from src.quant import qc, sources
from src.quant.session import LEVEL_NAMES, RESULTS, STATUS_LABELS
from src.render import comparison_plots, go_terms, insight_data
from src.render.postprocessing_plots import (
    IMPUTED_COLORS, IMPUTED_SERIES, MISSING_COLORS, MISSING_SERIES, STAGE_COLORS, STAGE_SERIES,
    imputed_value_preview, missingness_vs_intensity, pvalue_histogram, sample_distributions,
)
from src.render.render import RECLICK, MissingResult, _component, _fixed
from src.render.series import show_series


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


def _plot(session, name, spec, key, state_manager, height):
    """Draw a plot of src/render/comparison_plots.py (its settings and ``data``)."""
    show_series(spec.pop('data'), name, session.file_manager.cache_path, key, height=height,
                state_manager=state_manager, **spec)


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
    """(column definitions, fields hidden until the viewer turns them on).

    Shown by default: the entity, what changed (change, log2 fold change,
    adjusted p-value), the identification of a feature when there is one and
    in how many compared samples it was detected. The rest (mass, RT,
    sequence details, PrSMs, the status in each sample) is in the column
    menu and in the row's overview."""
    level = session.level
    noun = _noun(session).capitalize()
    names = session.annotations().collect_schema().names()
    hidden = []
    columns = [{'title': noun, 'field': 'Proteoform', 'sorter': 'string',
                'headerTooltip': {
                    'feature': 'Monoisotopic mass and apex retention time of the consensus feature.',
                    'proteoform': 'The protein, the residues the proteoform covers and its mass shifts. The ProForma is in its own column.',
                    'protein': 'The protein the proteoforms belong to.',
                }[level]}]
    if level == 'proteoform':
        columns.append({'title': 'Protein', 'field': 'protein_label', 'sorter': 'string'})
        if 'start' in names:
            columns.append({'title': 'Residues', 'field': 'residues', 'sorter': 'string',
                            'headerTooltip': 'First and last residue of the proteoform in the protein sequence.'})
            hidden.append('residues')
        if 'shifts' in names:
            columns.append({'title': 'Mass shifts', 'field': 'shifts', 'sorter': 'string',
                            'headerTooltip': 'Mass shifts (Da) at a residue, or within a range when not localized.'})
            hidden.append('shifts')
        columns.append({'title': 'ProForma', 'field': 'proforma', 'sorter': 'string'})
        hidden.append('proforma')
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
        hidden.append('mass')
    if level == 'feature':
        columns.append({'title': 'RT (min)', 'field': 'rt', 'sorter': 'number', **_fixed(2)})
        hidden.append('rt')
        columns.append({'title': 'Identified as', 'field': 'parent', 'sorter': 'string',
                        'headerTooltip': 'The proteoform of the PrSMs assigned to this feature (FLASHTnT).'})
        # Empty in every row of a comparison without identification
        if 'parent' not in names or session.annotations().select(pl.col('parent').is_not_null().any()).collect().item() is not True:
            hidden.append('parent')
    if level == 'protein':
        columns.append({'title': 'Proteoforms', 'field': 'proteoforms', 'sorter': 'number'})
        hidden.append('proteoforms')
    if 'PrSMs' in names:
        columns.append({'title': 'PrSMs', 'field': 'PrSMs', 'sorter': 'number',
                        'headerTooltip': 'PrSMs (FLASHTnT, q-value ≤ 1%) assigned to it in all runs.'})
        hidden.append('PrSMs')
    sample = _sample(session).lower()
    columns.append({
        'title': 'Detected in', 'field': 'Detected', 'sorter': 'number',
        # "6/8": detected in 6 of the 8 compared samples
        'formatter': 'money', 'formatterParams': {'symbol': f'/{len(runs)}', 'symbolAfter': True, 'precision': False},
        'headerTooltip': f'In how many of the {len(runs)} compared {sample}s it was detected; in the others its '
                         'value was transferred (match between runs) or is missing. '
                         f'Each {sample} has a column of its own in the column menu.',
    })
    for run in runs:
        columns.append({'title': run, 'field': run, 'sorter': 'string',
                        'headerTooltip': 'detected; transferred, identified (match between runs, identified in the run); transferred (match between runs); missing'})
        hidden.append(run)
    return columns, hidden


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
    """Every entity of the level, most significant first: what changed by
    default, the rest (its status in every compared sample among it) in the
    column menu and the row's overview."""
    if not _ready(session):
        return
    columns, hidden = _table_columns(session, session.compared)
    _component(
        Table, f'quant_table_{session.level}', session.file_manager, session.dataset_id,
        _sources(session, 'statistics') + [f'quant_{session.level}', f'quant_{session.level}_annotations'],
        lambda: dict(
            # Hidden as `visible: False` in the definitions, which the cache
            # keeps: _component loads a cached table with its cache_id only,
            # so a hidden_columns argument would apply to the first render only
            data=_comparison_table(session),
            column_definitions=[{**c, 'visible': False} if c['field'] in hidden else c for c in columns],
            # The hidden columns of a row, the status in each sample among them
            row_overview='{Proteoform}',
            title=f'{LEVEL_NAMES[session.level]} of the compared {_sample(session).lower()}s', index_field='ConsensusIndex',
            interactivity={'consensus': 'ConsensusIndex'}, go_to_fields=['Proteoform'],
            # Check boxes pick entries for the GO enrichment cell
            multi_select={go_terms.COMPARISON_SET: 'Proteoform'},
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
            # A box or lasso, or All changed / Higher / Lower, picks entries
            # for the GO enrichment cell
            multi_select={go_terms.COMPARISON_SET: 'Proteoform'},
            selection_presets=True,
            preset_labels={'up': f'Higher in {compared}', 'down': f'Lower in {compared}'},
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
            intensity_label='Z-score', group_colors=colors,
            multi_select={go_terms.COMPARISON_SET: 'Proteoform'},
        ),
    )(key=key, state_manager=state_manager, height=height)


def _frame(session, name):
    return session.result(name).collect().to_pandas()


# -- identification ------------------------------------------------------------
# What each run identified: these need no conditions, so they show as soon as
# the runs are compared. They count proteoforms when the runs were identified
# (the level shown otherwise).

def _id_level(session):
    return qc.identification_level(session.levels, session.level)


def _id_noun(session, level):
    if level == 'feature' and 'proteoform' not in session.levels:
        return 'features'
    return LEVEL_NAMES[level].lower()



def build_id_counts(session, state_manager, key, height=None):
    """Entities each run detected itself and those matched from other runs."""
    level = _id_level(session)
    counts = qc.counts_per_run(session.long(level).collect(), session.runs)
    _plot(session, f'quant-id-counts-{level}',
          comparison_plots.counts_per_run(counts, _id_noun(session, level), _sample(session)),
          key, state_manager, height)


def build_id_frequency(session, state_manager, key, height=None):
    """How many entities were found in 1, 2, ... all runs."""
    level = _id_level(session)
    table = qc.frequency(session.long(level).collect(), session.runs)
    _plot(session, f'quant-id-frequency-{level}',
          comparison_plots.frequency(table, _id_noun(session, level), _sample(session)), key, state_manager, height)


def build_id_overlap(session, state_manager, key, height):
    """What every two runs share (Jaccard index of the detected entities),
    clustered so that runs alike sit together."""
    level = _id_level(session)
    runs = session.runs
    if len(runs) < 2:
        st.info(f'Shared identifications need at least two {_sample(session).lower()}s.')
        return
    conditions = session.conditions()
    _component(
        ClusteredHeatmap, f'quant_id_overlap_{level}', session.file_manager, session.dataset_id,
        [f'quant_{level}', 'quant_design', 'quant_settings'],
        lambda: dict(
            data=qc.overlap(session.long(level).collect(), runs).lazy(), id_col='Sample',
            metadata=_metadata(conditions) if conditions else None,
            title=f'Shared {_id_noun(session, level)} (% of those either {_sample(session).lower()} detected)',
            x_label=f'{_sample(session)}s', y_label=f'{_sample(session)}s', intensity_label='Shared (%)',
            colorscale='Blues', reversescale=True,
            group_colors=condition_colors(set(conditions.values())) if conditions else None,
        ),
        required=[f'quant_{level}'],
    )(key=key, state_manager=state_manager, height=height)


def build_id_table(session, state_manager, key, height):
    """Every entity any run found, with how each run found it."""
    level = _id_level(session)
    runs = session.runs
    table = qc.identification_table(session.long(level).collect(), session.annotations(level).collect(), runs)
    noun = _id_noun(session, level)
    columns = [{'title': noun[:-1].capitalize() if noun.endswith('s') else noun, 'field': 'label', 'sorter': 'string'}]
    if 'protein_label' in table.columns and level == 'proteoform':
        columns.append({'title': 'Protein', 'field': 'protein_label', 'sorter': 'string'})
    if 'mass' in table.columns:
        columns.append({'title': 'Mass (Da)', 'field': 'mass', 'sorter': 'number', **_fixed(2)})
    if 'proteoforms' in table.columns:
        columns.append({'title': 'Proteoforms', 'field': 'proteoforms', 'sorter': 'number'})
    if 'PrSMs' in table.columns:
        columns.append({'title': 'PrSMs', 'field': 'PrSMs', 'sorter': 'number',
                        'headerTooltip': 'Spectra (PrSMs) that identify it, over all runs.'})
    columns += [
        {'title': 'Detected', 'field': 'Detected', 'sorter': 'number',
         'headerTooltip': f'{_sample(session)}s that found it themselves.'},
        {'title': 'Found', 'field': 'Found', 'sorter': 'number',
         'headerTooltip': f'{_sample(session)}s that found it themselves or matched it from another one.'},
    ]
    for run in runs:
        columns.append({'title': run, 'field': run, 'sorter': 'string',
                        'headerTooltip': 'detected; matched, identified (match between runs, identified in the run); matched (match between runs); missing'})
    _component(
        Table, f'quant_id_table_{level}', session.file_manager, session.dataset_id,
        [f'quant_{level}', f'quant_{level}_annotations', 'quant_design', 'quant_settings'],
        lambda: dict(
            data=table.lazy(), column_definitions=columns, index_field='RowIndex',
            title=f'{LEVEL_NAMES[level]} across the {_sample(session).lower()}s',
            go_to_fields=['label'], row_overview='{label}', interactivity={'quant_form': 'RowIndex'},
        ),
        required=[f'quant_{level}', f'quant_{level}_annotations'],
    )(key=key, state_manager=state_manager, height=height)


# -- quality control -------------------------------------------------------

def build_totals(session, state_manager, key, height=None):
    """Summed intensity of every run or labelled channel: unequal totals
    mean the samples were not loaded alike (normalization corrects for it)."""
    conditions = session.conditions()
    frame = qc.totals(session.long().collect(), session.runs)
    _plot(session, 'quant-totals', comparison_plots.totals(
        frame, condition_colors(set(conditions.values())), conditions, _sample(session)), key, state_manager, height)
    unlike = qc.unlike_replicates(frame, conditions)
    if unlike:
        st.warning(f'Under half or over twice the total of their replicates: {", ".join(unlike)}. '
                   'Check how they were loaded; normalization corrects for it.')


def build_correlation(session, state_manager, key, height):
    """Pearson correlation of the log intensities of every two samples,
    clustered: replicates should correlate best."""
    runs = session.runs
    if len(runs) < 2:
        st.info(f'A correlation needs at least two {_sample(session).lower()}s.')
        return
    conditions = session.conditions()
    level = session.level
    _component(
        ClusteredHeatmap, f'quant_correlation_{level}', session.file_manager, session.dataset_id,
        [f'quant_{level}', 'quant_design', 'quant_settings'],
        lambda: dict(
            data=qc.correlation(session.long(level).collect(), runs).lazy(), id_col='Sample',
            metadata=_metadata(conditions) if conditions else None,
            title=f'Correlation of the {_sample(session).lower()}s (log2 intensities of {_noun(session, True)})',
            x_label=f'{_sample(session)}s', y_label=f'{_sample(session)}s', intensity_label='Pearson r',
            colorscale='Viridis', reversescale=False,
            group_colors=condition_colors(set(conditions.values())) if conditions else None,
        ),
        required=[f'quant_{level}'],
    )(key=key, state_manager=state_manager, height=height)


def build_cv(session, state_manager, key, height=None):
    """CV of every entity within each condition, from the values measured
    (before imputation and normalization)."""
    if not _ready(session):
        return
    conditions = session.conditions()
    frame = qc.coefficients_of_variation(session.result('filtered').collect(), conditions)
    if frame.height == 0:
        st.info('A CV needs at least two replicates of a condition.')
        return
    _plot(session, 'quant-cv', comparison_plots.coefficients_of_variation(
        frame, condition_colors(set(conditions.values()))), key, state_manager, height)


def build_match_quality(session, state_manager, key, height=None):
    """How well the values matched between runs (match between runs) scored,
    against the matches that were rejected."""
    frame = qc.match_quality(session.long('feature').collect(), session.runs)
    if session.isobaric or frame.height == 0:
        st.info('No values were matched between runs: the samples were measured together '
                '(isobaric labels) or every run found every feature itself.')
        return
    _plot(session, 'quant-match-quality', comparison_plots.match_quality(frame, session.runs, _sample(session)),
          key, state_manager, height)


# The steps of the comparison chain, each a plot of its own in the layout.
# Saved layouts written before they were split name them together as
# 'quant_diagnostics' (quant_view.split_diagnostics).

def build_missing_values(session, state_manager, key, height=None):
    """Where missing values occur: mean intensity of complete entries and of
    entries with gaps."""
    if not _ready(session):
        return
    runs = list(session.conditions())
    data, gaps = missingness_vs_intensity(_frame(session, 'filtered'), runs)
    show_series(
        data, 'quant-missing-values', session.file_manager.cache_path, key, height=height,
        state_manager=state_manager, x_column='x', y_column='y', series_column='series',
        series_names=MISSING_SERIES, series_colors=MISSING_COLORS,
        title=f'Where missing values occur ({gaps} of {len(_frame(session, "filtered"))} have gaps)',
        x_label='Mean log2 intensity', y_label='Density',
    )


def build_imputed_values(session, state_manager, key, height=None):
    """The values imputed for the gaps against the measured ones."""
    if not _ready(session):
        return
    runs = list(session.conditions())
    show_series(
        imputed_value_preview(_frame(session, 'filtered'), _frame(session, 'imputed'), runs),
        'quant-imputed-values', session.file_manager.cache_path, key, height=height,
        state_manager=state_manager, x_column='x', y_column='y', series_column='series',
        series_names=IMPUTED_SERIES, series_colors=IMPUTED_COLORS,
        title='Observed values and the values imputation fills in', x_label='log2 intensity', y_label='Values',
    )


def _normalization_height(samples):
    """A box per sample, laid out horizontally."""
    return 44 * samples + 140


def build_normalization(session, state_manager, key, height=None):
    """Intensity distribution of each sample before and after normalization."""
    if not _ready(session):
        return
    samples = list(session.conditions())
    data = sample_distributions(
        _frame(session, 'imputed'), _frame(session, 'normalized_wide'), samples,
        after_is_log=session.settings['postproc-transform'] != 'None',
    )
    show_series(
        data, 'quant-run-distributions', session.file_manager.cache_path, key,
        height=max(height or 0, _normalization_height(len(samples))), state_manager=state_manager,
        x_column='Sample', y_column='Intensity', series_column='Stage', series_names=STAGE_SERIES,
        series_colors=STAGE_COLORS, mode='box', orientation='h',
        title='Intensity distributions before and after normalization', x_label=_sample(session),
        y_label='log2 intensity', y_range_mode='normal',
    )


def build_pvalues(session, state_manager, key, height=None):
    """Distribution of the p-values of the group comparison."""
    if not _ready(session):
        return
    data, expected = pvalue_histogram(_frame(session, 'statistics'))
    show_series(
        data, 'quant-p-values', session.file_manager.cache_path, key, height=height,
        state_manager=state_manager, x_column='p-value', y_column='Proteoforms', mode='bars', show_legend=False,
        title=f'p-value distribution ({expected:.0f} per bin if nothing changes)', x_label='p-value',
        y_label=_noun(session, True).capitalize(),
    )


def build_go(session, state_manager, key, height=None):
    """GO terms over-represented in the entries picked in the volcano plot,
    the table or the heatmap, against every tested entry."""
    del height
    if not _ready(session):
        return
    from src.render import go_view
    go_view.render_comparison(session, state_manager, key)


# -- from a proteoform to its PrSMs and their spectra ------------------------
# The identification table selects an entity (`quant_form`), the PrSM table
# lists the PrSMs of every run that identify it and selects one
# (`quant_prsm`), and the sequence view and spectrum of that PrSM are drawn
# from the run's own results with a state of their own.

def _prsm_frame(session):
    """The PrSMs of the comparison, with the RowIndex of the entity they
    identify in the identification table."""
    level = _id_level(session)
    if level == 'feature' and 'proteoform' not in session.levels:
        return level, None
    table = qc.identification_table(session.long(level).collect(), session.annotations(level).collect(),
                                    session.runs)
    prsms = sources.comparison_prsms(session.file_manager, session.dataset_id, level)
    frame = prsms.join(table.select(pl.col('id').alias('key'), 'RowIndex', 'label'), on='key', how='inner')
    # Which PrSMs can be shown: their run's results are in this workspace and
    # hold the spectrum. They come first, then the runs in their order.
    shown = []
    for run, run_ds in sources.run_datasets(session.file_manager, session.manifest).items():
        try:
            scans = insight_data.proteins(session.file_manager, run_ds).select('Scan').collect()['Scan']
        except Exception:  # a run without identification results
            continue
        shown.append(pl.DataFrame({'Run': run, 'Scan': scans.drop_nulls().cast(pl.Int64).unique()}))
    here = pl.concat(shown) if shown else pl.DataFrame(schema={'Run': pl.Utf8, 'Scan': pl.Int64})
    order = {run: i for i, run in enumerate(session.runs)}
    frame = (
        frame.join(here.with_columns(pl.lit('yes').alias('Spectrum')), on=['Run', 'Scan'], how='left')
        .with_columns(pl.col('Spectrum').fill_null('no'),
                      pl.col('Run').replace_strict(order, default=len(order), return_dtype=pl.Int64).alias('_run'))
        .sort(['RowIndex', 'Spectrum', '_run', 'Scan'], descending=[False, True, False, False]).drop('_run')
    )
    return level, frame


def _no_prsms_note(session):
    st.info('This comparison has no PrSMs: the runs were not identified, or were compared before the '
            'PrSMs were kept with the comparison. Compare identified runs to follow a proteoform to its spectra.')


def build_id_prsms(session, state_manager, key, height):
    """Every PrSM of the proteoform picked in the identification table, run by run."""
    level, frame = _prsm_frame(session)
    if frame is None or frame.height == 0:
        _no_prsms_note(session)
        return
    columns = [
        {'title': _sample(session), 'field': 'Run', 'sorter': 'string'},
        {'title': 'Scan', 'field': 'Scan', 'sorter': 'number'},
        {'title': 'RT (min)', 'field': 'RT', 'sorter': 'number', **_fixed(2)},
        {'title': 'Proteoform mass (Da)', 'field': 'ProteoformMass', 'sorter': 'number', **_fixed(2)},
        {'title': 'Precursor mass (Da)', 'field': 'PrecursorMass', 'sorter': 'number', **_fixed(2)},
        {'title': 'Error (ppm)', 'field': 'MassErrorPPM', 'sorter': 'number', **_fixed(2)},
        {'title': 'q-value', 'field': 'qvalue', 'sorter': 'number', **_P},
        {'title': 'Spectrum', 'field': 'Spectrum', 'sorter': 'string',
         'headerTooltip': "yes: the run's results are in this workspace, pick the PrSM to see its spectrum."},
        {'title': 'Proteoform', 'field': 'Proforma', 'sorter': 'string'},
    ]
    _component(
        Table, f'quant_id_prsms_{level}', session.file_manager, session.dataset_id,
        [f'quant_{level}', f'quant_{level}_annotations', 'identifications_tsv'],
        lambda: dict(
            data=frame.lazy(), column_definitions=columns, index_field='prsm_row',
            title='PrSMs of the selected identification', go_to_fields=['Run', 'Scan'],
            filters={'quant_form': 'RowIndex'}, interactivity={'quant_prsm': 'prsm_row'},
            reclick=RECLICK, reset_sort=True,
        ),
        extra=sorted(sources.run_datasets(session.file_manager, session.manifest).items()),
        required=[f'quant_{level}', 'identifications_tsv'],
    )(key=key, state_manager=state_manager, height=height)


def _run_state(session, run_ds, extra=''):
    from openms_insight import StateManager
    return StateManager(session_key=f'flash_quant_run_state_{session.dataset_id}_{run_ds}{extra}')


def _set(state, identifier, value):
    if state.get_selection(identifier) != value:
        state.set_selection(identifier, value)


def _selected_prsm(session, state_manager):
    """(run name, dataset id, state, scan) of the PrSM picked in the PrSM
    table, its run's state pointing at it; None (with a note) otherwise."""
    from src.render import render
    _, frame = _prsm_frame(session)
    if frame is None or frame.height == 0:
        _no_prsms_note(session)
        return None
    picked = state_manager.get_selection('quant_prsm')
    row = frame.filter(pl.col('prsm_row') == picked) if picked is not None else frame.head(0)
    form = state_manager.get_selection('quant_form')
    if row.height == 0 or (form is not None and row['RowIndex'][0] != form):
        st.info('Pick a PrSM in **PrSMs of the selected identification** to see its spectrum here.')
        return None
    prsm = row.row(0, named=True)
    run_ds = sources.run_datasets(session.file_manager, session.manifest).get(prsm['Run'])
    if run_ds is None:
        st.info(f"The results of {prsm['Run']} are not in this workspace, so its spectra cannot be shown.")
        return None
    hit = sources.prsm_in_run(render.insight_data.proteins(session.file_manager, run_ds).collect(), prsm['Scan'])
    if hit is None:
        st.info(f"Scan {prsm['Scan']} of {prsm['Run']} has no PrSM in the run's own results.")
        return None
    state = _run_state(session, run_ds)
    index, scan_idx = hit
    _set(state, 'protein', index)
    _set(state, 'scan', scan_idx)
    _set(state, 'pscan', scan_idx)
    st.caption(f"{prsm['label']}: scan {prsm['Scan']} of {prsm['Run']}")
    return run_ds, state


def _run_cell(comp_name, session, state_manager, key):
    from src.render import render
    picked = _selected_prsm(session, state_manager)
    if picked is None:
        return
    run_ds, state = picked
    render._render_cell(comp_name, session.file_manager, run_ds, 'topdown', state, key)


def build_prsm_sequence(session, state_manager, key, height=None):
    """Sequence view of the picked PrSM, from its run's results."""
    del height
    _run_cell('sequence_view', session, state_manager, key)


def build_prsm_spectrum(session, state_manager, key, height=None):
    """Annotated spectrum of the picked PrSM, from its run's results."""
    del height
    _run_cell('anno_spectrum', session, state_manager, key)


# -- reporter ions (isobaric labels) ----------------------------------------

def _fractions(session):
    """Fraction (run name) -> dataset id of the labelled runs of an isobaric
    comparison that are in this workspace; None with a note otherwise."""
    if not session.isobaric:
        st.info('This comparison has no reporter ions: its runs were not labelled (TMT, iTRAQ).')
        return None
    found = sources.run_datasets(session.file_manager, session.manifest)
    found = {name: ds for name, ds in found.items() if session.file_manager.result_exists(ds, 'isobaric_dfs')}
    if not found:
        st.info('The labelled runs of this comparison are not in this workspace, so their reporter ions '
                'cannot be shown.')
        return None
    return found


def _reporter_ions(session):
    fractions = _fractions(session)
    if fractions is None:
        return None
    frames = {}
    for name, ds in fractions.items():
        frame = session.file_manager.get_results(ds, ['isobaric_dfs'])['isobaric_dfs']
        frames[name] = frame if isinstance(frame, pl.DataFrame) else pl.from_pandas(frame)
    channels = list(session.manifest.get('channels') or [])
    return sources.reporter_ions(frames, channels, session.runs), channels


def _plex_label(session):
    from src.quant import isobaric
    return isobaric.label(session.manifest.get('plex', ''))


def build_reporter_intensities(session, state_manager, key, height=None):
    """Intensity of every label's reporter ions, and how often each label was found."""
    got = _reporter_ions(session)
    if got is None:
        return
    ions, channels = got
    conditions = session.conditions()
    presence = sources.channel_presence(ions, channels)
    _plot(session, 'quant-reporter-intensities', comparison_plots.reporter_intensities(
        presence, ions, condition_colors(set(conditions.values())), conditions), key, state_manager, height)
    low = [f"{r['channel']} ({r['sample']})" for r in presence.iter_rows(named=True) if r['present'] < 50]
    if low:
        st.warning(f'Found in under half of the labelled spectra: {", ".join(low)}. Check the labelling '
                   'and the channel assignment.')


def build_reporter_ratios(session, state_manager, key, height=None):
    """Ratio of every label to the first label in the same spectrum."""
    got = _reporter_ions(session)
    if got is None:
        return
    ions, channels = got
    if len(channels) < 2:
        st.info('Ratios need at least two labels.')
        return
    reference = channels[0]
    ratios = sources.reporter_ratios(ions, reference)
    names = dict(zip(channels, session.runs))
    _plot(session, 'quant-reporter-ratios', comparison_plots.reporter_ratios(
        ratios, f'{reference} ({names.get(reference, reference)})', [(c, names.get(c, c)) for c in channels[1:]]),
        key, state_manager, height)


def build_reporter_labels(session, state_manager, key, height=None):
    """How many spectra carry 0, 1, ... all reporter ions."""
    got = _reporter_ions(session)
    if got is None:
        return
    ions, channels = got
    _plot(session, 'quant-reporter-labels', comparison_plots.labels_found(
        sources.labels_found(ions, channels), _plex_label(session)), key, state_manager, height)


def _fraction_state(session, key, picker):
    """(dataset id, state) of the labelled run shown by the reporter
    components; the first cell asks which one when there are several."""
    fractions = _fractions(session)
    if fractions is None:
        return None
    names = list(fractions)
    pick_key = f'quant_reporter_fraction_{session.dataset_id}'
    if st.session_state.get(pick_key) not in names:
        st.session_state[pick_key] = names[0]
    if picker and len(names) > 1:
        st.selectbox('Run', names, key=pick_key)
    ds = fractions[st.session_state[pick_key]]
    return ds, _run_state(session, ds, '_reporters')


def build_reporter_spectra(session, state_manager, key, height=None):
    """The reporter ions of every labelled spectrum of a run."""
    del state_manager, height
    from src.render import render
    picked = _fraction_state(session, key, picker=True)
    if picked is not None:
        render._render_cell('reporter_table', session.file_manager, picked[0], 'topdown', picked[1], key)


def build_reporter_ions(session, state_manager, key, height=None):
    """The reporter ions of the spectrum picked in the reporter table."""
    del state_manager, height
    from src.render import render
    picked = _fraction_state(session, key, picker=False)
    if picked is not None:
        render._render_cell('reporter_ions', session.file_manager, picked[0], 'topdown', picked[1], key)


# -- consensus features (label-free) ----------------------------------------

def _features(session):
    if 'feature' not in session.levels:
        raise MissingResult('consensus features')
    long = session.long('feature').collect()
    overview = sources.feature_overview(long, session.annotations('feature').collect(), session.runs)
    return long, overview


def build_feature_map(session, state_manager, key, height):
    """Every consensus feature by retention time and mass, coloured by the
    number of runs that detected it."""
    long, overview = _features(session)
    if 'rt' not in overview.columns or 'mass' not in overview.columns:
        st.info('The features of this comparison have no mass or retention time.')
        return
    _component(
        Heatmap, 'quant_feature_map', session.file_manager, session.dataset_id,
        ['quant_feature', 'quant_feature_annotations'],
        lambda: dict(
            data=overview.filter(pl.col('rt').is_not_null() & pl.col('mass').is_not_null()).lazy(),
            x_column='rt', y_column='mass', intensity_column='Detected', log_scale=False,
            interactivity={'quant_feature': 'RowIndex'}, highlight={'quant_feature': 'RowIndex'},
            title=f'Consensus features (colour: {_sample(session).lower()}s that detected them)',
            x_label='Retention time (min)', y_label='Monoisotopic mass (Da)',
            intensity_label=f'Detected in ({_sample(session).lower()}s)', colorscale='Viridis',
        ),
    )(key=key, state_manager=state_manager, height=height)
    del long


def build_feature_intensity(session, state_manager, key, height=None):
    """Intensity of the features each run detected and of those matched into it."""
    long, overview = _features(session)
    values = sources.feature_values(long, overview, session.runs)
    _plot(session, 'quant-feature-intensity', comparison_plots.feature_intensity(
        values, session.runs, _sample(session)), key, state_manager, height)


def build_feature_table(session, state_manager, key, height):
    """Every consensus feature: mass, retention time, charges and how many
    runs detected it, matched it or miss it."""
    _, overview = _features(session)
    columns = [
        {'title': 'Feature', 'field': 'label', 'sorter': 'string'},
        {'title': 'Mass (Da)', 'field': 'mass', 'sorter': 'number', **_fixed(2)},
        {'title': 'RT (min)', 'field': 'rt', 'sorter': 'number', **_fixed(2)},
        {'title': 'Charges', 'field': 'charges', 'sorter': 'string'},
        {'title': 'Detected', 'field': 'Detected', 'sorter': 'number',
         'headerTooltip': f'{_sample(session)}s that detected it themselves.'},
        {'title': 'Matched', 'field': 'Matched', 'sorter': 'number',
         'headerTooltip': f'{_sample(session)}s it was matched into (match between runs).'},
        {'title': 'Missing', 'field': 'Missing', 'sorter': 'number'},
        {'title': 'Proteoform', 'field': 'proteoform', 'sorter': 'string'},
    ]
    frame = overview
    if 'MinCharge' in frame.columns and 'MaxCharge' in frame.columns:
        frame = frame.with_columns(
            pl.when(pl.col('MinCharge') == pl.col('MaxCharge')).then(pl.col('MinCharge').cast(pl.Utf8))
            .otherwise(pl.col('MinCharge').cast(pl.Utf8) + '–' + pl.col('MaxCharge').cast(pl.Utf8))
            .alias('charges'),
        )
    else:
        frame = frame.with_columns(pl.lit(None, pl.Utf8).alias('charges'))
    if 'label' not in frame.columns:
        frame = frame.with_columns(pl.col('id').alias('label'))
    _component(
        Table, 'quant_feature_table', session.file_manager, session.dataset_id,
        ['quant_feature', 'quant_feature_annotations'],
        lambda: dict(
            data=frame.lazy(), column_definitions=columns, index_field='RowIndex',
            title='Consensus features', go_to_fields=['label'], interactivity={'quant_feature': 'RowIndex'},
            reclick=RECLICK, reset_sort=True,
        ),
    )(key=key, state_manager=state_manager, height=height)


def build_feature_profile(session, state_manager, key, height):
    """The picked consensus feature in every run: detected, matched or missing."""
    long, overview = _features(session)
    values = sources.feature_values(long, overview, session.runs).with_columns(
        (2 ** pl.col('log2_intensity')).alias('intensity'),
    )
    columns = [
        {'title': _sample(session), 'field': 'Run', 'sorter': 'string'},
        {'title': 'Found', 'field': 'Found', 'sorter': 'string'},
        {'title': 'Intensity', 'field': 'intensity', 'sorter': 'number', **_P},
        {'title': 'log2 intensity', 'field': 'log2_intensity', 'sorter': 'number', **_fixed(2)},
    ]
    _component(
        Table, 'quant_feature_profile', session.file_manager, session.dataset_id,
        ['quant_feature', 'quant_feature_annotations'],
        lambda: dict(
            data=values.with_row_index('ValueIndex').with_columns(pl.col('ValueIndex').cast(pl.Int64)).lazy(),
            column_definitions=columns, index_field='ValueIndex',
            title=f'Selected feature per {_sample(session).lower()}', filters={'quant_feature': 'RowIndex'},
            pagination=False, default_row_selects=False,
        ),
    )(key=key, state_manager=state_manager, height=height)


BUILDERS = {
    'quant_summary': (build_summary, None),
    'quant_table': (build_table, 376),
    'quant_volcano': (build_volcano, 480),
    'quant_profile': (build_profile, 480),
    'quant_runs': (build_runs, 300),
    'quant_pca': (build_pca, 480),
    'quant_heatmap': (build_heatmap, 600),
    'quant_missing_values': (build_missing_values, None),
    'quant_imputed_values': (build_imputed_values, None),
    'quant_normalization': (build_normalization, None),
    'quant_pvalues': (build_pvalues, None),
    'quant_go': (build_go, None),
    'quant_id_counts': (build_id_counts, None),
    'quant_id_frequency': (build_id_frequency, None),
    'quant_id_overlap': (build_id_overlap, 480),
    'quant_id_table': (build_id_table, 376),
    'quant_totals': (build_totals, None),
    'quant_correlation': (build_correlation, 480),
    'quant_cv': (build_cv, None),
    'quant_match_quality': (build_match_quality, None),
    'quant_id_prsms': (build_id_prsms, 300),
    'quant_prsm_sequence': (build_prsm_sequence, None),
    'quant_prsm_spectrum': (build_prsm_spectrum, None),
    'quant_reporter_intensities': (build_reporter_intensities, None),
    'quant_reporter_ratios': (build_reporter_ratios, None),
    'quant_reporter_labels': (build_reporter_labels, None),
    'quant_reporter_spectra': (build_reporter_spectra, None),
    'quant_reporter_ions': (build_reporter_ions, None),
    'quant_feature_map': (build_feature_map, 480),
    'quant_feature_intensity': (build_feature_intensity, None),
    'quant_feature_table': (build_feature_table, 376),
    'quant_feature_profile': (build_feature_profile, 300),
}


# Views that size themselves (metrics, GO, a run's own sequence and spectra)
FREE_HEIGHT = {'quant_summary', 'quant_go', 'quant_prsm_sequence', 'quant_prsm_spectrum',
               'quant_reporter_spectra', 'quant_reporter_ions'}
DEFAULT_HEIGHT = 400


def row_height(row, samples=0):
    """One height for the views of a grid row, so their panels line up: the
    tallest view of the row sets it. ``samples``: how many samples the
    comparison has (the normalization view grows with them)."""
    heights = [BUILDERS[name][1] or DEFAULT_HEIGHT for name in row
               if name in BUILDERS and name not in FREE_HEIGHT]
    if 'quant_normalization' in row:
        heights.append(_normalization_height(samples))
    return max(heights) if len(row) > 1 and heights else None


def render_cell(comp_name, session, state_manager, key, height=None):
    builder, default = BUILDERS[comp_name]
    height = height or default
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
        st.caption(f'{_sample(session)}s with the same condition are replicates. Untick **Compared** to leave '
                   f'a {_sample(session).lower()} out; a {_sample(session).lower()} without a condition is left out too.')
        hidden = set(session.hidden)
        table = session.design().to_pandas()
        table.insert(0, 'compared', [s not in hidden for s in table['sample']])
        design = st.data_editor(
            table,
            column_config={
                'compared': st.column_config.CheckboxColumn('Compared', width='small'),
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
            compared = [str(s) for s, keep in zip(design['sample'], design['compared']) if keep]
            if len(compared) < 2:
                st.error(f'A comparison needs at least two {_sample(session).lower()}s.')
                return False
            frame = pl.from_pandas(design.drop(columns='compared')).with_columns(
                pl.col('condition').fill_null('').cast(pl.Utf8).str.strip_chars(),
                pl.col('replicate').fill_null(1).cast(pl.Int64),
            )
            session.save_design(frame)
            session.save_settings({**used, **chosen, 'level': level, 'fc-cutoff': fc, 'padj-cutoff': padj,
                                   'hidden': [s for s in session.all_runs if s not in compared]})
            return True
    return False
