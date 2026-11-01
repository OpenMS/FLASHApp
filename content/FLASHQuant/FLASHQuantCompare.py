"""Compare FLASHQuant results of several runs: link their proteoforms and test
abundance differences between conditions, as OpenDDA compares groups."""

from pathlib import Path

import polars as pl
import streamlit as st

from openms_insight import StateManager

from src.common.common import page_setup, save_params
from src.parse import differential as diff
from src.parse import quant_compare
from src.render.compare import DATASET, RESULTS, render_comparison
from src.workflow.FileManager import FileManager

params = page_setup()

st.title('Compare runs')
st.markdown(
    'Which proteoforms are more or less abundant in one condition than in another? Pick the '
    'runs, give each its condition and compare. The comparison runs with settings that suit '
    'most experiments; change them under **Comparison settings**.'
)
with st.expander('What happens here, in plain words'):
    st.markdown("""
1. **Link runs.** Feature groups of the runs are matched by monoisotopic mass and apex
   retention time with OpenMS's FeatureLinkerUnlabeledKD algorithm, which also corrects
   retention time shifts between the runs. Each linked proteoform gets one abundance per run.
2. **Filter proteoforms.** Proteoforms seen in only a few runs cannot be compared reliably, so
   a proteoform is kept only if it was measured in at least half of the runs of at least one
   condition.
3. **Fill in missing values.** A proteoform missing in a run was usually too faint to detect
   there. It gets that proteoform's lowest measured value.
4. **Normalize.** Abundances are log2-transformed (a difference of 1 means twice as much) and
   each run is shifted so its median matches the others.
5. **Test.** A moderated t-test (similar to *limma*) compares the conditions proteoform by
   proteoform, and the p-values are adjusted with the Benjamini-Hochberg method.

Steps 2 to 5 are those of OpenDDA's *Compare groups*, with the same settings.
""")

runs_cache = FileManager(
    st.session_state['workspace'], Path(st.session_state['workspace'], 'flashquant', 'cache')
)
compare_cache = FileManager(
    st.session_state['workspace'], Path(st.session_state['workspace'], 'flashquant', 'compare-cache')
)
results = runs_cache.get_results_list(['quant_dfs'])

if len(results) < 2:
    st.info('Comparing runs needs FLASHQuant results of at least two runs. Upload them on the File Upload page.')
    st.stop()

previous, saved = {}, {}
if compare_cache.result_exists(DATASET, 'design'):
    design = compare_cache.get_results(DATASET, ['design'], use_polars=True)['design'].collect()
    previous = dict(zip(design['Run'], design['Condition']))
if compare_cache.result_exists(DATASET, 'settings'):
    saved = compare_cache.get_results(DATASET, ['settings'])['settings']
used = {**diff.default_settings(), **saved}

st.subheader('Runs and conditions')
editor = st.data_editor(
    pl.DataFrame({
        'Compare': [r in previous or not previous for r in results],
        'Run': results,
        'Condition': [previous.get(r, '') for r in results],
    }).to_pandas(),
    column_config={
        'Compare': st.column_config.CheckboxColumn('Compare'),
        'Run': st.column_config.TextColumn('Run', disabled=True),
        'Condition': st.column_config.TextColumn(
            'Condition', help='Runs with the same condition are replicates.'
        ),
    },
    hide_index=True, use_container_width=True, key='flashquant_compare_design',
)
chosen = editor[editor['Compare']]
design = {
    run: str(condition).strip()
    for run, condition in zip(chosen['Run'], chosen['Condition'])
    if condition is not None and str(condition).strip()
}
n_conditions = len(set(design.values()))


def _choice(label, options, setting, help=None, column=st):
    value = used[setting]
    return column.selectbox(
        label, options, index=options.index(value) if value in options else 0, help=help,
        key=f'flashquant_compare_{setting}',
    )


with st.expander('Comparison settings'):
    st.markdown('**Linking**')
    mass_column, rt_column = st.columns(2)
    mass_tolerance = mass_column.number_input(
        'Mass tolerance (ppm)', 0.1, 100.0, float(used.get('link-mass-ppm', quant_compare.MASS_TOLERANCE_PPM)), 1.0,
        key='flashquant_compare_mass_tolerance',
    )
    rt_tolerance = rt_column.number_input(
        'Retention time tolerance (s)', 1.0, 1000.0, float(used.get('link-rt-s', quant_compare.RT_TOLERANCE_S)), 5.0,
        key='flashquant_compare_rt_tolerance',
    )
    chosen_settings = {}
    st.markdown('**1. Filter proteoforms**')
    filter_column, threshold_column = st.columns(2)
    method = _choice('Filter by', diff.FILTER_METHODS, 'postproc-filter-method', column=filter_column)
    chosen_settings['postproc-filter-method'] = method
    threshold_setting, threshold_label, threshold_help = {
        'Low Repeatability': (
            'postproc-filter-max-missing', 'Max missing ratio (%)',
            'Keeps proteoforms where at least one condition has at most this share of runs without a value.',
        ),
        'Low Abundance': (
            'postproc-filter-abundance-percentile', 'Threshold percentile (%)',
            "Keeps proteoforms where at least one condition's median is above this percentile.",
        ),
        'Low Variance': (
            'postproc-filter-variance-percentile', 'Threshold percentile (%)',
            "Keeps proteoforms where at least one condition's variance is above this percentile.",
        ),
    }[method]
    chosen_settings[threshold_setting] = threshold_column.slider(
        threshold_label, 0.0, 100.0, float(used[threshold_setting]), 5.0, help=threshold_help,
        key=f'flashquant_compare_{threshold_setting}',
    )
    st.markdown('**2. Missing values**')
    class_column, strategy_column = st.columns(2)
    impute_class = _choice('Why are values missing?', diff.IMPUTE_CLASSES, 'postproc-impute-class', column=class_column)
    chosen_settings['postproc-impute-class'] = impute_class
    if impute_class.startswith('MAR'):
        chosen_settings['postproc-impute-mar-strategy'] = _choice(
            'Fill with the condition', diff.MAR_STRATEGIES, 'postproc-impute-mar-strategy', column=strategy_column,
        )
    else:
        chosen_settings['postproc-impute-mnar-scope'] = _choice(
            'Smallest value of', diff.MNAR_SCOPES, 'postproc-impute-mnar-scope', column=strategy_column,
            help='row: the proteoform\'s lowest value; global: the lowest value of all proteoforms.',
        )
    st.markdown('**3. Normalize**')
    transform_column, norm_column, scaling_column = st.columns(3)
    chosen_settings['postproc-transform'] = _choice('Transformation', diff.TRANSFORMS, 'postproc-transform', column=transform_column)
    chosen_settings['postproc-normalization'] = _choice('Normalization', diff.NORMALIZATIONS, 'postproc-normalization', column=norm_column)
    chosen_settings['postproc-scaling'] = _choice('Scaling', diff.SCALINGS, 'postproc-scaling', column=scaling_column)
    st.markdown('**4. Statistical test**')
    test_column, fdr_column = st.columns(2)
    chosen_settings['postproc-stat-test'] = _choice(
        'Test', diff.test_options(max(n_conditions, 2)), 'postproc-stat-test', column=test_column,
        help='limma_like: moderated t-test (recommended). welch: plain t-test (two conditions). '
        'paired: matched runs (two conditions, same order). anova: any of several conditions differs.',
    )
    chosen_settings['postproc-fdr'] = _choice(
        'Correction', diff.FDR_METHODS, 'postproc-fdr', column=fdr_column,
        help='Benjamini-Hochberg (BH) controls the false discovery rate; Bonferroni is much stricter.',
    )

if st.button('Compare runs', type='primary'):
    try:
        with st.spinner('Linking runs and testing conditions...'):
            runs = {
                run: runs_cache.get_results(run, ['quant_dfs'], use_polars=True)['quant_dfs'].collect()
                for run in chosen['Run']
            }
            consensus, members = quant_compare.link_runs(runs, mass_tolerance, rt_tolerance)
            settings = {**diff.default_settings(), **saved, **chosen_settings}
            statistics, normalized, steps = quant_compare.differential(consensus, design, settings)
            zscores = quant_compare.zscores(statistics, list(design))
        stored = {
            'consensus': consensus, 'members': members, 'statistics': statistics,
            'normalized': normalized, 'filtered': steps['filtered'], 'imputed': steps['imputed'],
            'normalized_wide': steps['normalized'], 'zscores': zscores,
            'design': pl.DataFrame({'Run': list(design), 'Condition': list(design.values())}),
            'settings': {
                **settings, 'postproc-stat-test': steps['stat_test'],
                'link-mass-ppm': mass_tolerance, 'link-rt-s': rt_tolerance,
            },
        }
        for name in RESULTS:
            compare_cache.store_data(DATASET, name, stored[name])
        StateManager(session_key='flashquant_compare_state').clear()
    except ValueError as e:
        st.error(str(e))

if all(compare_cache.result_exists(DATASET, name) for name in RESULTS):
    st.divider()
    render_comparison(compare_cache, StateManager(session_key='flashquant_compare_state'))

save_params(params)
