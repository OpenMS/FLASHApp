"""Compare FLASHQuant results of several runs: link their proteoforms and test
abundance differences between conditions."""

from pathlib import Path

import polars as pl
import streamlit as st

from openms_insight import StateManager

from src.common.common import page_setup, save_params
from src.parse import quant_compare
from src.render.compare import DATASET, RESULTS, render_comparison
from src.workflow.FileManager import FileManager

params = page_setup()

st.title('Compare runs')

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

previous = {}
if compare_cache.result_exists(DATASET, 'design'):
    design = compare_cache.get_results(DATASET, ['design'], use_polars=True)['design'].collect()
    previous = dict(zip(design['Run'], design['Condition']))

with st.form('compare_runs'):
    st.markdown(
        'Pick the runs and give each one its condition. Proteoforms are linked across runs '
        'by monoisotopic mass and apex retention time (OpenMS FeatureLinkerUnlabeledKD, which '
        'also corrects retention time shifts), then compared between conditions with a '
        'moderated t-test (two conditions) or F-test (more), Benjamini-Hochberg adjusted.'
    )
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
    mass_column, rt_column, present_column = st.columns(3)
    mass_tolerance = mass_column.number_input(
        'Mass tolerance (ppm)', 0.1, 100.0, quant_compare.MASS_TOLERANCE_PPM, 1.0,
    )
    rt_tolerance = rt_column.number_input(
        'Retention time tolerance (s)', 1.0, 1000.0, quant_compare.RT_TOLERANCE_S, 5.0,
    )
    min_present = present_column.number_input(
        'Quantified in at least (share of a condition\'s runs)', 0.0, 1.0,
        quant_compare.MIN_PRESENT_RATIO, 0.05,
    )
    submitted = st.form_submit_button('Compare runs', type='primary')

if submitted:
    chosen = editor[editor['Compare']]
    design = {
        run: str(condition).strip()
        for run, condition in zip(chosen['Run'], chosen['Condition'])
        if condition is not None and str(condition).strip()
    }
    try:
        with st.spinner('Linking runs and testing conditions...'):
            runs = {
                run: runs_cache.get_results(run, ['quant_dfs'], use_polars=True)['quant_dfs'].collect()
                for run in chosen['Run']
            }
            consensus, members = quant_compare.link_runs(runs, mass_tolerance, rt_tolerance)
            statistics, normalized = quant_compare.differential(consensus, design, min_present)
            zscores = quant_compare.zscores(statistics, list(design))
        for name, data in zip(RESULTS, [
            consensus, members, statistics, normalized, zscores,
            pl.DataFrame({'Run': list(design), 'Condition': list(design.values())}),
        ]):
            compare_cache.store_data(DATASET, name, data)
        StateManager(session_key='flashquant_compare_state').clear()
    except ValueError as e:
        st.error(str(e))

if compare_cache.result_exists(DATASET, 'zscores'):
    render_comparison(compare_cache, StateManager(session_key='flashquant_compare_state'))

save_params(params)
