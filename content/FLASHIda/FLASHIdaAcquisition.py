import time

import pandas as pd
import streamlit as st

from src import flashida
from src.Workflow import IdaWorkflow
from src.parse.deconv import parseDeconv
from src.parse.ida import parseIda
from src.common.common import page_setup, save_params, local_folder_files


params = page_setup()
if not flashida.available():
    st.info('FLASHIda is a developer feature of the desktop app.')
    st.stop()

wf = IdaWorkflow()

st.title('FLASHIda - Intelligent Data Acquisition')
st.caption('🛠️ Developer feature')

t = st.tabs(["⚙️ **Configure**", "🚀 **Run**", "💡 **Load Results**"])
with t[0]:
    wf.show_parameter_section()

with t[1]:
    st.info(
        'While running, FLASHIda is started for every new raw file '
        '`<name>FLASHIda_<method>.raw` in the raw file folder, with `<method>.xml` '
        'from the method folder (`..._2.raw`: the secondary FLASHIda). '
        'Stop it when the acquisitions are done.'
    )
    wf.show_execution_section()

with t[2]:
    # (file name suffix, name tag, prefix FLASHDeconv uses in its own output)
    UPLOAD_FILE_TYPES = (
        ('deconv.mzML', 'out_deconv_mzML', 'out'),
        ('annotated.mzML', 'anno_annotated_mzML', 'anno'),
        ('spec2.tsv', 'spec2_tsv', ''),
        ('prsm_single.tsv', 'id_tsv', ''),
    )

    def experiment_name(file_name, suffix, tool_prefix, default):
        name = file_name[:-len(suffix)].rstrip('_')
        # TopPIC names its output <run>_ms2_toppic_prsm_single.tsv
        for tail in ('_ms2_toppic', '_toppic'):
            if suffix == 'prsm_single.tsv' and name.endswith(tail):
                name = name[:-len(tail)]
        return default if name in ('', tool_prefix) else name

    def process_files(files):
        default_dataset = time.strftime('uploaded_%Y%m%d-%H%M%S')
        for file in files:
            for suffix, name_tag, tool_prefix in UPLOAD_FILE_TYPES:
                if file.name.endswith(suffix):
                    experiment = experiment_name(file.name, suffix, tool_prefix, default_dataset)
                    wf.file_manager.store_file(experiment, name_tag, file, remove=False)
                    break
            else:
                st.warning(f'Invalid file : {file.name}')

        complete = set(wf.file_manager.get_results_list(
            ['out_deconv_mzML', 'anno_annotated_mzML', 'spec2_tsv'], partial=False
        ))
        for dataset in sorted(complete):
            results = wf.file_manager.get_results(
                dataset, ['out_deconv_mzML', 'anno_annotated_mzML', 'spec2_tsv', 'id_tsv'],
                partial=True,
            )
            with st.spinner(f"Processing '{dataset}'..."):
                if not wf.file_manager.result_exists(dataset, 'deconv_dfs'):
                    parseDeconv(
                        wf.file_manager, dataset, results['out_deconv_mzML'],
                        results['anno_annotated_mzML'], None, results['spec2_tsv'],
                        logger=wf.logger,
                    )
                parseIda(wf.file_manager, dataset, results['spec2_tsv'], results.get('id_tsv'))

    st.subheader('Load the results of FLASHIda acquisitions')
    st.info(
        """
        **💡 Files of one acquisition**

        FLASHDeconv's `*_deconv.mzML`, `*_annotated.mzML` and `*_spec2.tsv`, and
        optionally TopPIC's `*_prsm_single.tsv` to mark the identified MS2 scans.
        Files of one acquisition share the part of the name in front of these endings.

        1. Enter the folder with the files (or pick it with 📁)
        2. Click **Add files from folder**
        """
    )
    local_files = local_folder_files(
        'ida-results', 'FLASHIda results', tuple(s for s, _, _ in UPLOAD_FILE_TYPES)
    )
    if local_files:
        process_files(local_files)
        st.success('Added the files.')

    experiments = sorted(
        set(wf.file_manager.get_results_list(['out_deconv_mzML']))
        | set(wf.file_manager.get_results_list(['anno_annotated_mzML']))
        | set(wf.file_manager.get_results_list(['spec2_tsv']))
        | set(wf.file_manager.get_results_list(['id_tsv']))
    )
    table = pd.DataFrame({
        'Experiment': experiments,
        'Deconvolved mzML': [wf.file_manager.result_exists(e, 'out_deconv_mzML') for e in experiments],
        'Annotated mzML': [wf.file_manager.result_exists(e, 'anno_annotated_mzML') for e in experiments],
        'MS2 TSV': [wf.file_manager.result_exists(e, 'spec2_tsv') for e in experiments],
        '(TopPIC PrSMs)': [wf.file_manager.result_exists(e, 'id_tsv') for e in experiments],
        'In Viewer': [wf.file_manager.result_exists(e, 'ida_ms2_table') for e in experiments],
    })
    st.markdown('**Acquisitions in this workspace**')
    st.dataframe(table, hide_index=True)

    with st.expander('🗑️ Remove acquisitions'):
        to_remove = st.multiselect('select acquisitions', options=experiments)
        c1, c2 = st.columns(2)
        if c2.button('Remove **selected**', type='primary', disabled=not to_remove):
            for dataset_id in to_remove:
                wf.file_manager.remove_results(dataset_id)
            st.rerun()
        if c1.button('⚠️ Remove **all**'):
            wf.file_manager.clear_cache()
            st.rerun()

save_params(params)
