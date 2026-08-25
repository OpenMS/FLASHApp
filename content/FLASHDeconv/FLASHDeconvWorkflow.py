import time

import pandas as pd
import streamlit as st

from pathlib import Path

from src.Workflow import DeconvWorkflow
from src.parse.deconv import parseDeconv
from src.common.common import page_setup, save_params


params = page_setup()

wf = DeconvWorkflow()

st.title('FLASHDeconv - Ultrafast Deconvolution')

t = st.tabs(["📁 **File Upload**", "⚙️ **Configure**", "🚀 **Run**", "💡 **Manual Result Upload**"])
with t[0]:
    wf.show_file_upload_section()

with t[1]:
    wf.show_parameter_section()

with t[2]:
    wf.show_execution_section()
with t[3]:

    # (file name suffix, name tag, prefix FLASHDeconv uses in its own output)
    UPLOAD_FILE_TYPES = (
        ('deconv.mzML', 'out_deconv_mzML', 'out'),
        ('annotated.mzML', 'anno_annotated_mzML', 'anno'),
        ('spec1.tsv', 'spec1_tsv', ''),
        ('spec2.tsv', 'spec2_tsv', ''),
    )

    def process_uploaded_files(uploaded_files):

        # FLASHDeconv names its output files 'out_deconv.mzML' and
        # 'anno_annotated.mzML', so the part in front of the suffix ('out'
        # and 'anno') is not an experiment name: file those under a common
        # dataset, otherwise the two halves of a run never meet.
        default_dataset = time.strftime('uploaded_%Y%m%d-%H%M%S')

        # Store all uploaded files
        for file in uploaded_files:
            for suffix, name_tag, tool_prefix in UPLOAD_FILE_TYPES:
                if not file.name.endswith(suffix):
                    continue
                experiment = file.name[:-len(suffix)].rstrip('_')
                if experiment in ('', tool_prefix):
                    experiment = default_dataset
                wf.file_manager.store_file(experiment, name_tag, file)
                break
            else:
                st.warning(f'Invalid file : {file.name}')

        # Get the unparsed files
        input_files = set(wf.file_manager.get_results_list(['out_deconv_mzML', 'anno_annotated_mzML']))
        parsed_files = set(wf.file_manager.get_results_list(['deconv_dfs', 'anno_dfs']))
        unparsed_files = input_files - parsed_files

        # Process unparsed datasets
        for unparsed_dataset in unparsed_files:
            results = wf.file_manager.get_results(
                unparsed_dataset, 
                ['out_deconv_mzML', 'anno_annotated_mzML', 
                 'spec1_tsv', 'spec2_tsv'],
                 partial=True
            )
            if not ('out_deconv_mzML' in results and 'anno_annotated_mzML' in results):
                st.warning(
                    f"Experiment '{unparsed_dataset}' needs both the "
                    "deconvolved and the annotated mzML file."
                )
                continue

            with st.spinner(f"Processing '{unparsed_dataset}'..."):
                parseDeconv(
                    wf.file_manager, unparsed_dataset,
                    results['out_deconv_mzML'], results['anno_annotated_mzML'],
                    results.get('spec1_tsv'), results.get('spec2_tsv'),
                    logger=wf.logger
                )

    st.subheader("**Upload FLASHDeconv output files (\*_annotated.mzML & \*_deconv.mzML) or spec1/2 TSV files (Qscore Density Plot only)**")
    st.info(
        """
        **💡 How to upload files**

        1. Browse files on your computer or drag and drops files
        2. Click the **Add the uploaded files** button to use them in the workflows

        Select data for analysis from the uploaded files shown below.

        **💡 Make sure that the same number of deconvolved and annotated mzML files are uploaded!**
        """
    )
    with st.form('input_files', clear_on_submit=True):
        uploaded_files = st.file_uploader(
            "FLASHDeconv output mzML files or TSV files", accept_multiple_files=True, type=["mzML", "tsv"]
        )
        _, c2, _ = st.columns(3)
        if c2.form_submit_button("Add files to workspace", type="primary"):
            if uploaded_files:
                # A list of files is required, since online allows only single upload, create a list
                if type(uploaded_files) != list:
                    uploaded_files = [uploaded_files]

                # opening file dialog and closing without choosing a file results in None upload
                process_uploaded_files(uploaded_files)
                st.success("Successfully added uploaded files!")
            else:
                st.warning("Upload some files before adding them.")

    # File Upload Table
    experiments = (
        set(wf.file_manager.get_results_list(['spec1_tsv']))
        | set(wf.file_manager.get_results_list(['spec2_tsv']))
        | set(wf.file_manager.get_results_list(['out_deconv_mzML']))
        | set(wf.file_manager.get_results_list(['anno_annotated_mzML']))
    )
    table = {
        'Experiment Name' : [],
        'Deconvolved Files' : [],
        'Annotated Files' : [],
        '(MS1 TSV Files)' : [],
        '(MS2 TSV Files)' : [],
    }
    for experiment in experiments:
        table['Experiment Name'].append(experiment)

        if wf.file_manager.result_exists(experiment, 'out_deconv_mzML'):
            table['Deconvolved Files'].append(True)
        else:
            table['Deconvolved Files'].append(False)

        if wf.file_manager.result_exists(experiment, 'anno_annotated_mzML'):
            table['Annotated Files'].append(True)
        else:
            table['Annotated Files'].append(False)

        if wf.file_manager.result_exists(experiment, 'spec1_tsv'):
            table['(MS1 TSV Files)'].append(True)
        else:
            table['(MS1 TSV Files)'].append(False)
        if wf.file_manager.result_exists(experiment, 'spec2_tsv'):
            table['(MS2 TSV Files)'].append(True)
        else:
            table['(MS2 TSV Files)'].append(False)

    st.markdown('**Uploaded experiments in current workspace**')
    st.dataframe(pd.DataFrame(table))

    # Remove files
    with st.expander("🗑️ Remove mzML files"):
        to_remove = st.multiselect(
            "select files", options=experiments
        )
        c1, c2 = st.columns(2)
        if c2.button(
                "Remove **selected**", type="primary", disabled=not any(to_remove)
        ):
            for dataset_id in to_remove:
                wf.file_manager.remove_results(dataset_id)
            st.rerun()

        if c1.button("⚠️ Remove **all**"):
            wf.file_manager.clear_cache()
            st.success("All files removed!")
            st.rerun()

save_params(params)