import time

import pandas as pd
import streamlit as st

from pathlib import Path

from src.parse.tnt import parseTnT
from src.parse.deconv import parseDeconv
from src.Workflow import TagWorkflow
from src.common.common import page_setup, save_params


params = page_setup()

wf = TagWorkflow()

st.title('FLASHTnT - Tag and Extend')

t = st.tabs(["📁 **File Upload**", "⚙️ **Configure**", "🚀 **Run**", "💡 **Manual Result Upload**"])
with t[0]:
    wf.show_file_upload_section()

with t[1]:
    wf.show_parameter_section()

with t[2]:
    wf.show_execution_section()
with t[3]:

    # (file name suffix, name tag, prefix the tools use in their own output)
    UPLOAD_FILE_TYPES = (
        ('deconv.mzML', 'out_deconv_mzML', 'out'),
        ('annotated.mzML', 'anno_annotated_mzML', 'anno'),
        ('tags.tsv', 'tags_tsv', ''),
        ('tagged.tsv', 'tags_tsv', ''),
        ('protein.tsv', 'protein_tsv', ''),
    )

    REQUIRED_FILES = (
        ('deconvolved mzML', 'out_deconv_mzML'),
        ('annotated mzML', 'anno_annotated_mzML'),
        ('tags.tsv', 'tags_tsv'),
        ('protein.tsv', 'protein_tsv'),
    )

    def process_uploaded_files(uploaded_files):

        # FLASHDeconv and FLASHTnT name their output files 'out_deconv.mzML',
        # 'anno_annotated.mzML', 'tags.tsv' and 'protein.tsv', so the part in
        # front of the suffix is not an experiment name: file those under a
        # common dataset, otherwise the files of a run never meet.
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
        input_files = set(wf.file_manager.get_results_list(
            ['out_deconv_mzML', 'anno_annotated_mzML', 'tags_tsv', 'protein_tsv']
        ))
        parsed_files = set(wf.file_manager.get_results_list(
            ['deconv_dfs', 'anno_dfs', 'tag_dfs', 'protein_dfs']
        ))
        unparsed_files = input_files - parsed_files

        # Process unparsed datasets
        for unparsed_dataset in unparsed_files:
            results = wf.file_manager.get_results(
                unparsed_dataset, 
                ['out_deconv_mzML', 'anno_annotated_mzML', 'tags_tsv', 'protein_tsv'],
                partial=True
            )
            missing = [n for n, tag in REQUIRED_FILES if tag not in results]
            if missing:
                st.warning(
                    f"Experiment '{unparsed_dataset}' is missing the "
                    f"{', '.join(missing)} file(s)."
                )
                continue

            with st.spinner(f"Processing '{unparsed_dataset}'..."):
                # The tags are matched against the deconvolved masses, so the
                # mzML files have to be parsed first
                if not wf.file_manager.result_exists(
                    unparsed_dataset, 'deconv_tolerance'
                ):
                    parseDeconv(
                        wf.file_manager, unparsed_dataset,
                        results['out_deconv_mzML'], results['anno_annotated_mzML'],
                        logger=wf.logger
                    )
                parseTnT(
                    wf.file_manager, unparsed_dataset,
                    results['out_deconv_mzML'], results['anno_annotated_mzML'],
                    results['tags_tsv'], results['protein_tsv'],
                    logger=wf.logger
                )

    # Upload files via upload widget
    st.subheader("**Upload FLASHDeconv & FLASHTnT output files (anno_annotated.mzML, out_deconv.mzML, tags.tsv & protein.tsv)**")
    # Display info how to upload files
    st.info(
        """
    **💡 How to upload files**
    
    1. Browse files on your computer or drag and drops files
    2. Click the **Add files to workspace** button to use them in the viewer
    
    Select data for analysis from the uploaded files shown below.
    
    **💡 Make sure that the same number of deconvolved and annotated mzML and FLASHTagger output files files are uploaded!**
    """
    )
    with st.form('input_mzML', clear_on_submit=True):
        uploaded_file = st.file_uploader(
            "FLASHDeconv & FLASHTagger output files", accept_multiple_files=True, type=["mzML", "tsv"]
        )
        _, c2, _ = st.columns(3)
        # User needs to click button to upload selected files
        if c2.form_submit_button("Add files to workspace", type="primary"):
            if uploaded_file:
                # A list of files is required, since online allows only single upload, create a list
                if type(uploaded_file) != list:
                    uploaded_file = [uploaded_file]

                # opening file dialog and closing without choosing a file results in None upload
                process_uploaded_files(uploaded_file)
                st.success("Successfully added uploaded files!")
            else:
                st.warning("Upload some files before adding them.")

    # File Upload Table
    experiments = (
        set(wf.file_manager.get_results_list(['tags_tsv']))
        | set(wf.file_manager.get_results_list(['protein_tsv']))
        | set(wf.file_manager.get_results_list(['out_deconv_mzML']))
        | set(wf.file_manager.get_results_list(['anno_annotated_mzML']))
    )
    table = {
        'Experiment Name' : [],
        'Deconvolved Files' : [],
        'Annotated Files' : [],
        'Protein Files' : [],
        'Tag Files' : [],
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

        if wf.file_manager.result_exists(experiment, 'protein_tsv'):
            table['Protein Files'].append(True)
        else:
            table['Protein Files'].append(False)

        if wf.file_manager.result_exists(experiment, 'tags_tsv'):
            table['Tag Files'].append(True)
        else:
            table['Tag Files'].append(False)

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