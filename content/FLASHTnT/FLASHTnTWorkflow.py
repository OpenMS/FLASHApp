import pandas as pd
import streamlit as st

from pathlib import Path

from src.parse.tnt import parseTnT
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

    def process_uploaded_files(uploaded_files):
        
        # Store all uploaded files
        for file in uploaded_files:
            if file.name.endswith("mzML"):
                if file.name.endswith('_deconv.mzML'):
                    wf.file_manager.store_file(
                        file.name.split('_deconv.mzML')[0], 'out_deconv_mzML', file
                    )
                elif file.name.endswith('_annotated.mzML'):
                    wf.file_manager.store_file(
                        file.name.split('_annotated.mzML')[0], 'anno_annotated_mzML', file
                    )
                else:
                    st.warning(f'Invalid file : {file.name}')
            elif file.name.endswith("tsv"):
                if file.name.endswith('_tagged.tsv'):
                    wf.file_manager.store_file(
                        file.name.split('_tagged.tsv')[0], 'tags_tsv', file
                    )
                elif file.name.endswith('_protein.tsv'):
                    wf.file_manager.store_file(
                        file.name.split('_protein.tsv')[0], 'protein_tsv', file
                    )
                else:
                    st.warning(f'Invalid file : {file.name}')
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
        for unparsed_dataset in (unparsed_files):
            results = wf.file_manager.get_results(
                unparsed_dataset, 
                ['out_deconv_mzML', 'anno_annotated_mzML', 'tags_tsv', 'protein_tsv']
            )
            
            parsed_data = parseTnT(
                results['out_deconv_mzML'], results['anno_annotated_mzML'], 
                results['tags_tsv'], results['protein_tsv']
            )

            for k, v in parsed_data.items():
                wf.file_manager.store_data(unparsed_dataset, k, v)

    # Upload files via upload widget
    st.subheader("**Upload FLASHDeconv & FLASHTagger output files (\*_annotated.mzML, \*_deconv.mzML, \*_tagged.tsv & \*_protein.tsv)**")
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