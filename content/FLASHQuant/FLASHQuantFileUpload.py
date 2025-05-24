import pandas as pd
import streamlit as st

from pathlib import Path

from src.Workflow import QuantWorkflow
from src.parse.quant import parseQuant
from src.common.common import page_setup


# page initialization
params = page_setup()

wf = QuantWorkflow()

st.title("File Upload")

def process_uploaded_files(uploaded_files):
        
        # Store all uploaded files
        for file in uploaded_files:
            
            if file.name.endswith("tsv"):
                if file.name.endswith('.mts.tsv'):
                    wf.file_manager.store_file(
                        file.name.split('.mts.tsv')[0], 'trace_tsv', file
                    )
                elif file.name.endswith('_shared.tsv'):
                    wf.file_manager.store_file(
                        file.name.split('_shared.tsv')[0], 'conflict_tsv', file
                    )
                elif file.name.endswith('.tsv'):
                    wf.file_manager.store_file(
                        file.name.split('.tsv')[0], 'quant_tsv', file
                    )
                else:
                    st.warning(f'Invalid file : {file.name}')
            else:
                st.warning(f'Invalid file : {file.name}')
        
        # Get the unparsed files
        input_files = set(wf.file_manager.get_results_list(['quant_tsv', 'trace_tsv']))
        parsed_files = set(wf.file_manager.get_results_list(['quant_dfs']))
        unparsed_files = input_files - parsed_files

        # Get the unpared conflict resolution files
        tsv_files = set(wf.file_manager.get_results_list(['conflict_tsv']))
        parsed_tsv_files = set(wf.file_manager.get_results_list(['conflict_resolution_dfs']))
        unparsed_tsv_files = (tsv_files - parsed_tsv_files) & input_files

        # Process unparsed datasets
        for unparsed_dataset in (unparsed_files | unparsed_tsv_files):
            results = wf.file_manager.get_results(
                unparsed_dataset, ['quant_tsv', 'trace_tsv']
            )

            conflict_results = None
            if wf.file_manager.result_exists(unparsed_dataset, 'conflict_tsv'):
                conflict_results = wf.file_manager.get_results(
                    unparsed_dataset, ['conflict_tsv']
                )['conflict_tsv']
            
            parsed_data = parseQuant(
                results['quant_tsv'], results['trace_tsv'], conflict_results
            )

            for k, v in parsed_data.items():
                wf.file_manager.store_data(unparsed_dataset, k, v)


# Upload files via upload widget
st.subheader("**Upload FLASHQuant output files**")

# Display info how to upload files
st.info(
    """
**💡 How to upload files**

1. Browse files on your computer or drag and drops files
2. Click the **Add the uploaded quant files** button to use them in the workflows

Select data for analysis from the uploaded files shown below.

**💡 Make sure that the same number of FLASHQuant result files (\*fq.tsv and \*fq.mts.tsv) are uploaded!**

**💡 To visualize conflict resolution, \*fq_shared.tsv files should be uploaded**
"""
)
with st.form('files_uploader_form', clear_on_submit=True):
    uploaded_files = st.file_uploader(
        "FLASHQuant output files", accept_multiple_files=True
    )
    _, c2, _ = st.columns(3)
    # User needs to click button to upload selected files
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
    set(wf.file_manager.get_results_list(['quant_tsv']))
    | set(wf.file_manager.get_results_list(['trace_tsv']))
    | set(wf.file_manager.get_results_list(['conflict_tsv']))
)
table = {
    'Experiment Name' : [],
    'Quant Result Files' : [],
    'Mass Trace Files' : [],
    '(Conflict resolution Files)' : [],
}
for experiment in experiments:
    table['Experiment Name'].append(experiment)

    if wf.file_manager.result_exists(experiment, 'quant_tsv'):
        table['Quant Result Files'].append(True)
    else:
        table['Quant Result Files'].append(False)

    if wf.file_manager.result_exists(experiment, 'trace_tsv'):
        table['Mass Trace Files'].append(True)
    else:
        table['Mass Trace Files'].append(False)

    if wf.file_manager.result_exists(experiment, 'conflict_tsv'):
        table['(Conflict resolution Files)'].append(True)
    else:
        table['(Conflict resolution Files)'].append(False)

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
