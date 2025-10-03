import streamlit as st

from pathlib import Path
from io import BytesIO
from zipfile import ZipFile, ZIP_DEFLATED

from src.common.common import page_setup
from src.workflow.FileManager import FileManager


page_setup()

st.title('Download')

file_manager = FileManager(
    st.session_state["workspace"],
    Path(st.session_state['workspace'], 'cache')
)

targets = [
    'out_tsv', 'spec1_tsv', 'spec2_tsv', 'spec3_tsv', 'spec4_tsv', 'quant_tsv', 
    'toppic_ms1_msalign', 'toppic_ms1_feature', 'toppic_ms2_msalign', 
    'toppic_ms2_feature', 'out_deconv_mzML', 'anno_annotated_mzML', 
    'FD_parameters_json'
]
experiments = file_manager.get_results_list(targets, partial=True)

# Show error if no content is available for download
if len(experiments) == 0:
    st.error('No results to show yet. Please run a workflow first!')
else:
    # Table Header
    columns = st.columns([1.1, 1, 1])
    columns[0].write('**Name**')
    columns[1].write('**Download**')
    columns[2].write('**Delete Result Set**')

    # Table Body
    for i, experiment in enumerate(experiments):
        st.divider()
        columns = st.columns([0.1, 1, 1, 1])
        current_name = file_manager.get_display_name(experiment)
        
        # Initialize edit mode session state for this experiment
        edit_mode_key = f"edit_mode_{experiment}"
        if edit_mode_key not in st.session_state:
            st.session_state[edit_mode_key] = False
        
        # Display Name or Edit Input
        with columns[1]:
            if st.session_state[edit_mode_key]:
                # Edit mode: Show text input with current display name
                new_name = st.text_input(
                    "New name",
                    value=current_name,
                    key=f"input_{experiment}",
                    label_visibility="collapsed"
                )
            else:
                st.write(current_name)
        
        # Edit/Save Button
        with columns[0]:
            if st.session_state[edit_mode_key]:
                # Show save button in edit mode
                if st.button("💾", key=f"save_{experiment}", help="Save new name", use_container_width=True):
                    new_name = st.session_state.get(f"input_{experiment}", "").strip()
                    
                    # Validate input
                    if not new_name:
                        st.error("Name cannot be empty")
                    elif len(new_name) > 100:
                        st.error("Name is too long (max 100 characters)")
                    else:
                        # Attempt to rename
                        success = file_manager.rename_dataset(experiment, new_name)
                        if success:
                            st.success(f"Renamed to: {new_name}")
                            st.session_state[edit_mode_key] = False
                            st.rerun()
                        else:
                            st.error("Failed to rename dataset")
            else:
                # Show edit button in normal mode
                if st.button("✏️", key=f"edit_{experiment}", help="Edit name", use_container_width=True):
                    st.session_state[edit_mode_key] = True
                    st.rerun()
        
        # Download
        with columns[2]:
            button_placeholder = st.empty()
            
            # Show placeholder button before download is prepared
            clicked = button_placeholder.button('Prepare Download', key=i, use_container_width=True)
            if clicked:
                button_placeholder.empty()
                with st.spinner():
                    # Create ZIP file
                    if not file_manager.result_exists(
                        experiment, 'download_archive'
                    ):
                        zip_buffer = BytesIO()
                        with ZipFile(zip_buffer, 'w', ZIP_DEFLATED) as f:
                            for filepath in file_manager.get_all_files_except(
                                experiment, ['download_archive']
                            ).values():
                                f.write(filepath, arcname=Path(filepath).name)
                        zip_buffer.seek(0)
                        file_manager.store_file(
                            experiment, 'download_archive', zip_buffer, 
                            file_name='download_archive.zip'
                        )
                    out_zip = file_manager.get_results(
                        experiment, ['download_archive']
                    )['download_archive']
                    # Show download button after ZIP file was created
                    with open(out_zip, 'rb') as f:
                        button_placeholder.download_button(
                            "Download ⬇️", f, 
                            file_name = f'{current_name}.zip',
                            use_container_width=True
                        )

        # Delete
        with columns[3]:
            if st.button(f"🗑️ {current_name}", use_container_width=True):
                file_manager.remove_results(experiment)
                st.rerun()