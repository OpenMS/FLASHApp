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
    Path(st.session_state['workspace'], 'flashdeconv', 'cache')
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
    columns = st.columns(3)
    columns[0].write('**Run**')
    columns[1].write('**Download**')
    columns[2].write('**Delete Result Set**')

    # Table Body
    for i, experiment in enumerate(experiments):
        st.divider()
        columns = st.columns(3)
        columns[0].empty().write(experiment)
        
        with columns[1]:
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
                            for filepath in file_manager.get_results(
                                experiment, targets, partial=True
                            ).values():
                                f.write(filepath, arcname=Path(filepath).name)
                        zip_buffer.seek(0)
                        file_manager.store_file(
                            experiment, 'download_archive', zip_buffer, 
                            file_name=f'{experiment}.zip'
                        )
                    out_zip = file_manager.get_results(
                        experiment, ['download_archive']
                    )['download_archive']
                    # Show download button after ZIP file was created
                    with open(out_zip, 'rb') as f:
                        button_placeholder.download_button(
                            "Download ⬇️", f, 
                            file_name = f'{experiment}.zip',
                            use_container_width=True
                        )

        with columns[2]:
            if st.button(f"🗑️ {experiment}", use_container_width=True):
                file_manager.remove_results(experiment)
                st.rerun()