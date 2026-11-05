import streamlit as st

from pathlib import Path

from src.common.common import page_setup, export_folder_input
from src.common.guide import next_steps, page_intro
from src.common import result_download as rd
from src.workflow.FileManager import FileManager
from src import former_examples


page_setup()

# A local install exports each result set to a folder on disk instead of
# downloading it as a ZIP.
local = rd.is_local()
st.title('Export results' if local else 'Download results')
page_intro(
    ("Copy the result files of your runs to a folder on this computer: one folder per run."
     if local else
     "Download the result files of your runs: one ZIP per run, or every run at once.")
    + " A run holds the deconvolved masses and spectra (TSV, mzML, msalign and feature files), "
    "with **Identify proteoforms** on also the proteins, PrSMs, tags and annotated fragments, "
    "and the parameters it ran with. The run shown in the **Viewer** can also be "
    f"{'exported' if local else 'downloaded'} there with its **{rd.verb()}** button."
)

# A workspace from an earlier release gets today's example runs (once), as in the Viewer
former_examples.add_example_results(Path(st.session_state['workspace'], 'cache'))
file_manager = FileManager(
    st.session_state["workspace"],
    Path(st.session_state['workspace'], 'cache')
)
experiments = rd.downloadable_runs(file_manager)

if len(experiments) == 0:
    st.error('No results yet. Run the analysis first!')
    next_steps('topdown', ('run', 'Start the analysis of your runs.'))
    st.stop()

export_dir = export_folder_input('results-export') if local else None
rd.all_runs_button(file_manager, experiments, export_dir)

for i, experiment in enumerate(experiments):
    current_name = file_manager.get_display_name(experiment)
    files = rd.result_files(file_manager, experiment)
    edit_mode_key = f"edit_mode_{experiment}"
    if edit_mode_key not in st.session_state:
        st.session_state[edit_mode_key] = False

    with st.container(border=True):
        name_col, download_col, delete_col = st.columns([3, 1.2, 0.5], vertical_alignment='center')
        with name_col:
            if st.session_state[edit_mode_key]:
                with st.container(horizontal=True, vertical_alignment='bottom'):
                    st.text_input("New name", value=current_name, key=f"input_{experiment}")
                    if st.button("Save", key=f"save_{experiment}", icon=':material/save:'):
                        new_name = st.session_state.get(f"input_{experiment}", "").strip()
                        if not new_name:
                            st.error("Name cannot be empty")
                        elif len(new_name) > 100:
                            st.error("Name is too long (max 100 characters)")
                        elif file_manager.rename_dataset(experiment, new_name):
                            st.session_state[edit_mode_key] = False
                            st.rerun()
                        else:
                            st.error("Failed to rename dataset")
            else:
                with st.container(horizontal=True, vertical_alignment='center'):
                    st.markdown(f"**{current_name}**")
                    if st.button("", key=f"edit_{experiment}", icon=':material/edit:',
                                 help="Rename", type='tertiary'):
                        st.session_state[edit_mode_key] = True
                        st.rerun()
            st.caption(f"{len(files)} files, {rd.format_size(rd.total_size(files))}")
        message_area = st.container()
        with download_col:
            rd.run_button(file_manager, experiment, i, export_dir,
                          label=None if local else 'Download ZIP', message_area=message_area)
        with delete_col:
            if st.button("", key=f"delete_{experiment}", icon=':material/delete:',
                         help=f"Delete {current_name} and its results", width='stretch'):
                file_manager.remove_results(experiment)
                st.rerun()
        with st.expander('Files'):
            rd.file_table(files)
