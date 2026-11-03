import time

import pandas as pd
import streamlit as st

from pathlib import Path

from src.parse.tnt import parseTnT
from src.parse.deconv import parseDeconv
from src.Workflow import TopDownWorkflow
from src.common.common import page_setup, save_params, local_folder_files
from src.common.guide import next_steps, page_intro


params = page_setup()

wf = TopDownWorkflow()

st.title('Upload your data')
page_intro(
    """
FLASHDeconv turns the peaks of your spectra into the neutral masses of the
proteoforms behind them; with a protein database, FLASHTnT then identifies
them: it reads short sequence tags from the deconvolved fragment masses and
extends them along the proteins of the database.

- **MS data**: one file per run, mzML, Thermo *.raw* or a Bruker timsTOF *.d*
  folder (vendor files are converted for you).
- **Database** (optional): the protein sequences to search, FASTA or PEFF.
  Decoys are added when the FASTA file has none. Without a database the runs
  are deconvolved only.

No data at hand? Example runs and a database are used when nothing is uploaded.
Already ran the tools yourself? Add their output under **Tool results**.
""",
    background="""
- **mzML** is the open file format for raw MS data. Other vendor formats convert
  to it with [ProteoWizard msConvert](https://proteowizard.sourceforge.io/)
  (keep *peak picking* on).
- **Deconvolution only** suits MS1-only or intact-mass runs, and samples no
  database describes (heavily modified proteins, antibodies, unknown samples):
  the masses are the result.
- A reference proteome as FASTA comes from [UniProt](https://www.uniprot.org/proteomes).
  Only proteins in the database can be identified.
- To **compare runs** (which masses or proteoforms go up or down between
  conditions), upload every run of the experiment; you assign them to
  conditions on the next page.
""",
    background_title='Input files explained',
)

t = st.tabs(["📁 **MS data and database**", "💡 **Tool results**"])
with t[0]:
    wf.show_file_upload_section()
    next_steps(
        'topdown',
        ('configure', 'Choose the runs, whether to identify them, their conditions and the settings.'),
    )
with t[1]:

    # (file name suffix, name tag, prefix the tools use in their own output)
    UPLOAD_FILE_TYPES = (
        ('deconv.mzML', 'out_deconv_mzML', 'out'),
        ('annotated.mzML', 'anno_annotated_mzML', 'anno'),
        # FLASHDeconv's spectrum tables (QScore density plot only)
        ('spec1.tsv', 'spec1_tsv', ''),
        ('spec2.tsv', 'spec2_tsv', ''),
        # FLASHTnT's output, for identified runs
        ('tags.tsv', 'tags_tsv', ''),
        ('tagged.tsv', 'tags_tsv', ''),
        ('protein.tsv', 'protein_tsv', ''),
        # optional FLASHAnnotator output (-out, -out_annotation of the workflow)
        ('annotated.tsv', 'annotated_tsv', ''),
        ('fragments.tsv', 'fragments_tsv', ''),
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
                wf.file_manager.store_file(experiment, name_tag, file, remove=False)
                break
            else:
                st.warning(f'Invalid file : {file.name}')

        # Every dataset with something new to parse: its deconvolution, or
        # identifications added to a deconvolved run
        uploaded = set(wf.file_manager.get_results_list(['out_deconv_mzML', 'anno_annotated_mzML'], partial=True)) \
            | set(wf.file_manager.get_results_list(['tags_tsv', 'protein_tsv'], partial=True))
        for dataset in sorted(uploaded):
            results = wf.file_manager.get_results(
                dataset,
                ['out_deconv_mzML', 'anno_annotated_mzML', 'spec1_tsv', 'spec2_tsv',
                 'tags_tsv', 'protein_tsv', 'annotated_tsv', 'fragments_tsv'],
                partial=True
            )
            if not ('out_deconv_mzML' in results and 'anno_annotated_mzML' in results):
                st.warning(
                    f"Experiment '{dataset}' needs both the deconvolved and the annotated mzML file."
                )
                continue
            identified = 'tags_tsv' in results and 'protein_tsv' in results
            if ('tags_tsv' in results) != ('protein_tsv' in results):
                st.warning(
                    f"Experiment '{dataset}' needs both tags.tsv and protein.tsv to show "
                    "identifications; it is shown deconvolved only."
                )
            deconvolved = wf.file_manager.result_exists(dataset, 'deconv_dfs')
            if deconvolved and (not identified or wf.file_manager.result_exists(dataset, 'protein_dfs')):
                continue

            with st.spinner(f"Processing '{dataset}'..."):
                # The tags are matched against the deconvolved masses, so the
                # mzML files are parsed first
                if not deconvolved:
                    parseDeconv(
                        wf.file_manager, dataset,
                        results['out_deconv_mzML'], results['anno_annotated_mzML'],
                        results.get('spec1_tsv'), results.get('spec2_tsv'),
                        logger=wf.logger
                    )
                if identified:
                    parseTnT(
                        wf.file_manager, dataset,
                        results['out_deconv_mzML'], results['anno_annotated_mzML'],
                        results['tags_tsv'], results['protein_tsv'],
                        annotation_tsv=results.get('fragments_tsv'),
                        annotated_protein_tsv=results.get('annotated_tsv'),
                        logger=wf.logger
                    )

    # Upload files via upload widget
    st.subheader("**Upload the output of the tools: FLASHDeconv's anno_annotated.mzML & out_deconv.mzML (optionally spec1/2.tsv), and for identified runs FLASHTnT's tags.tsv & protein.tsv (optionally FLASHAnnotator's annotated.tsv & fragments.tsv)**")
    # Display info how to upload files
    st.info(
        """
    **💡 How to upload files**
    
    {step1}
    {step2} to use them in the viewer
    
    Select data for analysis from the uploaded files shown below.
    
    **💡 Upload a deconvolved and an annotated mzML file for every run; add tags.tsv and protein.tsv to see its identifications.**
    """.format(
        step1="1. Enter the folder with the files (or pick it with 📁)" if st.session_state.location == "local"
        else "1. Browse files on your computer or drag and drops files",
        step2="2. Click the **Add files from folder** button" if st.session_state.location == "local"
        else "2. Click the **Add files to workspace** button",
    )
    )
    if st.session_state.location == "local":
        # A local install reads the files from a folder instead of uploading them.
        local_files = local_folder_files(
            "tnt-results", "FLASHDeconv & FLASHTnT output files",
            tuple(suffix for suffix, _, _ in UPLOAD_FILE_TYPES),
        )
        if local_files:
            process_uploaded_files(local_files)
            st.success("Successfully added files!")
    else:
        with st.form('input_mzML', clear_on_submit=True):
            uploaded_file = st.file_uploader(
                "FLASHDeconv & FLASHTnT output files", accept_multiple_files=True, type=["mzML", "tsv"]
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