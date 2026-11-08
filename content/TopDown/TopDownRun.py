import streamlit as st

from src import analysis
from src.Workflow import TopDownWorkflow
from src.common.common import page_setup, save_params
from src.common.guide import PAGES, next_steps, page_intro
from src.workflow.RunManager import RUN_NAME_KEY


params = page_setup()

wf = TopDownWorkflow()

st.title('Run the analysis')
page_intro(
    "Start the analysis and follow its progress. Every selected run is deconvolved and, with "
    "**Identify proteoforms** on, searched in turn; with a **Quantification** chosen on Describe experiment, "
    "the samples are then compared. You can close this page and come back later with the same workspace.\n\n"
    "Each run is kept as an experiment on **Explore results**, and a comparison as an entry "
    "**Comparison: ...** in its experiment selector."
)

wf_params = wf.parameter_manager.get_parameters_from_json()
selected = wf_params.get('mzML-files') or []
selected = [selected] if isinstance(selected, str) else selected
identify = bool(wf_params.get(analysis.IDENTIFY_KEY, True))
# The name given beside Run analysis: before the start in the settings, then kept with it
name = wf_params.get(RUN_NAME_KEY) or (wf.run_name() if wf.get_workflow_status().get('running') else '')
st.info(
    (f"**{name}**: " if name else '') +
    f"**{analysis.name(identify, analysis.quantification(wf_params), analysis.plex(wf_params))}** "
    f"of {len(selected)} run{'s' if len(selected) != 1 else ''}. Change it on **Describe experiment**.",
    icon='🧭',
)

wf.show_execution_section(results_page=PAGES['topdown']['viewer'][0], results_label='Explore results')

next_steps(
    'topdown',
    ('viewer', 'Look at each run, or pick a comparison of runs in the experiment selector.'),
    ('download', 'Take the result files of every run with you: TSV, mzML and msalign, one ZIP per run.'),
)

save_params(params)
