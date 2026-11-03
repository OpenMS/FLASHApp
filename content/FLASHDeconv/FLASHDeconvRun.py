import streamlit as st

from src.Workflow import DeconvWorkflow
from src.common.common import page_setup, save_params
from src.common.guide import next_steps, page_intro


params = page_setup()

wf = DeconvWorkflow()

st.title('Run the analysis')
page_intro(
    "Start the analysis and follow its progress. Every selected run is processed in turn; "
    "with **Compare the runs** on (Describe experiment), the runs are then linked and compared. "
    "You can close this page and come back later with the same workspace.\n\n"
    "Each run is kept as an experiment of the **Viewer**, and a comparison as an entry "
    "**Comparison: ...** in its experiment selector."
)

wf.show_execution_section()

next_steps(
    'flashdeconv',
    ('viewer', 'Look at each run, or pick a comparison of runs in the experiment selector.'),
)

save_params(params)
