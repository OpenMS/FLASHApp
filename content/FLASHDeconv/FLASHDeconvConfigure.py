import streamlit as st

from src.Workflow import DeconvWorkflow
from src.common.common import page_setup, save_params
from src.common.guide import next_steps, page_intro


params = page_setup()

wf = DeconvWorkflow()

st.title('Describe your experiment')
page_intro(
    """
For most experiments only three things need your attention:

1. **MS data** (below): the runs to deconvolve.
2. **Samples and groups** (first tab): whether to compare the runs, and the condition of each run.
3. **Deconvolution** (second tab): the charge and mass range of your proteoforms, if they differ
   from the defaults.

Everything else starts at values that suit a typical experiment and is saved as you go.
Switch on **Show advanced parameters** only if you need to change a specific setting.
""",
    background="""
- **Deconvolution** turns the m/z peaks of a spectrum, which carry many charges and isotopes
  of each proteoform, into one neutral monoisotopic mass per proteoform.
- A **feature** is one mass followed over the retention time of a run; its **abundance** is
  the area of its raw signal.
- **Comparing runs** links the features of the runs by mass and retention time, fills a run
  where a feature was not detected from its signal there (*match between runs*), and tests
  which masses differ between two conditions.
""",
    background_title='Terms explained',
)

wf.show_parameter_section()

next_steps(
    'flashdeconv',
    ('run', 'Start the deconvolution of the selected runs.'),
    ('sequence', 'Optional: give a proteoform sequence to see its fragment coverage in the Viewer.'),
)

save_params(params)
