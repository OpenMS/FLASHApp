import streamlit as st

from src.Workflow import TopDownWorkflow
from src.common.common import page_setup, save_params
from src.common.guide import next_steps, page_intro


params = page_setup()

wf = TopDownWorkflow()

st.title('Describe your experiment')
page_intro(
    """
For most experiments only a few things need your attention:

1. **MS data** (below): the runs to analyse.
2. **Identify proteoforms** (below): on to search the runs against a protein database, off to
   deconvolve them only (MS1-only or intact-mass runs, or no database for the sample).
3. **Sample type** (buttons below, when identifying): native or targeted runs, or lysate LC-MS/MS.
   A recommendation follows the selected runs.
4. **Samples and groups** (first tab): whether to compare the runs, and the condition of each run.

Everything else starts at values that suit a typical experiment and is saved as you go.
Switch on **Show advanced parameters** only if you need to change a specific setting.
""",
    background="""
- **Deconvolution** turns the m/z peaks of a spectrum, which carry many charges and isotopes
  of each proteoform, into one neutral monoisotopic mass per proteoform: the precursor masses
  and the fragment masses the search reads.
- **Tags** are short amino acid sequences read from the mass differences of fragments. FLASHTnT
  looks them up in the database and extends them to whole proteoforms (**PrSMs**: proteoform-spectrum
  matches), with q-values from decoy proteins.
- A **feature** is one mass followed over the retention time of a run; its **abundance** is
  the area of its raw signal.
- **Comparing runs** links the features of the runs by mass and retention time, fills a run
  where a feature was not detected from its signal there (*match between runs*) and tests which
  masses differ between two conditions. With identification on, the PrSMs of the runs say which
  proteoform each linked feature is, so proteoforms and proteins are compared too.
""",
    background_title='Terms explained',
)

wf.show_parameter_section()

next_steps(
    'topdown',
    ('run', 'Start the analysis of the selected runs.'),
    ('sequence', 'Optional: give a proteoform sequence to see its fragment coverage in the Viewer.'),
)

save_params(params)
