import streamlit as st

from src.Workflow import TagWorkflow
from src.common.common import page_setup, save_params
from src.common.guide import next_steps, page_intro


params = page_setup()

wf = TagWorkflow()

st.title('Describe your experiment')
page_intro(
    """
For most experiments only three things need your attention:

1. **Sample type** (buttons below): native or targeted runs, or lysate LC-MS/MS. A
   recommendation follows the selected runs.
2. **MS data and database** (below): the runs to search and the proteins to search them against.
3. **Samples and groups** (first tab): whether to compare the runs, and the condition of each run.

Everything else starts at values that suit a typical experiment and is saved as you go.
Switch on **Show advanced parameters** only if you need to change a specific setting.
""",
    background="""
- **Deconvolution** turns the m/z peaks of each spectrum into neutral masses: the precursor
  masses and the fragment masses the search reads.
- **Tags** are short amino acid sequences read from the mass differences of fragments. FLASHTnT
  looks them up in the database and extends them to whole proteoforms (**PrSMs**: proteoform-spectrum
  matches), with q-values from decoy proteins.
- **Comparing runs** links the features of the runs, says which proteoform each linked feature is
  (from the PrSMs of the runs), and tests which proteoforms and proteins differ between two conditions.
""",
    background_title='Terms explained',
)

wf.show_parameter_section()

next_steps(
    'flashtnt',
    ('run', 'Start the search of the selected runs.'),
)

save_params(params)
