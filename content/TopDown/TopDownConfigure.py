import streamlit as st

from src.Workflow import TopDownWorkflow
from src.common.common import page_setup, save_params
from src.common.guide import PAGES, page_intro


params = page_setup()

wf = TopDownWorkflow()

st.title('Describe your experiment')
page_intro(
    """
Choose your runs and your analysis below. All other settings start at values that suit a typical experiment.

- **MS data**: one file per run, mzML, Thermo *.raw* or a Bruker timsTOF *.d* folder
  (vendor files are converted for you). **📂 Browse** opens your files by folder and adds new ones.
- **Protein database** (to identify proteoforms): FASTA or PEFF. Decoys are added when the
  FASTA file has none. Without a database the runs are deconvolved only.

No data at hand? Example runs and a database are used until you add your own.
""",
    background="""
- **mzML** is the open file format for raw MS data. Other vendor formats convert to it with
  [ProteoWizard msConvert](https://proteowizard.sourceforge.io/) (keep *peak picking* on).
- A reference proteome as FASTA comes from [UniProt](https://www.uniprot.org/proteomes).
  Only proteins in the database can be identified.
- **Deconvolution** turns the m/z peaks of a spectrum, which carry many charges and isotopes
  of each proteoform, into one neutral monoisotopic mass per proteoform: the precursor masses
  and the fragment masses the search reads.
- **Tags** are short amino acid sequences read from the mass differences of fragments. FLASHTnT
  looks them up in the database and extends them to whole proteoforms (**PrSMs**: proteoform-spectrum
  matches), with q-values from decoy proteins.
- A **feature** is one mass followed over the retention time of a run; its **abundance** is
  the area of its raw signal.
- **Label-free quantification** links the features of the runs by mass and retention time, fills
  a run where a feature was not detected from its signal there (*match between runs*) and tests which
  masses differ between two conditions. With identification on, the PrSMs of the runs say which
  proteoform each linked feature is, so proteoforms and proteins are compared too.
- **Isobaric labels** (TMT, iTRAQ) give each sample a tag of the same mass that breaks into a
  reporter ion of its own mass on fragmentation (HCD). The reporter ions of a spectrum give the
  abundance of its precursor in every sample at once.
""",
    background_title='Terms explained',
)

wf.analysis_section()
wf.show_parameter_section()

# Opens Run analysis and starts the analysis there at once
wf.show_run_button(PAGES['topdown']['run'][0], label='Run analysis', summary=wf.run_bar_summary(),
                   help='Save these settings, open Run analysis and start the analysis of the selected runs.')

save_params(params)
