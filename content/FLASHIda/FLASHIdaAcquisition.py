import streamlit as st

from src import flashida
from src.common.common import page_setup, save_params


params = page_setup()
if not flashida.available():
    st.info('FLASHIda is a developer feature of the desktop app.')
    st.stop()

from src.Workflow import IdaWorkflow  # noqa: E402

wf = IdaWorkflow()

st.title('FLASHIda - Intelligent Data Acquisition')
st.caption('🛠️ Developer feature')

t = st.tabs(["⚙️ **Configure**", "🚀 **Run**"])
with t[0]:
    wf.show_parameter_section()

with t[1]:
    st.info(
        'While running, FLASHIda is started for every new raw file '
        '`<name>FLASHIda_<method>.raw` in the raw file folder, with `<method>.json` '
        'from the method folder. Follow its decisions on the **Viewer** page. '
        'Stop it when the acquisitions are done.'
    )
    wf.show_execution_section()

save_params(params)
