import streamlit as st

from src import flashida
from src.Workflow import IdaSimulatorWorkflow
from src.common.common import page_setup

params = page_setup()
if not flashida.available():
    st.info('FLASHIda is a developer feature of the desktop app.')
    st.stop()

wf = IdaSimulatorWorkflow()

st.title('FLASHIda Simulator')
st.caption('🛠️ Developer feature. Simulates which precursors FLASHIda would choose in an mzML.')

t = st.tabs(["📁 **File Upload**", "⚙️ **Configure**", "🚀 **Run**"])
with t[0]:
    wf.show_file_upload_section()
with t[1]:
    wf.show_parameter_section()
with t[2]:
    wf.show_execution_section()
