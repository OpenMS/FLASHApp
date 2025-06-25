import streamlit as st

from src.Workflow import IdaSimulatorWorkflow
from src.common.common import page_setup


params = page_setup()

wf = IdaSimulatorWorkflow()

st.title('FLASHIda - Intelligent Data Acquisition')

t = st.tabs(["📁 **File Upload**", "⚙️ **Configure**", "🚀 **Run**"])
with t[0]:
    wf.show_file_upload_section()

with t[1]:
    wf.show_parameter_section()

with t[2]:
    wf.show_execution_section()