import streamlit as st

from src.Workflow import IdaWorkflow
from src.common.common import page_setup


params = page_setup()

wf = IdaWorkflow()

st.title('FLASHIda - Intelligent Data Acquisition')

t = st.tabs(["⚙️ **Configure**", "🚀 **Run**"])


with t[0]:
    wf.show_parameter_section()

with t[1]:
    wf.show_execution_section()