import streamlit as st

from src import flashida
from src.common.common import page_setup

page_setup()
if not flashida.available():
    st.info('FLASHIda is a developer feature of the desktop app.')
    st.stop()

from src import flashida_folder as ff  # noqa: E402
from src.Workflow import IdaWorkflow  # noqa: E402
from src.render import ida_runs  # noqa: E402

wf = IdaWorkflow()
folder = str(wf.parameter_manager.get_parameters_from_json().get('method-folder') or '').strip()

st.title('FLASHIda method editor')
st.caption(
    '🛠️ Developer feature. Makes a new method from any method in the method folder, or changes one, '
    'by step of the acquisition. ⚠️ marks what FLASHIda would refuse or what does not fit together; '
    'a greyed-out setting is one FLASHIda does not read with the method as it is.'
)

if not ff.usable(folder):
    st.info('Set the method folder on the **Acquisition** page (Instrument PC) first.')
    st.page_link('content/FLASHIda/FLASHIdaAcquisition.py', label='Acquisition', icon='⚙️')
    st.stop()
ida_runs.method_editor(wf.workflow_dir, folder)
