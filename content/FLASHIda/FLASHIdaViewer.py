import streamlit as st

from src import flashida
from src.common.common import page_setup
from src.render.ida_viewer import show

LAYOUT = [
    ['ida_ms2_table'],
    ['ida_precursor_spectrum'],
    ['anno_spectrum', 'deconv_spectrum'],
]

page_setup()
if not flashida.available():
    st.info('FLASHIda is a developer feature of the desktop app.')
    st.stop()

st.title('FLASHIda Viewer')
st.caption(
    '🛠️ Developer feature. Click an MS2 scan to see the MS1 spectrum its precursor '
    'was chosen from and the MS2 spectrum.'
)
show(
    'flashida', 'ida_ms2_table', LAYOUT, 'flashida',
    'No acquisitions yet. Load results on the FLASHIda Acquisition page.',
)
