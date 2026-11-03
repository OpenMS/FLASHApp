import streamlit as st

from src import flashida
from src.common.common import page_setup
from src.render.ida_viewer import show

LAYOUT = [['ida_simulation_heat_map'], ['ida_simulation_table']]

page_setup()
if not flashida.available():
    st.info('FLASHIda is a developer feature of the desktop app.')
    st.stop()

st.title('FLASHIda Simulator Viewer')
st.caption('🛠️ Developer feature.')
show(
    'flashidasimulator', 'simulation_dfs', LAYOUT, 'flashidasimulator',
    'No simulations yet. Run one on the FLASHIda Simulator page.',
)
