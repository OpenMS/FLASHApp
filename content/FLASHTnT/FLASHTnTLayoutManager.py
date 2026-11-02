from pathlib import Path

import streamlit as st

from src.common.common import page_setup, save_params
from src.workflow.FileManager import FileManager
from src.render.layout_editor import layout_editor
from src.render.layout_editor.catalog import TNT_COMPONENTS, TNT_DEFAULT_LAYOUT, TNT_DEFAULT_LAYOUTS

# page_setup() first: it initialises st.session_state['workspace'], which the
# cache access below needs when this page is opened by direct link.
params = page_setup()

# Setup cache access
file_manager = FileManager(
    st.session_state["workspace"],
    Path(st.session_state['workspace'], 'cache')
)

st.title("Layout Manager")
st.caption(
    "Arrange the components shown on the **👀 Viewer** page. Drag components "
    "into the grid (up to 3 per row), add experiments to compare several files, "
    "then press **Save**."
)

EDITOR_KEY = 'tnt_layout_editor'


def start_from(sample_type):
    """Put a default layout into the editor as an unsaved draft."""
    st.session_state['tnt_layout_start'] = sample_type
    st.session_state[f'{EDITOR_KEY}_revision'] = st.session_state.get(f'{EDITOR_KEY}_revision', 0) + 1


st.markdown(
    "Without a saved layout, the Viewer shows the **Native** default for native or "
    "targeted runs (few spectra) and the **Lysate** default for LC-MS/MS runs of a "
    "lysate, chosen per experiment. Start from either:"
)
columns = st.columns(len(TNT_DEFAULT_LAYOUTS) + 2)
for column, sample_type in zip(columns, TNT_DEFAULT_LAYOUTS):
    column.button(f'{sample_type} default', key=f'tnt_layout_start_{sample_type}',
                  on_click=start_from, args=(sample_type,), use_container_width=True)

# A save or reset ends the started draft
submitted = st.session_state.get(EDITOR_KEY)
if submitted and submitted.get('nonce') != st.session_state.get(f'{EDITOR_KEY}_handled'):
    st.session_state.pop('tnt_layout_start', None)
start = st.session_state.get('tnt_layout_start')

layout_editor(file_manager, 'flashtnt_layout', TNT_COMPONENTS, TNT_DEFAULT_LAYOUT, key=EDITOR_KEY,
              start_layout=[TNT_DEFAULT_LAYOUTS[start]] if start else None)

save_params(params)
