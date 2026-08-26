import streamlit as st

from pathlib import Path

from src.workflow.FileManager import FileManager
from src.common.common import page_setup, save_params

# NOTE (Stage D rewiring): FLASHQuant now renders through the reusable
# `openms_insight.FeatureView` component instead of the bespoke
# `flash_viewer_grid` / `FLASHQuantView` path. The old render path is left
# importable on purpose (do NOT delete) so it can be restored or compared.
#   from src.render.components import flash_viewer_grid_component, FlashViewerComponent, FLASHQuant
#   from src.render.render import render_grid
from openms_insight import FeatureView

# page initialization
params = page_setup()


# Get available results
workspace = st.session_state["workspace"]
file_manager = FileManager(
    workspace,
    Path(workspace, 'flashquant', 'cache')
)
results = file_manager.get_results_list(
    ['quant_dfs']
)

### if no input file is given, show blank page
if len(results) == 0:
    st.error('No results to show yet. Please run a workflow first!')
    st.stop()

# Map names to index
name_to_index = {n: i for i, n in enumerate(results)}


# FLASHQuant is a single-experiment, single-component page (no cross-linking,
# no configurable grid). Pick one experiment and render one FeatureView for it.
st.selectbox("choose experiment", results, key="selected_experiment0_quant")
selected_exp0 = st.session_state.selected_experiment0_quant

# Load the parsed feature-group frame produced by src/parse/flashquant.py
# (`connectTraceWithResult`): the 12 scalar columns plus the per-feature-group
# array columns Charges / IsotopeIndices / CentroidMzs / RTs / MZs / Intensities,
# where RTs/MZs/Intensities elements are comma-joined point strings. FeatureView
# consumes this frame directly (no transformation needed).
quant_df = file_manager.get_results(selected_exp0, ['quant_dfs'])['quant_dfs']

# Cache id is per-experiment so switching experiments yields an independent,
# correctly-scoped cache and selection. The cache lives under the FLASHQuant
# workspace cache directory.
feature_view = FeatureView(
    cache_id=f'flashquant_{selected_exp0}',
    data=quant_df,
    cache_path=str(Path(workspace, 'flashquant', 'cache', 'featureview')),
)
feature_view(key=f'flashquant_featureview_{selected_exp0}')

save_params(params)
