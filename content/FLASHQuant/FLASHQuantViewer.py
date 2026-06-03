import streamlit as st

from pathlib import Path

from src.common.common import page_setup, save_params, show_linked_grid
from src.workflow.FileManager import FileManager
from src.render.render import make_builders
from src.render.schema import build_insight_caches

# FLASHQuant recipe: a feature Table linked to a Plot3D of that feature's traces
# (Table click sets `feature`; Plot3D filters by `feature`).
DEFAULT_LAYOUT = [["quant_visualization", "quant_traces_3d"]]

# page initialization
params = page_setup()

# FLASHQuant keeps its own workspace-rooted cache (oracle parity).
file_manager = FileManager(
    st.session_state["workspace"],
    Path(st.session_state["workspace"], "flashquant", "cache"),
)

# Gate: need at least one processed FLASHQuant result.
results = file_manager.get_results_list(["quant_dfs"])
if len(results) == 0:
    st.error("No results to show yet. Please run a workflow first!")
    st.stop()

names = [file_manager.get_display_name(r) for r in results]
to_id = {file_manager.get_display_name(r): r for r in results}

sel = st.selectbox("choose experiment", names, key="flashquant_exp_0")
ds = to_id[sel]

# Lazily build the Insight tidy caches for this dataset (idempotent).
build_insight_caches(file_manager, ds, "flashquant")
builders = make_builders(file_manager, ds, "flashquant")
show_linked_grid([DEFAULT_LAYOUT], builders, tool=f"flashquant_{ds}")

save_params(params)
