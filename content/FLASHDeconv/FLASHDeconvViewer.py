import streamlit as st

from pathlib import Path

from src.common.common import page_setup, save_params, show_linked_grid
from src.workflow.FileManager import FileManager
from src.render.render import make_builders
from src.render.schema import build_insight_caches

# Default panel layout (one experiment): heatmap on top, scan->mass tables,
# annotated + deconvolved spectra, then the precursor-signal 3D plot. Cross-links
# (scan -> mass -> spectrum -> 3D) are carried by the shared StateManager via each
# component's filters/interactivity.
DEFAULT_LAYOUT = [
    ["ms1_deconv_heat_map"],
    ["scan_table", "mass_table"],
    ["anno_spectrum", "deconv_spectrum"],
    ["3D_SN_plot"],
]

# page initialization
params = page_setup()

file_manager = FileManager(
    st.session_state["workspace"],
    Path(st.session_state["workspace"], "cache"),
)

# Gate: need at least one processed FLASHDeconv result.
results = file_manager.get_results_list(["threedim_SN_plot"])
if len(results) == 0:
    st.error("No results to show yet. Please run a workflow first!")
    st.stop()

# A global input sequence enables the Sequence View panel (oracle parity).
has_sequence = file_manager.result_exists("sequence", "sequence")

# Saved layout (trimmed nested list + side_by_side) or the default.
if file_manager.result_exists("layout", "layout"):
    saved = file_manager.get_results("layout", "layout")["layout"]
    layout, side_by_side = saved["layout"], saved["side_by_side"]
else:
    default = DEFAULT_LAYOUT + [["sequence_view"]] if has_sequence else DEFAULT_LAYOUT
    layout, side_by_side = [default], False

# Display-name <-> id mappings for the experiment selectors.
names = [file_manager.get_display_name(r) for r in results]
to_id = {file_manager.get_display_name(r): r for r in results}


def _render_experiment(exp_idx, exp_layout, container):
    """One experiment selector + its linked grid (tool/data-specific, so in-page)."""
    with container:
        sel = st.selectbox("choose experiment", names, key=f"deconv_exp_{exp_idx}")
        ds = to_id[sel]
        # Lazily build the Insight tidy caches for this dataset (idempotent).
        build_insight_caches(file_manager, ds, "flashdeconv")
        builders = make_builders(file_manager, ds, "flashdeconv")
        show_linked_grid([exp_layout], builders, tool=f"flashdeconv_{ds}")


if len(layout) == 2 and side_by_side:
    c1, c2 = st.columns(2)
    _render_experiment(0, layout[0], c1)
    _render_experiment(1, layout[1], c2)
else:
    for i, exp_layout in enumerate(layout):
        if i:
            st.divider()
        _render_experiment(i, exp_layout, st.container())

save_params(params)
