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

# Experiments are selected by their stable dataset id; the display name is shown
# via format_func so duplicate display names can't collapse distinct datasets.


def _render_experiment(exp_idx, exp_layout, container):
    """One experiment selector + its linked grid (tool/data-specific, so in-page)."""
    with container:
        # Oracle parity: start blank (nothing selected) and render nothing until the
        # user picks an experiment -- the old viewer used validate_selected_index
        # (initially None), which also avoided eagerly building caches on page load.
        sel = st.selectbox(
            "choose experiment", results, index=None,
            format_func=file_manager.get_display_name,
            placeholder="Choose an experiment", key=f"deconv_exp_{exp_idx}",
        )
        if sel is None:
            return
        ds = sel
        # Lazily build the Insight tidy caches for this dataset (idempotent).
        build_insight_caches(file_manager, ds, "flashdeconv")
        builders = make_builders(file_manager, ds, "flashdeconv")
        show_linked_grid([exp_layout], builders, tool=f"flashdeconv_{exp_idx}_{ds}")


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
