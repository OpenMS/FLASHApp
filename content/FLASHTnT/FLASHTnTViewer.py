import streamlit as st

from pathlib import Path

from src.common.common import page_setup, save_params, show_linked_grid
from src.workflow.FileManager import FileManager
from src.render.render import make_builders
from src.render.schema import build_insight_caches

# Default panel layout (one experiment): protein table -> sequence view ->
# tag table -> augmented spectrum. Cross-links (protein -> tag -> sequence;
# tag/peak -> mass) are carried by the shared StateManager.
DEFAULT_LAYOUT = [
    ["protein_table"],
    ["sequence_view"],
    ["tag_table"],
    ["combined_spectrum"],
]

# page initialization
params = page_setup("TaggerViewer")

file_manager = FileManager(
    st.session_state["workspace"],
    Path(st.session_state["workspace"], "cache"),
)

# Gate: need at least one processed FLASHTnT result.
results = file_manager.get_results_list(["protein_dfs"])
if len(results) == 0:
    st.error("No results to show yet. Please run a workflow first!")
    st.stop()

# Saved layout (trimmed nested list + side_by_side) or the default.
if file_manager.result_exists("flashtnt_layout", "layout"):
    saved = file_manager.get_results("flashtnt_layout", "layout")["layout"]
    layout, side_by_side = saved["layout"], saved["side_by_side"]
else:
    layout, side_by_side = [DEFAULT_LAYOUT], False

# Display-name <-> id mappings for the experiment selectors.
names = [file_manager.get_display_name(r) for r in results]
to_id = {file_manager.get_display_name(r): r for r in results}


def _render_experiment(exp_idx, exp_layout, container):
    """One experiment selector + its linked grid (tool/data-specific, so in-page)."""
    with container:
        # Oracle parity: blank until the user picks (no eager cache build on load).
        sel = st.selectbox(
            "choose experiment", names, index=None,
            placeholder="Choose an experiment", key=f"tnt_exp_{exp_idx}",
        )
        if sel is None:
            return
        ds = to_id[sel]
        # Lazily build the Insight tidy caches for this dataset (idempotent).
        build_insight_caches(file_manager, ds, "flashtnt")
        # round-8 finding 3-tables-002: per-experiment "Best per spectrum" toggle
        # (oracle ProteinTable ``bestPerSpectrumOnly``, default ON). Placed ABOVE the
        # grid so it governs this experiment's protein table; its value selects the
        # filtered vs full protein-table row set (+ cache_id) in make_builders.
        best_per_spectrum = st.checkbox(
            "Best per spectrum", value=True, key=f"tnt_best_{exp_idx}",
            help="Show only the highest-scoring proteoform per spectrum (scan). "
                 "Uncheck to show all proteoforms.",
        )
        # SequenceView ion-types / tolerance come from the oracle settings cache.
        settings = None
        if file_manager.result_exists(ds, "settings"):
            settings = file_manager.get_results(ds, ["settings"])["settings"]
        builders = make_builders(
            file_manager, ds, "flashtnt", settings=settings,
            best_per_spectrum=best_per_spectrum,
        )
        show_linked_grid([exp_layout], builders, tool=f"flashtnt_{ds}")


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
