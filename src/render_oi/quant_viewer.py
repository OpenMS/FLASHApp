"""OpenMS-Insight rendering engine for FLASHQuant (migration Phase 3).

FLASHQuant has a single visualization: the feature-group view (mass traces /
XIC / isotope pattern). Unlike Deconv/TnT there is no layout manager — the page
fixes ``[['quant_visualization']]`` — and no cross-component linking; a feature
group selector drives the view's internal selection.

This composes:
  - a feature-group ``Table`` (one row per feature group; click sets
    ``featureGroup``), and
  - a ``FeatureView`` filtered by ``{featureGroup: feature_group}`` over the
    long-format trace points (see ``explode_quant_traces_long``).
"""

from __future__ import annotations

import logging
from typing import Any, Callable, Optional

import polars as pl

from .deconv_viewer import _oi_cache_dir, _load_polars

logger = logging.getLogger(__name__)

FEATURE_GROUP = "featureGroup"

# Feature-group summary columns to show in the selector table (when present).
_FG_SUMMARY_COLUMNS = [
    "FeatureGroupIndex",
    "MonoisotopicMass",
    "AverageMass",
    "StartRetentionTime(FWHM)",
    "EndRetentionTime(FWHM)",
    "HighestApexRetentionTime",
    "FeatureGroupQuantity",
    "MinCharge",
    "MaxCharge",
    "MostAbundantFeatureCharge",
    "IsotopeCosineScore",
]


def build_quant_components(
    dataset_id: str,
    file_manager,
    state_manager,
    key_prefix: str,
) -> Optional[Callable[[], Any]]:
    """Build a render callable for the FLASHQuant feature-group view.

    Returns a zero-arg callable that renders a feature-group selector Table
    above a FeatureView (both sharing ``state_manager``), or None if the quant
    cache is unavailable.
    """
    import streamlit as st
    from openms_insight import FeatureView, Table

    from src.parse.long_format import explode_quant_traces_long

    cache_dir = _oi_cache_dir(file_manager, dataset_id)
    cid = lambda name: f"{dataset_id}__quant_{name}"  # noqa: E731
    skey = lambda name: f"{key_prefix}_{name}"  # noqa: E731

    quant = _load_polars(file_manager, dataset_id, "quant_dfs")
    schema = quant.collect_schema().names()

    # Feature-group selector: one row per group (the wide quant frame already is
    # one row per group), click sets featureGroup.
    summary_cols = [c for c in _FG_SUMMARY_COLUMNS if c in schema]
    # Drop the array columns from the table (keep only scalar summary columns).
    fg_table = Table(
        cache_id=cid("feature_table"),
        data=quant.select(summary_cols) if summary_cols else quant,
        interactivity={FEATURE_GROUP: "FeatureGroupIndex"},
        index_field="FeatureGroupIndex",
        title="Feature Groups",
        cache_path=cache_dir,
    )

    # Long-format trace points, filtered by the selected feature group.
    traces_long = explode_quant_traces_long(quant)
    feature_view = FeatureView(
        cache_id=cid("feature_view"),
        data=traces_long,
        filters={FEATURE_GROUP: "feature_group"},
        charge_column="charge",
        mz_column="mz",
        rt_column="rt",
        intensity_column="intensity",
        isotope_column="isotope",
        # Break the polyline between isotope traces within a charge (matching the
        # legacy per-isotope-trace breaks) instead of one connected line per charge.
        trace_key_column="isotope",
        title="Feature Group Visualization",
        cache_path=cache_dir,
    )

    def _render() -> None:
        fg_table(key=skey("feature_table"), state_manager=state_manager)
        feature_view(key=skey("feature_view"), state_manager=state_manager)

    return _render


def render_experiment_quant(
    dataset_id: str,
    file_manager,
    panel_key: str,
) -> None:
    """Render the FLASHQuant feature-group view with a dedicated StateManager."""
    import streamlit as st
    from openms_insight import StateManager

    state_manager = StateManager(session_key=f"oi_quant_state_{panel_key}")
    try:
        render = build_quant_components(
            dataset_id, file_manager, state_manager, key_prefix=panel_key
        )
        if render is not None:
            render()
        else:
            st.warning("FLASHQuant visualization unavailable")
    except Exception as exc:  # pragma: no cover - defensive UI guard
        logger.exception("Failed to render FLASHQuant view")
        st.error(f"Error rendering FLASHQuant view: {exc}")
