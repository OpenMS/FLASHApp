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

from .deconv_viewer import _column_definitions, _load_polars, _oi_cache_dir

logger = logging.getLogger(__name__)

FEATURE_GROUP = "featureGroup"

# Feature-group selector columns mirroring the legacy
# ``FLASHQuantView.featureGroupTableColumnDefinitions`` as (field, title, tooltip,
# is_float). The explicit titles match the legacy human-readable labels -- without
# column_definitions the OI Table auto-titles the CamelCase fields ("Monoisotopicmass",
# "Startretentiontime(Fwhm)", ...). Legacy did not surface HighestApexRetentionTime,
# so it is intentionally omitted.
_FG_TABLE_COLUMNS = [
    (
        "FeatureGroupIndex",
        "Index",
        "The sequential index of the feature group in the dataset.",
        False,
    ),
    (
        "MonoisotopicMass",
        "Monoisotopic Mass",
        "The monoisotopic mass of the feature group in Daltons.",
        True,
    ),
    (
        "AverageMass",
        "Average Mass",
        "The average mass of the feature group in Daltons.",
        True,
    ),
    (
        "StartRetentionTime(FWHM)",
        "Start Retention Time (FWHM)",
        "The start of the feature group's elution window (full width at half "
        "maximum) in seconds.",
        True,
    ),
    (
        "EndRetentionTime(FWHM)",
        "End Retention Time (FWHM)",
        "The end of the feature group's elution window (full width at half "
        "maximum) in seconds.",
        True,
    ),
    (
        "FeatureGroupQuantity",
        "Feature Group Quantity",
        "The integrated abundance (quantity) of the feature group.",
        True,
    ),
    (
        "MinCharge",
        "Min Charge",
        "The minimum charge state observed for the feature group.",
        False,
    ),
    (
        "MaxCharge",
        "Max Charge",
        "The maximum charge state observed for the feature group.",
        False,
    ),
    (
        "MostAbundantFeatureCharge",
        "Most Abundant Charge",
        "The charge state of the most abundant feature in the group.",
        False,
    ),
    (
        "IsotopeCosineScore",
        "Isotope Cosine Score",
        "The cosine similarity between the observed and theoretical isotope "
        "patterns.",
        True,
    ),
]

# Field projection for the selector table (column order follows the spec above).
_FG_SUMMARY_COLUMNS = [field for field, *_ in _FG_TABLE_COLUMNS]


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
    table_data = quant.select(summary_cols) if summary_cols else quant
    fg_table = Table(
        cache_id=cid("feature_table"),
        data=table_data,
        interactivity={FEATURE_GROUP: "FeatureGroupIndex"},
        index_field="FeatureGroupIndex",
        # Legacy column titles/tooltips + number sorter/precision (without these the
        # OI Table auto-titles the CamelCase fields, e.g. "Monoisotopicmass").
        column_definitions=_column_definitions(summary_cols, _FG_TABLE_COLUMNS),
        title="Feature groups",
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
        # Match the legacy 3D plot title ("Feature group signals").
        title="Feature group signals",
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

    # Scope the StateManager to the dataset so switching the selected experiment
    # starts from a clean selection (OI Table default-row-0) instead of inheriting
    # the previous dataset's featureGroup -- the legacy render_grid reset selections
    # on dataset change (src/render/render.py:80-82).
    state_manager = StateManager(
        session_key=f"oi_quant_state_{panel_key}_{dataset_id}"
    )
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
