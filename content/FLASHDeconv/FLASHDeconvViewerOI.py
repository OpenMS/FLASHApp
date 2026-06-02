"""FLASHDeconv viewer rendered entirely with OpenMS-Insight components (Stage B).

This is the NEW viewer for the FLASHApp -> OpenMS-Insight visualization migration.
It renders the FLASHDeconv workflow using the reusable ``openms_insight`` component
library (``Table``, ``LinePlot``, ``Heatmap``, ``Scatter3D``, ``DensityPlot``,
``SequenceView``) instead of the bespoke ``flash_viewer_grid`` Vue grid in
``src/render/*``.

Design goals (see ``/home/user/parity/STRATEGY.md`` §4/§5):

* ONE shared ``StateManager`` per rendered experiment panel, keyed by a DISTINCT
  ``session_key`` (``svc_state_deconv_<experiment_id>``) so that selections never
  leak between side-by-side experiment panels (HARD edge #6).
* Layout parity: the ``[experiment][row][col]`` nested grid is reproduced with
  ``st.columns`` per row (<=3 cols), rows stacked; multi-experiment side-by-side
  uses a top-level ``st.columns`` (<=5 panels).
* The component->frame->filters/interactivity wiring exactly mirrors the schema
  from the long-format parse producers in ``src/parse/deconv.py``.

The OLD render path (``src/render/render.py`` / ``flash_viewer_grid``) is left
intact and importable; the page chooses which path to use.

NOTE ON CACHES: every OpenMS-Insight component persists a preprocessed cache under
``{cache_path}/{cache_id}/``. We derive a per-experiment cache directory inside the
workspace so the caches live next to the FLASHApp parquet cache and are stable
across reruns. ``cache_id`` is suffixed with the experiment id to keep experiments
isolated on disk as well as in session state.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import polars as pl
import streamlit as st

from openms_insight import (
    DensityPlot,
    Heatmap,
    LinePlot,
    Scatter3D,
    SequenceView,
    StateManager,
    Table,
)

# Map the layout COMPONENT_NAMES (FLASHDeconvLayoutManager) to a builder. Every
# builder returns a *callable* OpenMS-Insight component already wired with the
# shared filters/interactivity identifiers. The identifiers below are the FLASHApp
# StateTracker keys (scanIndex / massIndex / heatmap zoom ids) so that state flows
# across components exactly like the legacy grid.

SCAN_KEY = "scanIndex"
MASS_KEY = "massIndex"

# Curated column definitions mirroring the LEGACY Vue tables (titles / order /
# field selection). The OI Table's ``_get_columns_to_select`` projects to ONLY the
# fields named here (plus index / interactivity / filter columns), so any internal
# frame column not listed is hidden -- the visual-parity goal.

# TabulatorScanTable.vue columns -> scan_table fields. Legacy "Index" (id) maps to
# the frame's `index` (row position == scan index).
_SCAN_COLUMN_DEFINITIONS = [
    {"title": "Index", "field": "index", "sorter": "number"},
    {"title": "Scan Number", "field": "Scan", "sorter": "number"},
    {"title": "MS Level", "field": "MSLevel", "sorter": "number"},
    {"title": "Retention time", "field": "RT", "sorter": "number",
     "formatter": "money", "formatterParams": {"precision": 2, "symbol": ""}},
    {"title": "Precursor Mass", "field": "PrecursorMass", "sorter": "number",
     "formatter": "money", "formatterParams": {"precision": 2, "symbol": ""}},
    {"title": "#Masses", "field": "#Masses", "sorter": "number"},
]

# TabulatorMassTable.vue columns -> mass_table_long fields. Legacy "Index" (id) maps
# to the long frame's `mass_id` (0-based mass position within the scan).
_MASS_COLUMN_DEFINITIONS = [
    {"title": "Index", "field": "mass_id", "sorter": "number"},
    {"title": "Monoisotopic mass", "field": "MonoMass", "sorter": "number",
     "formatter": "money", "formatterParams": {"precision": 2, "symbol": ""}},
    {"title": "Sum intensity", "field": "SumIntensity", "sorter": "number",
     "formatter": "money", "formatterParams": {"precision": 2, "symbol": ""}},
    {"title": "Min charge", "field": "MinCharges", "sorter": "number"},
    {"title": "Max charge", "field": "MaxCharges", "sorter": "number"},
    {"title": "Min isotope", "field": "MinIsotopes", "sorter": "number"},
    {"title": "Max isotope", "field": "MaxIsotopes", "sorter": "number"},
    {"title": "Cosine score", "field": "CosineScore", "sorter": "number",
     "formatter": "money", "formatterParams": {"precision": 2, "symbol": ""}},
    {"title": "SNR", "field": "SNR", "sorter": "number",
     "formatter": "money", "formatterParams": {"precision": 2, "symbol": ""}},
    {"title": "QScore", "field": "QScore", "sorter": "number",
     "formatter": "money", "formatterParams": {"precision": 2, "symbol": ""}},
]


def _component_cache_dir(file_manager, experiment_id: str) -> str:
    """Directory under the workspace cache where OI component caches are written."""
    cache_root = Path(file_manager.cache_path, "oi_components", str(experiment_id))
    cache_root.mkdir(parents=True, exist_ok=True)
    return str(cache_root)


def _data_path(file_manager, experiment_id: str, name_tag: str) -> Optional[str]:
    """Resolve the on-disk parquet path for a stored frame, or None if absent."""
    if not file_manager.result_exists(experiment_id, name_tag):
        return None
    res = file_manager.get_results(experiment_id, [name_tag], partial=True)
    path = res.get(name_tag)
    return str(path) if path is not None else None


def _lazy(file_manager, experiment_id: str, name_tag: str) -> Optional[pl.LazyFrame]:
    """Load a stored frame as a polars LazyFrame, or None if absent."""
    if not file_manager.result_exists(experiment_id, name_tag):
        return None
    return file_manager.get_results(
        experiment_id, [name_tag], use_polars=True
    )[name_tag]


# ---------------------------------------------------------------------------
# Per-component builders. Each returns an OpenMS-Insight component instance, or
# None when the underlying data frame is missing (component is silently skipped).
# ---------------------------------------------------------------------------

def _build_heatmap(
    file_manager, experiment_id: str, cache_dir: str, frame_tag: str,
    zoom_id: str, title: str,
):
    data = _lazy(file_manager, experiment_id, frame_tag)
    if data is None:
        return None
    # Long heatmap frames carry columns: mass, rt, intensity, scan_idx, mass_idx.
    # Axes per Heatmap.md: x = Retention Time (rt), y = Mass (mass).
    return Heatmap(
        cache_id=f"{frame_tag}_{experiment_id}",
        data=data,
        x_column="rt",
        y_column="mass",
        intensity_column="intensity",
        zoom_identifier=zoom_id,
        title=title,
        x_label="Retention Time",
        y_label="Mass",
        cache_path=cache_dir,
    )


def _build_scan_table(file_manager, experiment_id: str, cache_dir: str):
    data = _lazy(file_manager, experiment_id, "scan_table")
    if data is None:
        return None
    # Scan table: clicking a row sets scanIndex to the row's `index`.
    return Table(
        cache_id=f"scan_table_{experiment_id}",
        data=data,
        interactivity={SCAN_KEY: "index"},
        index_field="index",
        column_definitions=_SCAN_COLUMN_DEFINITIONS,
        go_to_fields=["index", "Scan"],
        title="Scan Table",
        cache_path=cache_dir,
    )


def _build_mass_table(file_manager, experiment_id: str, cache_dir: str):
    data = _lazy(file_manager, experiment_id, "mass_table_long")
    if data is None:
        return None
    # Mass table (long): filtered to the selected scan via `index`; clicking a row
    # sets massIndex to the row's `mass_id`.
    return Table(
        cache_id=f"mass_table_{experiment_id}",
        data=data,
        filters={SCAN_KEY: "index"},
        interactivity={MASS_KEY: "mass_id"},
        index_field="mass_id",
        column_definitions=_MASS_COLUMN_DEFINITIONS,
        go_to_fields=["mass_id"],
        title="Mass Table",
        cache_path=cache_dir,
    )


def _build_deconv_spectrum(file_manager, experiment_id: str, cache_dir: str):
    data = _lazy(file_manager, experiment_id, "deconv_spectrum_long")
    if data is None:
        return None
    # Deconvolved spectrum: filtered by scan; clicking a peak sets massIndex.
    return LinePlot(
        cache_id=f"deconv_spectrum_{experiment_id}",
        data=data,
        filters={SCAN_KEY: "index"},
        interactivity={MASS_KEY: "peak_id"},
        x_column="MonoMass",
        y_column="SumIntensity",
        title="Deconvolved Spectrum",
        x_label="Monoisotopic Mass",
        y_label="Intensity",
        cache_path=cache_dir,
    )


def _build_anno_spectrum(file_manager, experiment_id: str, cache_dir: str):
    data = _lazy(file_manager, experiment_id, "anno_spectrum_long")
    if data is None:
        return None
    # Annotated/raw spectrum: filtered by scan; consumer only (no interactivity).
    return LinePlot(
        cache_id=f"anno_spectrum_{experiment_id}",
        data=data,
        filters={SCAN_KEY: "index"},
        x_column="MonoMass_Anno",
        y_column="SumIntensity_Anno",
        title="Annotated Spectrum",
        x_label="m/z",
        y_label="Intensity",
        cache_path=cache_dir,
    )


def _build_combined_spectrum(file_manager, experiment_id: str, cache_dir: str):
    primary = _lazy(file_manager, experiment_id, "combined_spectrum_long")
    if primary is None:
        return None
    anno = _lazy(file_manager, experiment_id, "anno_spectrum_long")
    # Augmented/combined: primary deconv series + signal-peak markers, with the
    # annotated overlay supplied as the second series. The LinePlot Vue reads the
    # x2/y2 columns as INDEPENDENT column arrays (their own length), NOT row-aligned
    # with the primary series. Because the deconv peak axis and the anno peak axis
    # have different per-scan lengths, we must VERTICALLY STACK the two long frames
    # (diagonal concat) rather than relationally join them (a join would multiply
    # rows cartesian-style). After the scanIndex value-filter on `index`, the
    # primary columns are populated on the deconv rows and the anno columns on the
    # anno rows; each column array is then the correct length for its series.
    if anno is not None:
        primary = pl.concat([primary, anno], how="diagonal")
        x2, y2 = "MonoMass_Anno", "SumIntensity_Anno"
    else:
        x2 = y2 = None
    return LinePlot(
        cache_id=f"combined_spectrum_{experiment_id}",
        data=primary,
        filters={SCAN_KEY: "index"},
        interactivity={MASS_KEY: "peak_id"},
        x_column="MonoMass",
        y_column="SumIntensity",
        signal_peak_column="is_signal",
        x2_column=x2,
        y2_column=y2,
        title="Augmented Deconvolved Spectrum",
        x_label="Monoisotopic Mass",
        y_label="Intensity",
        cache_path=cache_dir,
    )


def _build_scatter3d(file_manager, experiment_id: str, cache_dir: str):
    data = _lazy(file_manager, experiment_id, "threedim_SN_plot")
    if data is None:
        return None
    # 3D S/N plot: scanIndex value-filters on `index`; massIndex handled internally
    # as an array subscript (NOT a value filter).
    return Scatter3D(
        cache_id=f"threedim_SN_plot_{experiment_id}",
        data=data,
        scan_filter="index",
        signal_column="SignalPeaks",
        noisy_column="NoisyPeaks",
        title="Precursor Signals",
        cache_path=cache_dir,
    )


def _build_fdr_plot(file_manager, experiment_id: str, cache_dir: str):
    # Precomputed {x,y} density frames stored by deconv.py. The TnT/Deconv literals
    # (axis "QScore", series "Target/Decoy QScores") are the DensityPlot defaults.
    target = _lazy(file_manager, experiment_id, "density_target")
    decoy = _lazy(file_manager, experiment_id, "density_decoy")
    if target is None and decoy is None:
        return None
    return DensityPlot(
        cache_id=f"fdr_plot_{experiment_id}",
        density_target=target,
        density_decoy=decoy,
        title="Score Distribution",
        cache_path=cache_dir,
    )


def _get_sequence(file_manager):
    """Return the submitted (sequence, fix_C, fix_M) tuple, or None."""
    if not file_manager.result_exists("sequence", "sequence"):
        return None
    sequence = file_manager.get_results("sequence", "sequence")["sequence"]
    return (
        sequence["input_sequence"],
        sequence["fixed_mod_cysteine"],
        sequence["fixed_mod_methionine"],
    )


def _build_sequence_view(file_manager, experiment_id: str, cache_dir: str):
    seq = _get_sequence(file_manager)
    if seq is None:
        return None
    sequence_string, _fix_c, _fix_m = seq
    # Deconv peaks are neutral masses (deconvolved=True). Wire the deconv long
    # spectrum as the peaks_data (renamed to the SequenceView schema: peak_id,
    # mass, intensity), filtered by the selected scan. C/M fixed mods are computed
    # from the sequence (compute_fixed_mods=True) for Deconv parity.
    peaks = _lazy(file_manager, experiment_id, "deconv_spectrum_long")
    if peaks is None:
        return None
    peaks = peaks.select(
        pl.col("index"),
        pl.col("peak_id"),
        pl.col("MonoMass").alias("mass"),
        pl.col("SumIntensity").alias("intensity"),
    )
    return SequenceView(
        cache_id=f"sequence_view_{experiment_id}",
        sequence_data=sequence_string,
        peaks_data=peaks,
        filters={SCAN_KEY: "index"},
        interactivity={MASS_KEY: "peak_id"},
        deconvolved=True,
        compute_fixed_mods=True,
        title="Sequence View",
        cache_path=cache_dir,
    )


# COMPONENT_NAMES (layout) -> builder. Mirrors FLASHDeconvLayoutManager COMPONENT_NAMES.
COMPONENT_BUILDERS = {
    "ms1_raw_heatmap": lambda fm, eid, cd: _build_heatmap(
        fm, eid, cd, "ms1_raw_heatmap", "heatmap_raw", "Raw MS1 Heatmap"),
    "ms2_raw_heatmap": lambda fm, eid, cd: _build_heatmap(
        fm, eid, cd, "ms2_raw_heatmap", "heatmap_raw2", "Raw MS2 Heatmap"),
    "ms1_deconv_heat_map": lambda fm, eid, cd: _build_heatmap(
        fm, eid, cd, "ms1_deconv_heatmap", "heatmap_deconv", "Deconvolved MS1 Heatmap"),
    "ms2_deconv_heat_map": lambda fm, eid, cd: _build_heatmap(
        fm, eid, cd, "ms2_deconv_heatmap", "heatmap_deconv2", "Deconvolved MS2 Heatmap"),
    "scan_table": _build_scan_table,
    "deconv_spectrum": _build_deconv_spectrum,
    "anno_spectrum": _build_anno_spectrum,
    "mass_table": _build_mass_table,
    "3D_SN_plot": _build_scatter3d,
    "fdr_plot": _build_fdr_plot,
    "sequence_view": _build_sequence_view,
    # internal_fragment_map: deferred (component disabled in the legacy path too).
}


def build_component(file_manager, experiment_id: str, cache_dir: str, comp_name: str):
    """Instantiate the OpenMS-Insight component for a layout cell, or None."""
    builder = COMPONENT_BUILDERS.get(comp_name)
    if builder is None:
        return None
    return builder(file_manager, experiment_id, cache_dir)


def render_experiment_panel(
    experiment_id: str,
    layout_info_per_exp: List[List[str]],
    file_manager,
    panel_index: int,
):
    """Render one experiment's [row][col] grid with its OWN isolated StateManager.

    The StateManager uses a DISTINCT session_key per experiment so selections made
    in this panel do not leak into other side-by-side panels.
    """
    session_key = f"svc_state_deconv_{experiment_id}_{panel_index}"
    state_manager = StateManager(session_key=session_key)
    cache_dir = _component_cache_dir(file_manager, experiment_id)

    for row_index, row in enumerate(layout_info_per_exp):
        columns = st.columns(len(row))
        for col, (col_index, comp_name) in zip(columns, enumerate(row)):
            with col:
                component = build_component(
                    file_manager, experiment_id, cache_dir, comp_name
                )
                # A builder returns None when its optional backing frame is
                # absent (e.g. no sequence submitted, or *_long not yet cached);
                # skip silently rather than warning on every rerun.
                if component is None:
                    continue
                key = f"deconv_oi_{panel_index}_{row_index}_{col_index}_{comp_name}"
                component(key=key, state_manager=state_manager)
