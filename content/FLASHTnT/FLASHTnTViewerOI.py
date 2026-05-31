"""FLASHTnT viewer rendered entirely with OpenMS-Insight components (Stage C).

This is the NEW viewer for the FLASHApp -> OpenMS-Insight visualization migration,
mirroring ``content/FLASHDeconv/FLASHDeconvViewerOI.py`` (Stage B). It renders the
FLASHTnT (tagger / top-down identification) workflow using the reusable
``openms_insight`` component library (``Table``, ``LinePlot``, ``SequenceView``,
``DensityPlot``, ``Heatmap``) instead of the bespoke ``flash_viewer_grid`` Vue grid
in ``src/render/*``.

Design goals (see ``/home/user/parity/STRATEGY.md`` §4/§5 and Stage C edges):

* ONE shared ``StateManager`` per rendered experiment panel, keyed by a DISTINCT
  ``session_key`` (``svc_state_tnt_<experiment_id>_<panel_index>``) so selections never
  leak between side-by-side experiment panels (HARD edge #6).
* Layout parity: the ``[experiment][row][col]`` nested grid is reproduced with
  ``st.columns`` per row (<=3 cols), rows stacked; multi-experiment side-by-side
  uses a top-level ``st.columns``.
* TnT-specific wiring (STRATEGY §2/§3):
  - ``protein_dfs`` is row-per-proteoform with ``index``; the protein Table sets
    ``proteinIndex`` on click.
  - ``tag_dfs`` is row-per-tag with ``Scan``/``ProteinIndex``/``StartPos``/``EndPos``/``mzs``.
  - The per-proteoform ``sequence_data`` store (``sequence_data_store.py``) carries
    ``coverage``/``maxCoverage`` keyed by ``proteoform_index``.
  - **Scan resolution (HARD edge #3):** a proteoform selection must resolve to the
    correct deconv scan. ``build_proteoform_scan_frame`` (additive helper in
    ``src/render/scan_resolution.py``, reproducing the legacy
    ``build_proteoform_scan_map`` PyArrow pushdown) surfaces ``proteoform_index ->
    (scan, deconv_index)`` as COLUMNS. We stamp a ``proteoform_index`` column onto
    the combined-spectrum / sequence-peak frames by joining on the deconv ``index``,
    so the OpenMS-Insight components value-filter
    (``filters={'proteinIndex': 'proteoform_index'}``) exactly like Deconv filters by
    scan.
  - **Tagger overlay (HARD edge #1):** the Tag Table sets ``tagData`` to the clicked
    tag's list of masses; the combined-spectrum ``LinePlot`` highlights peaks whose
    ``MonoMass`` matches a selected tag mass within ``abs(Δ) < 1e-5``.

NOTE: FLASHTnT runs BOTH ``parseDeconv`` and ``parseTnT`` on the same dataset
(``src/Workflow.py``), so the Deconv long-format frames (``combined_spectrum_long``,
``scan_table``, heatmaps) are present alongside the TnT frames (``protein_dfs``,
``tag_dfs``, ``sequence_data``, ``settings``, ``density_id_target``/``density_id_decoy``).

The OLD render path (``src/render/render.py`` / ``flash_viewer_grid``) is left intact
and importable; the page chooses which path to use via ``use_openms_insight_viewer``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, List, Optional

import polars as pl
import streamlit as st

from openms_insight import (
    DensityPlot,
    Heatmap,
    LinePlot,
    SequenceView,
    StateManager,
    Table,
)

from src.render.scan_resolution import build_proteoform_scan_frame

# FLASHApp StateTracker keys reused as OpenMS-Insight identifiers so state flows
# across components exactly like the legacy grid.
PROTEIN_KEY = "proteinIndex"
# Tag selection: the Tag Table sets a SCALAR `tagData` to the clicked tag's
# `TagIndex` (a scalar — list-valued interactivity columns are not supported by
# the OpenMS-Insight Table, which calls `.item()` on the cell). The viewer then
# resolves that TagIndex to the tag's list of masses and publishes it under
# `TAG_MASSES_KEY`, which the combined-spectrum LinePlot consumes for the tagger
# overlay (`tag_filters={'tagMasses': 'MonoMass'}`).
TAG_KEY = "tagData"
TAG_MASSES_KEY = "tagMasses"
MASS_KEY = "massIndex"


def _component_cache_dir(file_manager, experiment_id: str) -> str:
    """Directory under the workspace cache where OI component caches are written."""
    cache_root = Path(file_manager.cache_path, "oi_components_tnt", str(experiment_id))
    cache_root.mkdir(parents=True, exist_ok=True)
    return str(cache_root)


def _lazy(file_manager, experiment_id: str, name_tag: str) -> Optional[pl.LazyFrame]:
    """Load a stored frame as a polars LazyFrame, or None if absent."""
    if not file_manager.result_exists(experiment_id, name_tag):
        return None
    return file_manager.get_results(
        experiment_id, [name_tag], use_polars=True
    )[name_tag]


def _pandas(file_manager, experiment_id: str, name_tag: str):
    """Load a stored frame as pandas (for the precomputed density frames), or None."""
    if not file_manager.result_exists(experiment_id, name_tag):
        return None
    return file_manager.get_results(experiment_id, [name_tag])[name_tag]


# ---------------------------------------------------------------------------
# Scan resolution: proteoform_index -> deconv index, exposed as a frame so the
# spectrum / sequence components can value-filter by proteoform.
# ---------------------------------------------------------------------------

def _proteoform_scan_frame(file_manager, experiment_id: str) -> Optional[pl.DataFrame]:
    """proteoform_index / scan / deconv_index frame for this experiment, or None.

    Reproduces the legacy ``build_proteoform_scan_map`` (PyArrow pushdown in
    ``src/render/update.py``) by reading the already-stored ``protein_dfs`` and
    ``scan_table`` frames. Cached in session state per experiment to avoid
    recomputing on every rerun.
    """
    protein = _lazy(file_manager, experiment_id, "protein_dfs")
    scan_table = _lazy(file_manager, experiment_id, "scan_table")
    if protein is None or scan_table is None:
        return None
    protein_df = protein.select(["index", "Scan"]).collect()
    scan_df = scan_table.select(["index", "Scan"]).collect()
    return build_proteoform_scan_frame(protein_df, scan_df)


def _stamp_proteoform_index(
    spectrum_lf: pl.LazyFrame, scan_frame: pl.DataFrame
) -> pl.LazyFrame:
    """Join a deconv-``index``-keyed long spectrum frame with the proteoform/scan
    frame so each peak row carries the ``proteoform_index`` whose scan it belongs
    to. This converts the proteoform selection into a plain value filter on the
    spectrum (``filters={'proteinIndex': 'proteoform_index'}``).

    A scan may map to multiple proteoforms; the inner join replicates the peak
    rows per proteoform so each proteoform selection sees its scan's peaks (the
    legacy path resolves proteoform->scan then pushes that single scan down, which
    is equivalent for the selected proteoform)."""
    map_lf = scan_frame.lazy().select(
        pl.col("deconv_index").alias("index"),
        pl.col("proteoform_index"),
    )
    return spectrum_lf.join(map_lf, on="index", how="inner")


# ---------------------------------------------------------------------------
# Per-component builders. Each returns an OpenMS-Insight component instance, or
# None when the underlying data frame is missing (component is silently skipped).
# ---------------------------------------------------------------------------

def _build_protein_table(file_manager, experiment_id: str, cache_dir: str):
    data = _lazy(file_manager, experiment_id, "protein_dfs")
    if data is None:
        return None
    # Protein table: clicking a row sets proteinIndex to the row's `index`.
    return Table(
        cache_id=f"protein_table_{experiment_id}",
        data=data,
        interactivity={PROTEIN_KEY: "index"},
        index_field="index",
        title="Protein Table",
        cache_path=cache_dir,
    )


def _build_tag_table(file_manager, experiment_id: str, cache_dir: str):
    data = _lazy(file_manager, experiment_id, "tag_dfs")
    if data is None:
        return None
    scan_frame = _proteoform_scan_frame(file_manager, experiment_id)
    if scan_frame is None:
        return None
    # Tags are scan (spectrum) data. To filter by the SELECTED PROTEOFORM we need
    # a proteoform_index column on each tag row; resolve via the proteoform's scan
    # (Scan column on the tag) so a proteoform selection shows its scan's tags
    # (parity with the legacy filter_data Tag-Table path stamping ProteinIndex).
    map_lf = scan_frame.lazy().select(
        pl.col("scan").alias("Scan"),
        pl.col("proteoform_index"),
    )
    tag_lf = data.join(map_lf, on="Scan", how="inner")
    # Clicking a tag row sets the SCALAR `tagData` to the row's `TagIndex`. The
    # viewer resolves that index to the tag's masses (see _resolve_tag_masses) and
    # publishes them for the combined-spectrum tagger overlay. A list-valued
    # interactivity column cannot be used here because the OI Table calls
    # `.item()` on the clicked cell.
    return Table(
        cache_id=f"tag_table_{experiment_id}",
        data=tag_lf,
        filters={PROTEIN_KEY: "proteoform_index"},
        interactivity={TAG_KEY: "TagIndex"},
        index_field="TagIndex",
        title="Tag Table",
        cache_path=cache_dir,
    )


def _tag_mass_lookup(file_manager, experiment_id: str) -> dict:
    """Map ``TagIndex`` -> list[float] of the tag's masses (parsed from the
    comma-joined ``mzs`` string with its trailing comma). Used to resolve a tag
    selection into the mass list the combined-spectrum tagger overlay matches."""
    tags = _lazy(file_manager, experiment_id, "tag_dfs")
    if tags is None:
        return {}
    df = (
        tags.select(["TagIndex", "mzs"])
        .with_columns(
            pl.col("mzs")
            .str.strip_chars(",")
            .str.split(",")
            .list.eval(pl.element().cast(pl.Float64, strict=False))
            .alias("tag_masses")
        )
        .collect()
    )
    return {
        int(ti): [m for m in masses if m is not None]
        for ti, masses in zip(df["TagIndex"], df["tag_masses"].to_list())
    }


def _resolve_tag_masses(file_manager, experiment_id: str, state_manager) -> None:
    """Resolve the selected ``tagData`` (a ``TagIndex``) to its list of masses and
    publish under ``tagMasses`` so the combined-spectrum LinePlot tagger overlay
    can read it. Clears ``tagMasses`` when no tag is selected."""
    tag_index = state_manager.get_selection(TAG_KEY)
    if tag_index is None:
        state_manager.clear_selection(TAG_MASSES_KEY)
        return
    lookup = _tag_mass_lookup(file_manager, experiment_id)
    masses = lookup.get(int(tag_index))
    if masses:
        state_manager.set_selection(TAG_MASSES_KEY, list(masses))
    else:
        state_manager.clear_selection(TAG_MASSES_KEY)


def _build_sequence_frame(
    file_manager, experiment_id: str
) -> Optional[pl.LazyFrame]:
    """Build the SequenceView-ready per-proteoform sequence frame.

    Source: the per-proteoform ``sequence_data`` store (keyed by
    ``proteoform_index``, carrying the per-residue ``coverage`` / ``maxCoverage``
    of the DISPLAYED proteoform substring and the full-protein ``sequence`` list +
    ``proteoform_start``/``proteoform_end``). We reconstruct the displayed
    proteoform sequence STRING (the substring the legacy SequenceView rendered)
    and attach coverage so OpenMS-Insight SequenceView can shade residues.

    Columns emitted: ``proteoform_index`` (filter key), ``sequence`` (str),
    ``precursor_charge`` (=1, neutral/deconvolved peaks), ``coverage`` (list[f64]),
    ``maxCoverage`` (f64), ``fixed_modifications`` (list[str])."""
    # ``sequence_data`` is a pickle-backed store: a dict keyed by
    # ``proteoform_index``, each value a dict with per-residue ``sequence`` /
    # ``coverage`` lists (full protein), ``maxCoverage``, ``proteoform_start`` /
    # ``proteoform_end`` and ``fixed_modifications``. It is NOT a tabular frame,
    # so it is loaded as the raw object (``_pandas`` returns the unpickled dict)
    # and iterated — loading it as a LazyFrame raises AttributeError and leaves
    # SequenceView blank.
    store = _pandas(file_manager, experiment_id, "sequence_data")
    if not isinstance(store, dict) or not store:
        return None

    proteoform_indices: List[int] = []
    sequences: List[str] = []
    coverages: List[list] = []
    max_coverages: List[float] = []
    fixed_mods: List[list] = []
    for pid in sorted(store.keys()):
        entry = store[pid] or {}
        full = list(entry.get("sequence") or [])
        cov = list(entry.get("coverage") or [])
        start = entry.get("proteoform_start")
        end = entry.get("proteoform_end")
        # Slice the displayed proteoform substring AND its coverage together so the
        # two stay aligned (the legacy SequenceView rendered the substring). A
        # negative/absent bound means render the full protein.
        if start is None or end is None or start < 0 or end < 0:
            sub_seq, sub_cov = full, cov
        else:
            sub_seq, sub_cov = full[start:end + 1], cov[start:end + 1]
        proteoform_indices.append(int(pid))
        sequences.append("".join(sub_seq))
        coverages.append([float(c) for c in sub_cov])
        mc = entry.get("maxCoverage")
        max_coverages.append(float(mc) if mc is not None else 0.0)
        fm = entry.get("fixed_modifications") or []
        fixed_mods.append([str(m) for m in fm])

    out = pl.DataFrame({
        "proteoform_index": proteoform_indices,
        "sequence": sequences,
        "precursor_charge": [1] * len(proteoform_indices),
        "coverage": coverages,
        "maxCoverage": max_coverages,
        "fixed_modifications": fixed_mods,
    })
    return out.lazy()


def _build_sequence_view(file_manager, experiment_id: str, cache_dir: str):
    seq_frame = _build_sequence_frame(file_manager, experiment_id)
    if seq_frame is None:
        return None
    scan_frame = _proteoform_scan_frame(file_manager, experiment_id)
    combined = _lazy(file_manager, experiment_id, "combined_spectrum_long")
    if combined is None:
        combined = _lazy(file_manager, experiment_id, "deconv_spectrum_long")
    peaks = None
    if combined is not None and scan_frame is not None:
        # Deconv peaks are neutral masses; filter by the proteoform's scan and
        # rename to the SequenceView peaks schema (peak_id, mass, intensity).
        peaks = (
            _stamp_proteoform_index(combined, scan_frame)
            .select(
                pl.col("proteoform_index"),
                pl.col("peak_id"),
                pl.col("MonoMass").alias("mass"),
                pl.col("SumIntensity").alias("intensity"),
            )
        )

    settings = _pandas(file_manager, experiment_id, "settings")
    settings = dict(settings) if isinstance(settings, dict) else None

    return SequenceView(
        cache_id=f"sequence_view_{experiment_id}",
        sequence_data=seq_frame,
        peaks_data=peaks,
        filters={PROTEIN_KEY: "proteoform_index"},
        interactivity={MASS_KEY: "peak_id"},
        deconvolved=True,
        compute_fixed_mods=True,
        settings=settings,
        title="Sequence View",
        cache_path=cache_dir,
    )


def _build_combined_spectrum(file_manager, experiment_id: str, cache_dir: str):
    primary = _lazy(file_manager, experiment_id, "combined_spectrum_long")
    if primary is None:
        return None
    scan_frame = _proteoform_scan_frame(file_manager, experiment_id)
    if scan_frame is None:
        return None
    primary = _stamp_proteoform_index(primary, scan_frame)

    # Annotated overlay (2nd series), stamped + filtered by the same proteoform.
    anno = _lazy(file_manager, experiment_id, "anno_spectrum_long")
    if anno is not None:
        anno = _stamp_proteoform_index(anno, scan_frame)
        primary = pl.concat([primary, anno], how="diagonal")
        x2, y2 = "MonoMass_Anno", "SumIntensity_Anno"
    else:
        x2 = y2 = None

    # Combined spectrum: filtered by proteoform (resolved to scan), clicking a
    # peak sets massIndex, signal-peak markers via is_signal, and the TAGGER
    # OVERLAY highlights peaks whose MonoMass matches a selected tag mass
    # (abs(Δ) < 1e-5, FLASHApp PlotlyLineplotTagger parity). The selected tag's
    # masses arrive via the `tagData` state value (a list) set by the Tag Table.
    return LinePlot(
        cache_id=f"combined_spectrum_{experiment_id}",
        data=primary,
        filters={PROTEIN_KEY: "proteoform_index"},
        interactivity={MASS_KEY: "peak_id"},
        x_column="MonoMass",
        y_column="SumIntensity",
        signal_peak_column="is_signal",
        x2_column=x2,
        y2_column=y2,
        tag_filters={TAG_MASSES_KEY: "MonoMass"},
        tag_mass_column="MonoMass",
        tag_tolerance=1e-5,
        title="Augmented Deconvolved Spectrum",
        x_label="Monoisotopic Mass",
        y_label="Intensity",
        cache_path=cache_dir,
    )


def _build_id_fdr_plot(file_manager, experiment_id: str, cache_dir: str):
    # Precomputed TnT id-FDR density frames (computed in parseTnT with the TnT
    # grouping: DECOY_ accession + ProteoformLevelQvalue>0). Literal labels stay
    # "QScore"/"Target QScores"/"Decoy QScores" (DensityPlot defaults).
    target = _lazy(file_manager, experiment_id, "density_id_target")
    decoy = _lazy(file_manager, experiment_id, "density_id_decoy")
    if target is None and decoy is None:
        return None
    return DensityPlot(
        cache_id=f"id_fdr_plot_{experiment_id}",
        density_target=target,
        density_decoy=decoy,
        title="Score Distribution",
        cache_path=cache_dir,
    )


def _build_heatmap(
    file_manager, experiment_id: str, cache_dir: str, frame_tag: str,
    zoom_id: str, title: str,
):
    data = _lazy(file_manager, experiment_id, frame_tag)
    if data is None:
        return None
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


# COMPONENT_NAMES (FLASHTnTLayoutManager) -> builder.
COMPONENT_BUILDERS = {
    "protein_table": _build_protein_table,
    "sequence_view": _build_sequence_view,
    "tag_table": _build_tag_table,
    "combined_spectrum": _build_combined_spectrum,
    "id_fdr_plot": _build_id_fdr_plot,
    "ms1_raw_heatmap": lambda fm, eid, cd: _build_heatmap(
        fm, eid, cd, "ms1_raw_heatmap", "heatmap_raw", "Raw MS1 Heatmap"),
    "ms1_deconv_heat_map": lambda fm, eid, cd: _build_heatmap(
        fm, eid, cd, "ms1_deconv_heatmap", "heatmap_deconv", "Deconvolved MS1 Heatmap"),
    # internal_fragment_map: deferred (disabled in the legacy path too; the
    # sequence_data store would need internal-fragment arrays — see module note).
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
    in this panel do not leak into other side-by-side panels (HARD edge #6).
    """
    session_key = f"svc_state_tnt_{experiment_id}_{panel_index}"
    state_manager = StateManager(session_key=session_key)
    cache_dir = _component_cache_dir(file_manager, experiment_id)

    # Resolve the selected tag (scalar TagIndex set by the Tag Table) into its
    # list of masses BEFORE rendering so the combined-spectrum tagger overlay
    # sees the up-to-date `tagMasses` selection this rerun.
    _resolve_tag_masses(file_manager, experiment_id, state_manager)

    for row_index, row in enumerate(layout_info_per_exp):
        columns = st.columns(len(row))
        for col, (col_index, comp_name) in zip(columns, enumerate(row)):
            with col:
                component = build_component(
                    file_manager, experiment_id, cache_dir, comp_name
                )
                if component is None:
                    st.warning(f"No data for '{comp_name}'.")
                    continue
                key = f"tnt_oi_{panel_index}_{row_index}_{col_index}_{comp_name}"
                component(key=key, state_manager=state_manager)
