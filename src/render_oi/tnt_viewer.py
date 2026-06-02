"""OpenMS-Insight rendering engine for FLASHTnT (migration Phase 2).

FLASHTnT's master selection is ``proteinIndex`` (set by the protein table). The
challenge the old ``update.py`` solved: spectrum / mass / sequence panels are
keyed by the deconv ROW INDEX, while the tag table is keyed by ``proteinIndex``
directly. We reproduce that here with a resolver:

    protein table click → proteinIndex
    resolver (build_proteoform_scan_map) → deconvIndex  (= scan's deconv row)
    spectra / mass / sequence  filter by {deconvIndex: index}
    tag table                  filter by {proteinIndex: ProteinIndex}
    sequence view              filter by {proteinIndex: proteoform_index}

The resolver runs after the protein table is rendered (its interactivity sets
``proteinIndex``) and stamps ``deconvIndex`` into the same StateManager so the
downstream value-filters reproduce the original PyArrow pushdown.
"""

from __future__ import annotations

import gzip
import logging
import pickle
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

import polars as pl

from .deconv_viewer import (
    _HEATMAP_SPEC,
    _load_pandas,
    _load_pandas_pl,
    _load_polars,
    _oi_cache_dir,
)

logger = logging.getLogger(__name__)

# State identifiers
PROTEIN = "proteinIndex"
DECONV = "deconvIndex"


def _load_pickle_gz(path: Path):
    with gzip.open(path, "rb") as f:
        return pickle.load(f)


def _build_proteoform_scan_map(file_manager, dataset_id: str) -> Dict[int, Dict[str, int]]:
    """proteinIndex → {scan, deconv_index} using the existing resolver."""
    from src.render.scan_resolution import build_proteoform_scan_map

    prot = _load_pandas(file_manager, dataset_id, "protein_dfs")
    scan = _load_pandas(file_manager, dataset_id, "scan_table")
    return build_proteoform_scan_map(prot[["index", "Scan"]], scan[["index", "Scan"]])


# Ion types whose precomputed per-residue fragment masses are forwarded to
# SequenceView (the proteoform entry carries a/b/c/x/y/z).
_ION_TYPES = ("a", "b", "c", "x", "y", "z")


def _observed_mass(value):
    """ProteoformMass, or None when FLASHTnT's -1.0 'unmatched' sentinel is set
    (so SequenceView omits the observed / Δ-mass header)."""
    try:
        mass = float(value)
    except (TypeError, ValueError):
        return None
    return mass if mass > 0 else None


def _frag_masses(value):
    """Per-residue precomputed fragment masses as list[list[float]] (the inner list
    holds modification-ambiguity variants); tolerates a flat list[float]."""
    if not value:
        return []
    out = []
    for residue in value:
        if isinstance(residue, (list, tuple)):
            out.append([float(x) for x in residue])
        else:
            out.append([float(residue)])
    return out


def _peaks_table(
    file_manager, dataset_id: str, scan_map: Dict[int, Dict[str, int]]
) -> Optional[pl.LazyFrame]:
    """Observed deconvolved peaks per proteoform, for SequenceView fragment matching.

    SequenceView matches the precomputed fragment masses against observed peaks. We
    supply each proteoform's scan peaks (neutral ``MonoMass`` + ``SumIntensity`` from
    ``deconv_spectrum``) stamped with ``proteoform_index`` so the component's
    ``filters`` select the right peaks for the current selection.
    """
    pairs = [
        {"proteoform_index": int(pid), "deconv_index": int(e["deconv_index"])}
        for pid, e in (scan_map or {}).items()
        if e.get("deconv_index") is not None
    ]
    if not pairs:
        return None
    spec = _load_polars(file_manager, dataset_id, "deconv_spectrum")
    long = (
        spec.select(["index", "MonoMass", "SumIntensity"])
        .explode(["MonoMass", "SumIntensity"])
        .rename({"MonoMass": "mass", "SumIntensity": "intensity"})
        .with_columns(pl.int_range(pl.len()).over("index").alias("peak_id"))
    )
    mapping = pl.DataFrame(pairs).lazy()
    return long.join(
        mapping, left_on="index", right_on="deconv_index", how="inner"
    ).select(["proteoform_index", "peak_id", "mass", "intensity"])


def _sequence_table(file_manager, dataset_id: str) -> Optional[pl.LazyFrame]:
    """Build a one-row-per-proteoform sequence frame for SequenceView.

    Reads the cached ``sequence_data`` ({pid: entry}); each entry already carries
    sequence/coverage/maxCoverage/fragment masses. We emit a LazyFrame with a
    ``proteoform_index`` filter column plus ``sequence`` (joined string),
    ``precursor_charge``, and per-proteoform coverage arrays so OI's extended
    SequenceView can shade residues.
    """
    res = file_manager.get_results(
        dataset_id, ["sequence_data"], partial=True, use_pyarrow=True
    )
    if "sequence_data" not in res:
        return None
    p = res["sequence_data"]
    # A fresh FLASHTnT run stores sequence_data as a parquet dataset; with
    # use_pyarrow=True FileManager hands back a pyarrow Dataset (without it we'd get
    # a pandas DataFrame that reconstruct_all can't read — the cause of the
    # "Component unavailable: sequence_view" regression). Older/example caches store
    # it as a pickled {pid: entry} dict (FileManager unpickles .pkl.gz for us, but a
    # path may also be handed back). Let failures propagate so they surface as a
    # render error rather than a silent None.
    if isinstance(p, Path) and p.suffix == ".gz":
        data = _load_pickle_gz(p)
    elif isinstance(p, dict):
        data = p
    else:
        from src.render.sequence_data_store import reconstruct_all

        data = reconstruct_all(p)

    rows = []
    for pid in sorted(data):
        entry = data[pid]
        seq = entry.get("sequence") or []
        row = {
            "proteoform_index": int(pid),
            "sequence": "".join(seq) if isinstance(seq, list) else str(seq),
            "precursor_charge": 1,
            "coverage": [float(c) for c in (entry.get("coverage") or [])],
            "max_coverage": float(entry.get("maxCoverage") or 0.0),
            # Header masses: theoretical from the sequence, observed = ProteoformMass
            # (FLASHTnT stores -1.0 when unmatched -> None so the header omits it).
            "theoretical_mass": float(entry.get("theoretical_mass") or 0.0),
            "observed_mass": _observed_mass(entry.get("computed_mass")),
        }
        # Precomputed per-residue fragment masses (account for proteoform mods);
        # SequenceView matches these against the scan peaks instead of recomputing.
        for ion in _ION_TYPES:
            row[f"fragment_masses_{ion}"] = _frag_masses(
                entry.get(f"fragment_masses_{ion}")
            )
        rows.append(row)
    if not rows:
        return None
    return (
        pl.DataFrame(rows)
        .with_columns(pl.col("observed_mass").cast(pl.Float64, strict=False))
        .lazy()
    )


def build_component_tnt(
    comp_name: str,
    dataset_id: str,
    file_manager,
    state_manager,
    key_prefix: str,
) -> Optional[Callable[[], Any]]:
    """Build a render callable for one FLASHTnT layout cell."""
    from openms_insight import DensityPlot, Heatmap, LinePlot, SequenceView, Table

    from src.parse.long_format import (
        density_series_long,
        explode_combined_spectrum_long,
    )

    cache_dir = _oi_cache_dir(file_manager, dataset_id)
    cid = lambda name: f"{dataset_id}__tnt_{name}"  # noqa: E731
    skey = lambda name: f"{key_prefix}_{name}"  # noqa: E731

    # ---- Protein table (master; click sets proteinIndex) ----
    if comp_name == "protein_table":
        data = _load_polars(file_manager, dataset_id, "protein_dfs")
        # Keep the informative columns the original Protein Table showed.
        keep = [
            "index",
            "accession",
            "description",
            "ProteoformMass",
            "Coverage(%)",
            "TagCount",
            "ProteoformLevelQvalue",
        ]
        schema = data.collect_schema().names()
        cols = [c for c in keep if c in schema]
        tbl = Table(
            cache_id=cid("protein_table"),
            data=data.select(cols) if cols else data,
            interactivity={PROTEIN: "index"},
            index_field="index",
            title="Protein Table",
            cache_path=cache_dir,
        )
        return lambda: tbl(key=skey("protein_table"), state_manager=state_manager)

    # ---- Tag table (filtered by proteinIndex) ----
    if comp_name == "tag_table":
        data = _load_polars(file_manager, dataset_id, "tag_dfs")
        keep = [
            "TagIndex",
            "TagSequence",
            "StartPos",
            "EndPos",
            "Length",
            "Score",
            "DeltaMass",
            "ProteinIndex",
        ]
        schema = data.collect_schema().names()
        cols = [c for c in keep if c in schema]
        tbl = Table(
            cache_id=cid("tag_table"),
            data=data.select(cols) if cols else data,
            filters={PROTEIN: "ProteinIndex"},
            index_field="TagIndex",
            title="Tag Table",
            cache_path=cache_dir,
        )
        return lambda: tbl(key=skey("tag_table"), state_manager=state_manager)

    # ---- Combined / augmented spectrum (deconv primary + annotated overlay) ----
    if comp_name == "combined_spectrum":
        per_scan = _load_polars(file_manager, dataset_id, "combined_spectrum")
        deconv_long, anno_long = explode_combined_spectrum_long(per_scan)
        lp = LinePlot(
            cache_id=cid("combined_spectrum"),
            data=deconv_long,
            overlay_data=anno_long,
            filters={DECONV: "index"},
            x_column="mass",
            y_column="intensity",
            overlay_x_column="mass",
            overlay_y_column="intensity",
            overlay_name="Annotated",
            title="Augmented Deconvolved Spectrum",
            x_label="Monoisotopic Mass",
            y_label="Intensity",
            cache_path=cache_dir,
        )
        return lambda: lp(key=skey("combined_spectrum"), state_manager=state_manager)

    # ---- Sequence view (filtered by proteinIndex; coverage coloring) ----
    if comp_name == "sequence_view":
        seq_tbl = _sequence_table(file_manager, dataset_id)
        if seq_tbl is None:
            return None
        settings = _tnt_settings(file_manager, dataset_id)
        scan_map = _build_proteoform_scan_map(file_manager, dataset_id)
        peaks_tbl = _peaks_table(file_manager, dataset_id, scan_map)
        sv = SequenceView(
            cache_id=cid("sequence_view"),
            sequence_data=seq_tbl,
            filters={PROTEIN: "proteoform_index"},
            deconvolved=True,
            coverage_column="coverage",
            max_coverage_column="max_coverage",
            # Header masses + precomputed fragment-ion masses (so b/y flags and the
            # matching-fragments table reflect the modified proteoform), matched
            # against the proteoform's observed scan peaks.
            theoretical_mass_column="theoretical_mass",
            observed_mass_column="observed_mass",
            fragment_mass_columns={
                ion: f"fragment_masses_{ion}" for ion in _ION_TYPES
            },
            peaks_data=peaks_tbl,
            interactivity={"peak": "peak_id"},
            annotation_config={
                "ion_types": settings.get("ion_types", ["b", "y"]),
                "tolerance": settings.get("tolerance", 10.0),
                "tolerance_ppm": True,
            },
            title="Sequence View",
            cache_path=cache_dir,
        )
        return lambda: sv(key=skey("sequence_view"), state_manager=state_manager)

    # ---- Identification FDR / score distribution (DensityPlot, precomputed) ----
    if comp_name == "id_fdr_plot":
        target = _load_pandas_pl(file_manager, dataset_id, "density_id_target")
        decoy = _load_pandas_pl(file_manager, dataset_id, "density_id_decoy")
        long = density_series_long(target, decoy)
        dp = DensityPlot(
            cache_id=cid("id_fdr_plot"),
            data=long.lazy(),
            precomputed=True,
            series_column="series",
            series_config={
                "Target": {"label": "Target", "color": "green"},
                "Decoy": {"label": "Decoy", "color": "red"},
            },
            title="Score Distribution",
            x_label="Proteoform-level q-value",
            cache_path=cache_dir,
        )
        return lambda: dp(key=skey("id_fdr_plot"), state_manager=state_manager)

    # ---- Heatmaps (reuse Deconv spec) ----
    if comp_name in _HEATMAP_SPEC:
        title, cache_name = _HEATMAP_SPEC[comp_name]
        data = _load_polars(file_manager, dataset_id, cache_name)
        hm = Heatmap(
            cache_id=cid(comp_name),
            data=data,
            x_column="rt",
            y_column="mass",
            intensity_column="intensity",
            title=title,
            x_label="Retention time",
            y_label="Monoisotopic mass",
            zoom_identifier=f"tnt_{comp_name}_zoom",
            cache_path=cache_dir,
        )
        return lambda: hm(key=skey(comp_name), state_manager=state_manager)

    logger.warning("Unknown / unavailable TnT component: %s", comp_name)
    return None


def _tnt_settings(file_manager, dataset_id: str) -> Dict[str, Any]:
    """Load the cached FLASHTnT settings ({tolerance, ion_types})."""
    res = file_manager.get_results(dataset_id, ["settings"], partial=True)
    s = res.get("settings")
    if isinstance(s, Path) and s.suffix == ".gz":
        try:
            return _load_pickle_gz(s)
        except Exception:
            return {}
    return s if isinstance(s, dict) else {}


def render_experiment_tnt(
    dataset_id: str,
    layout_rows: List[List[str]],
    file_manager,
    panel_key: str,
) -> None:
    """Render one FLASHTnT experiment panel with proteoform→scan resolution.

    A per-experiment StateManager keeps side-by-side panels isolated. After the
    protein table sets ``proteinIndex``, we resolve it to ``deconvIndex`` so the
    spectrum / mass / sequence panels filter by the deconv row index — exactly
    the proteoform→scan resolution the original update.py performed.
    """
    import streamlit as st
    from openms_insight import StateManager

    state_manager = StateManager(session_key=f"oi_tnt_state_{panel_key}")

    # Resolve proteinIndex → deconvIndex BEFORE rendering downstream panels so
    # the spectrum/mass/sequence filters see the right scan on this run.
    protein_index = state_manager.get_selection(PROTEIN)
    if protein_index is not None:
        scan_map = _build_proteoform_scan_map(file_manager, dataset_id)
        entry = scan_map.get(int(protein_index))
        deconv_index = entry["deconv_index"] if entry else None
        if state_manager.get_selection(DECONV) != deconv_index:
            state_manager.set_selection(DECONV, deconv_index)

    for row_index, row in enumerate(layout_rows):
        if not row:
            continue
        cols = st.columns(len(row))
        for col_index, comp_name in enumerate(row):
            with cols[col_index]:
                try:
                    render = build_component_tnt(
                        comp_name,
                        dataset_id,
                        file_manager,
                        state_manager,
                        key_prefix=f"{panel_key}_{row_index}_{col_index}",
                    )
                    if render is not None:
                        render()
                    else:
                        st.warning(f"Component unavailable: {comp_name}")
                except Exception as exc:  # pragma: no cover - defensive UI guard
                    logger.exception("Failed to render %s", comp_name)
                    st.error(f"Error rendering {comp_name}: {exc}")
