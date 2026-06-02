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
TAG = "tagIndex"
AAPOS = "AApos"


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


def _tag_data(file_manager, dataset_id: str, tag_index, aa_pos=None) -> Optional[dict]:
    """Build the tagger-overlay ``tagData`` for the selected tag row.

    Mirrors the legacy tag-table click: parse the comma-joined fragment ``mzs``,
    carry the tag sequence + span, and flag N-terminal tags (``Nmass == -1``).
    When a sequence-view residue is selected (``aa_pos``) and it falls within the
    tag span, ``selectedAA`` is its tag-relative offset (drives the gold
    selected-residue highlight); otherwise it stays unset (-1000).
    """
    tags = (
        _load_polars(file_manager, dataset_id, "tag_dfs")
        .filter(pl.col("TagIndex") == int(tag_index))
        .collect()
    )
    if tags.height == 0:
        return None
    r = tags.row(0, named=True)
    masses = [
        float(m)
        for m in str(r.get("mzs") or "").split(",")
        if m.strip() and float(m) != 0.0
    ]
    start = int(r.get("StartPos") or 0)
    end = int(r.get("EndPos") or 0)
    # Gold selected-residue highlight: when a sequence-view residue (aa_pos) is
    # selected and lies within this tag's span, selectedAA is its tag-relative
    # offset (legacy: selectedAApos - StartPos); otherwise unset (-1000).
    selected_aa = -1000
    if aa_pos is not None and start <= int(aa_pos) <= end:
        selected_aa = int(aa_pos) - start
    return {
        "masses": masses,
        "sequence": str(r.get("TagSequence") or ""),
        "nTerminal": float(r.get("Nmass", -1) or -1) == -1,
        "startPos": start,
        "endPos": end,
        "selectedAA": selected_aa,
    }


def _deconv_signal_peaks_long(per_scan: pl.LazyFrame) -> pl.LazyFrame:
    """Primary deconvolved spectrum with a per-peak ``signal_peaks`` column.

    Each row (one deconvolved mass) carries its constituent peaks as
    ``[mz, intensity, charge]`` triplets (``combined_spectrum.SignalPeaks`` stores
    ``[binIdx, mz, intensity, charge]``; we drop the bin index), aligned 1:1 with
    the primary sticks so LinePlot's tagger overlay can draw the per-charge buttons.
    """
    return (
        per_scan.select(["index", "MonoMass", "SumIntensity", "SignalPeaks"])
        .explode(["MonoMass", "SumIntensity", "SignalPeaks"])
        .with_columns(
            pl.col("SignalPeaks")
            .list.eval(pl.element().list.slice(1, 3))
            .alias("signal_peaks")
        )
        .rename({"MonoMass": "mass", "SumIntensity": "intensity"})
        .select(["index", "mass", "intensity", "signal_peaks"])
    )


# Tabulator float formatter: fixed-precision display standing in for the legacy
# ``toFixed(4)`` (matches the Deconv viewer's _FLOAT_FMT convention).
_FLOAT_FMT = {"formatter": "money", "formatterParams": {"precision": 4, "symbol": ""}}

# Protein-table column spec mirroring the legacy ``TabulatorProteinTable``
# (title / headerTooltip / sorter) for the REAL ``protein_dfs`` fields. Each entry
# is (field, title, tooltip, is_float, dash_sentinel). ``dash_sentinel`` flags the
# columns the legacy rendered with a "-1 -> '-'" formatter (ProteoformMass,
# ProteoformLevelQvalue) -- the FLASHTnT "unmatched" sentinel; we null those -1.0s
# in the data (see _apply_dash_sentinels) so the cell renders blank instead of a
# misleading -1, keeping the column numeric/sortable. ``Coverage(%)`` is added per
# the parity request (it was not in the legacy protein table).
_PROTEIN_TABLE_COLUMNS = [
    (
        "Scan",
        "Scan No.",
        "The identifier of the mass spectrometry scan associated with the "
        "identified proteoform.",
        False,
        False,
    ),
    (
        "accession",
        "Accession",
        "The unique identifier for the protein in the reference database.",
        False,
        False,
    ),
    (
        "description",
        "Description",
        "A human-readable description of the matched protein.",
        False,
        False,
    ),
    (
        "length",
        "Length",
        "The total number of amino acids in the matched protein.",
        False,
        False,
    ),
    (
        "ProteoformMass",
        "Mass",
        "The calculated mass of the proteoform in Daltons.",
        True,
        True,
    ),
    (
        "Coverage(%)",
        "Coverage (%)",
        "The percentage of the protein sequence covered by matched fragments.",
        True,
        False,
    ),
    (
        "MatchingFragments",
        "No. of Matched Fragments",
        "The number of fragment ions that match the protein sequence.",
        False,
        False,
    ),
    (
        "ModCount",
        "No. of Modifications",
        "The number of modifications identified in the protein.",
        False,
        False,
    ),
    (
        "TagCount",
        "No. of Tags",
        "The number of sequence tags associated with the proteoform match.",
        False,
        False,
    ),
    (
        "Score",
        "Score",
        "A score indicating the confidence of the protein match (higher is "
        "better).",
        False,
        False,
    ),
    (
        "ProteoformLevelQvalue",
        "Q-Value (Proteoform Level)",
        "The confidence value of the protein match at the proteoform level.",
        True,
        True,
    ),
]

# Tag-table column spec mirroring the legacy ``TabulatorTagTable``. ``Nmass`` and
# ``Cmass`` carry the legacy "-1 -> '-'" sentinel (N-/C-terminal offset absent).
_TAG_TABLE_COLUMNS = [
    (
        "Scan",
        "Scan Number",
        "The identifier of the mass spectrometry scan containing the sequence "
        "tag.",
        False,
        False,
    ),
    (
        "StartPos",
        "Start Position",
        "The position in the protein sequence where the sequence tag begins.",
        False,
        False,
    ),
    (
        "EndPos",
        "End Position",
        "The position in the protein sequence where the sequence tag ends.",
        False,
        False,
    ),
    (
        "TagSequence",
        "Sequence",
        "The amino acid sequence of the identified tag.",
        False,
        False,
    ),
    (
        "Length",
        "Length",
        "The number of amino acids in the sequence tag.",
        False,
        False,
    ),
    (
        "Score",
        "Tag Score",
        "A score indicating the confidence of the sequence tag identification "
        "(higher is better).",
        False,
        False,
    ),
    (
        "Nmass",
        "N mass",
        "The N-terminal mass offset from the start of the sequence tag in "
        "Daltons.",
        True,
        True,
    ),
    (
        "Cmass",
        "C mass",
        "The C-terminal mass offset from the end of the sequence tag in Daltons.",
        True,
        True,
    ),
    (
        "DeltaMass",
        "Δ mass",
        "Delta mass is the difference between the tag flanking mass and the "
        "(partial) proteoform mass, from its terminal to the tag boundary.",
        True,
        False,
    ),
]


def _tnt_column_definitions(present_fields, spec) -> List[Dict[str, Any]]:
    """Build Tabulator ``column_definitions`` from a (field,title,tooltip,float,
    dash) spec, emitting only fields present in the data.

    Numeric columns get a ``number`` sorter; float columns additionally get the
    fixed-precision ``money`` formatter (legacy ``toFixed(4)``). ``dash`` columns
    record ``_dashSentinel`` so the caller can null their -1.0 sentinel in the
    data (the legacy "-1 -> '-'" formatter); the key is ignored by Tabulator.
    """
    present = set(present_fields)
    defs: List[Dict[str, Any]] = []
    for field, title, tooltip, is_float, dash in spec:
        if field not in present:
            continue
        col: Dict[str, Any] = {
            "title": title,
            "field": field,
            "headerTooltip": tooltip,
            "sorter": "number",
        }
        if is_float:
            col.update(_FLOAT_FMT)
        if dash:
            col["_dashSentinel"] = True
        defs.append(col)
    return defs


def _apply_dash_sentinels(
    data: "pl.LazyFrame", column_defs: List[Dict[str, Any]]
) -> "pl.LazyFrame":
    """Null out the FLASHTnT -1.0 "unmatched" sentinel for the dash columns.

    The legacy tables rendered these cells with a ``-1 -> '-'`` formatter. The
    OpenMS-Insight Vue table only resolves named custom formatters (it cannot
    receive an inline JS function through the JSON-serialized column definitions),
    so we map the sentinel to null at the data layer: the cell renders blank
    instead of a misleading ``-1`` while the column stays numeric and sortable.
    """
    dash_fields = [c["field"] for c in column_defs if c.get("_dashSentinel")]
    # Strip the private marker so only valid Tabulator keys reach the frontend.
    for c in column_defs:
        c.pop("_dashSentinel", None)
    if not dash_fields:
        return data
    return data.with_columns(
        [
            pl.when(pl.col(f) == -1).then(None).otherwise(pl.col(f)).alias(f)
            for f in dash_fields
        ]
    )


def _max_score_per_scan(data: "pl.LazyFrame") -> "pl.LazyFrame":
    """Collapse the protein table to the single top-``Score`` row per ``Scan``.

    Reproduces the legacy "Best per spectrum" checkbox (default on): keep, for
    each scan, the row with the highest Score (ties resolved by first occurrence,
    matching the legacy Map insertion order); rows without a numeric Scan are kept
    as-is. Operates on the already-column-selected frame so the row identity (and
    ``index`` used for ``proteinIndex``) is preserved.
    """
    schema = data.collect_schema().names()
    if "Scan" not in schema or "Score" not in schema:
        return data
    # Window-based selection (no group_by/concat, so it survives the projection
    # pushdown the Table applies before collecting): for each Scan, keep the row
    # whose Score is the max for that Scan; ties resolved by first occurrence
    # (smallest row index, matching the legacy Map insertion order). Rows without a
    # Scan are kept verbatim (legacy pushes non-numeric Scan through). Column order
    # is preserved (with_columns/filter never reorder), so no realignment needed.
    return (
        data.with_row_index("_oi_row")
        .with_columns(
            pl.col("_oi_row")
            .filter(pl.col("Score") == pl.col("Score").max())
            .min()
            .over("Scan")
            .alias("_oi_best")
        )
        .filter(pl.col("Scan").is_null() | (pl.col("_oi_row") == pl.col("_oi_best")))
        .drop("_oi_row", "_oi_best")
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
        schema = data.collect_schema().names()
        col_defs = _tnt_column_definitions(schema, _PROTEIN_TABLE_COLUMNS)
        # index (proteinIndex) must travel through even though it has no column def;
        # plus Scan so "Best per spectrum" can collapse on it.
        keep = [c["field"] for c in col_defs] + [
            c for c in ("index", "Scan") if c in schema
        ]
        cols = list(dict.fromkeys(keep))  # de-dupe, preserve order
        data = _apply_dash_sentinels(data.select(cols), col_defs)

        def _render_protein_table():
            import streamlit as st

            # Legacy default: "Best per spectrum" is ON (bestPerSpectrumOnly: true).
            best_only = st.checkbox(
                "Best per spectrum",
                value=True,
                key=skey("protein_best_per_spectrum"),
                help="Show only the highest-scoring proteoform per spectrum (scan).",
            )
            shown = _max_score_per_scan(data) if best_only else data
            # Distinct cache_id per toggle state so the two row sets cache cleanly.
            suffix = "best" if best_only else "all"
            tbl = Table(
                cache_id=cid(f"protein_table_{suffix}"),
                data=shown,
                interactivity={PROTEIN: "index"},
                index_field="index",
                column_definitions=col_defs,
                go_to_fields=[f for f in ("Scan", "accession") if f in schema],
                initial_sort=[{"column": "Score", "dir": "desc"}],
                title="Protein Table",
                cache_path=cache_dir,
            )
            return tbl(key=skey("protein_table"), state_manager=state_manager)

        return _render_protein_table

    # ---- Tag table (filtered by proteinIndex) ----
    if comp_name == "tag_table":
        data = _load_polars(file_manager, dataset_id, "tag_dfs")
        schema = data.collect_schema().names()
        col_defs = _tnt_column_definitions(schema, _TAG_TABLE_COLUMNS)
        # TagIndex (index/interactivity) and ProteinIndex (filter) must travel
        # through even though they carry no column definition.
        keep = [c["field"] for c in col_defs] + [
            c for c in ("TagIndex", "ProteinIndex") if c in schema
        ]
        cols = list(dict.fromkeys(keep))  # de-dupe, preserve order
        data = _apply_dash_sentinels(data.select(cols), col_defs)
        has_span = {"StartPos", "EndPos"} <= set(schema)

        def _render_tag_table():
            # Residue-driven tag filter: when a sequence-view residue is selected
            # (AApos), restrict the tag table to tags spanning that residue
            # (StartPos <= AApos <= EndPos) -- the legacy residue-click tag filter.
            # It's a range predicate the identifier->column `filters` map can't
            # express, so apply it server-side at render time (AApos is read from
            # state, like the protein table's "best per spectrum" toggle) with a
            # cache_id that varies per residue (each filtered view caches cleanly).
            aa_pos = state_manager.get_selection(AAPOS)
            shown, tag_cid = data, cid("tag_table")
            if aa_pos is not None and has_span:
                shown = data.filter(
                    (pl.col("StartPos") <= int(aa_pos))
                    & (pl.col("EndPos") >= int(aa_pos))
                )
                tag_cid = f"{cid('tag_table')}_aa{int(aa_pos)}"
            tbl = Table(
                cache_id=tag_cid,
                data=shown,
                filters={PROTEIN: "ProteinIndex"},
                interactivity={TAG: "TagIndex"},
                index_field="TagIndex",
                column_definitions=col_defs,
                go_to_fields=[
                    f
                    for f in ("Scan", "StartPos", "EndPos", "TagSequence")
                    if f in schema
                ],
                initial_sort=[{"column": "Score", "dir": "desc"}],
                title="Tag Table",
                cache_path=cache_dir,
            )
            return tbl(key=skey("tag_table"), state_manager=state_manager)

        return _render_tag_table

    # ---- Combined / augmented spectrum (deconv primary + annotated overlay) ----
    if comp_name == "combined_spectrum":
        per_scan = _load_polars(file_manager, dataset_id, "combined_spectrum")
        _, anno_long = explode_combined_spectrum_long(per_scan)
        deconv_long = _deconv_signal_peaks_long(per_scan)
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
            # Tagger overlay: when a tag is selected (tagData pushed into state by
            # render_experiment_tnt), highlight matched sticks + draw per-charge
            # buttons and inter-residue amino-acid arrows over the spectrum.
            tag_overlay=True,
            signal_peaks_column="signal_peaks",
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
            # Residue click emits the peak (unused downstream) AND the residue's
            # 0-based index under AApos (POSITION_SENTINEL), driving the tag-table
            # residue filter + the gold selected-residue highlight in the overlay.
            interactivity={"peak": "peak_id", AAPOS: "<position>"},
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

    # When the selected protein changes, clear the tag + residue sub-selections so
    # a stale tag overlay / residue highlight from the previous proteoform doesn't
    # persist (the legacy reset these on protein-table click). A private state
    # marker tracks which protein the current tag/residue selection belongs to.
    if state_manager.get_selection("_tag_protein") != protein_index:
        state_manager.set_selection("_tag_protein", protein_index)
        if state_manager.get_selection(TAG) is not None:
            state_manager.set_selection(TAG, None)
        if state_manager.get_selection(AAPOS) is not None:
            state_manager.set_selection(AAPOS, None)

    # Resolve the selected tag -> tagData so the augmented spectrum's tagger
    # overlay can highlight tag-matched sticks and draw the charge buttons /
    # inter-residue amino-acid arrows. When a sequence-view residue (AApos) is
    # also selected, _tag_data marks the within-tag offset (gold highlight).
    # set_selection no-ops when unchanged.
    tag_index = state_manager.get_selection(TAG)
    aa_pos = state_manager.get_selection(AAPOS)
    tag_data = (
        _tag_data(file_manager, dataset_id, tag_index, aa_pos)
        if tag_index is not None
        else None
    )
    state_manager.set_selection("tagData", tag_data)

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
