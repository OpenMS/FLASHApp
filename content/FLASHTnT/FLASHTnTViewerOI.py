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
# Residue -> Tag-Table cross-link (legacy `selectionStore.selectedAApos`).
# Clicking a covered residue in the SequenceView sets this to the residue's
# PROTEIN-ABSOLUTE 0-based position; the Tag Table range-filters its rows to tags
# whose [StartPos, EndPos] span contains that position (StartPos <= pos <= EndPos),
# clearing on re-click (toggle). The SequenceView now renders the FULL protein
# sequence, so the residue grid index IS the protein-absolute position
# (`sequence_offset` == 0) and the emitted coordinate matches tag StartPos/EndPos
# for ALL proteoforms directly.
AA_KEY = "selectedAApos"
# Tag-span highlight on the SequenceView (legacy `selectedTag.{startPos,endPos}`).
# Published as {"start": StartPos, "end": EndPos, "nTerminal": bool} (protein-
# absolute indices) and consumed by the SequenceView `"tag_span"` interactivity
# sentinel, which reads this selection value to bracket-highlight the tag span.
TAG_SPAN_KEY = "tagSpan"


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

# Curated column definitions mirroring the LEGACY Vue tables (titles / order /
# field selection). The OI Table's ``_get_columns_to_select`` projects to ONLY the
# fields named here (plus the index / interactivity / filter columns), so any
# internal frame column not listed is hidden -- exactly the parity goal.

# TabulatorProteinTable.vue columns -> protein_dfs fields.
_PROTEIN_COLUMN_DEFINITIONS = [
    {"title": "Scan No.", "field": "Scan", "sorter": "number"},
    {"title": "Accession", "field": "accession", "sorter": "string"},
    {"title": "Description", "field": "description", "sorter": "string"},
    {"title": "Length", "field": "length", "sorter": "number"},
    # Legacy TabulatorProteinTable.vue renders the `-1` sentinel as "-" (the raw
    # value otherwise). The OI `dashNegativeOne` formatter reproduces the sentinel
    # rule; precision 4 matches the app-wide `toFixedFormatter()` default (4 dp,
    # used for every other numeric mass column, e.g. MonoMass in the mass table).
    {"title": "Mass", "field": "ProteoformMass", "sorter": "number",
     "formatter": "dashNegativeOne", "formatterParams": {"precision": 4}},
    {"title": "No. of Matched Fragments", "field": "MatchingFragments", "sorter": "number"},
    {"title": "No. of Modifications", "field": "ModCount", "sorter": "number"},
    {"title": "No. of Tags", "field": "TagCount", "sorter": "number"},
    {"title": "Score", "field": "Score", "sorter": "number"},
    # Q-Value also uses the `-1 -> "-"` sentinel rule (legacy formatter); 4 dp
    # matches the app-wide default decimal precision.
    {"title": "Q-Value (Proteoform Level)", "field": "ProteoformLevelQvalue", "sorter": "number",
     "formatter": "dashNegativeOne", "formatterParams": {"precision": 4}},
]

# TabulatorTagTable.vue columns -> tag_dfs fields.
_TAG_COLUMN_DEFINITIONS = [
    {"title": "Scan Number", "field": "Scan", "sorter": "number"},
    {"title": "Start Position", "field": "StartPos", "sorter": "number"},
    {"title": "End Position", "field": "EndPos", "sorter": "number"},
    {"title": "Sequence", "field": "TagSequence", "sorter": "string"},
    {"title": "Length", "field": "Length", "sorter": "number"},
    {"title": "Tag Score", "field": "Score", "sorter": "number"},
    # N/C mass use the legacy `-1 -> "-"` sentinel rule (TabulatorTagTable.vue
    # ~72-83); precision 4 matches the app-wide mass decimal default.
    {"title": "N mass", "field": "Nmass", "sorter": "number",
     "formatter": "dashNegativeOne", "formatterParams": {"precision": 4}},
    {"title": "C mass", "field": "Cmass", "sorter": "number",
     "formatter": "dashNegativeOne", "formatterParams": {"precision": 4}},
    {"title": "Δ mass", "field": "DeltaMass", "sorter": "number"},
]


def _filter_best_per_spectrum(protein_lf: pl.LazyFrame) -> pl.LazyFrame:
    """Collapse the protein frame to the highest-``Score`` proteoform per ``Scan``.

    Reproduces the legacy default-ON "Best per spectrum" toggle
    (TabulatorProteinTable.vue ~57-58, 116-198): keep, per ``Scan``, only the row
    with the maximum ``Score``; ties keep the first-seen row (lowest ``index``).
    Rows without a numeric ``Scan`` pass through unchanged (legacy passthrough).
    """
    # Rank within each Scan by descending Score, tie-broken by ascending index so a
    # deterministic single survivor is kept (mirrors the legacy first-seen tie rule
    # once the frame is read in index order).
    ranked = protein_lf.with_columns(
        pl.col("Score")
        .rank(method="ordinal", descending=True)
        .over("Scan")
        .alias("_score_rank")
    )
    # Keep the best row per Scan; rows with a null Scan are passed through (their
    # rank within the null group is irrelevant — keep them all, matching legacy).
    kept = ranked.filter(
        (pl.col("_score_rank") == 1) | pl.col("Scan").is_null()
    ).drop("_score_rank")
    return kept


def _build_protein_table(
    file_manager, experiment_id: str, cache_dir: str,
    best_per_spectrum: bool = True,
):
    data = _lazy(file_manager, experiment_id, "protein_dfs")
    if data is None:
        return None
    # Best-per-spectrum (default ON): pre-filter to the max-Score proteoform per
    # Scan BEFORE building the Table so the displayed rows / default-selected best
    # hit / pagination all operate on the collapsed set (legacy parity). The
    # checkbox in render_experiment_panel toggles this off to show every hit.
    if best_per_spectrum:
        data = _filter_best_per_spectrum(data)
    # The cache_id encodes the toggle so the ON / OFF frames get distinct caches
    # (the Table caches its preprocessed parquet by cache_id).
    suffix = "best" if best_per_spectrum else "all"
    # Protein table: clicking a row sets proteinIndex to the row's `index`.
    # Curated columns/titles match TabulatorProteinTable.vue (the `index`
    # interactivity column is auto-included by the Table but stays hidden).
    return Table(
        cache_id=f"protein_table_{experiment_id}_{suffix}",
        data=data,
        interactivity={PROTEIN_KEY: "index"},
        index_field="index",
        column_definitions=_PROTEIN_COLUMN_DEFINITIONS,
        initial_sort=[{"column": "Score", "dir": "desc"}],
        go_to_fields=["Scan", "accession"],
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
        # Residue -> Tag-Table cross-link (legacy StartPos<=selectedAApos<=EndPos):
        # when a covered residue is clicked in the SequenceView, AA_KEY holds its
        # protein-absolute position and the tags are narrowed to those whose span
        # contains it; cleared (no-op) when no residue is selected.
        range_filters={AA_KEY: ("StartPos", "EndPos")},
        interactivity={TAG_KEY: "TagIndex"},
        index_field="TagIndex",
        column_definitions=_TAG_COLUMN_DEFINITIONS,
        initial_sort=[{"column": "Score", "dir": "desc"}],
        go_to_fields=["Scan", "StartPos", "EndPos", "TagSequence"],
        title="Tag Table",
        cache_path=cache_dir,
    )


def _resolve_tag_masses(file_manager, experiment_id: str, state_manager) -> None:
    """Resolve the selected ``tagData`` (a ``TagIndex``) to its masses + residue
    walk and publish under ``tagMasses`` so the combined-spectrum LinePlot tagger
    overlay renders the tag walk (residue letters between consecutive masses, with
    the x-axis auto-zoomed to the tag span). Clears ``tagMasses`` when no tag is
    selected.

    Only the selected tag's row is collected (filtered by ``TagIndex``). The tag
    ``mzs`` are a comma-joined string (trailing comma); parse and drop non-numeric
    entries, keeping the STORED order (ascending for C-term tags, descending for
    N-term tags). ``TagSequence`` gives the residue letters; the legacy walks
    consecutive stored masses labelling gap ``i`` with ``sequence[len-1-i]`` —
    i.e. the REVERSED sequence aligns to the stored-order gaps regardless of
    anchoring (verified against both an ascending C-term and a descending N-term
    tag). Do NOT sort the masses: sorting breaks the alignment for descending
    (N-term) tags. The published value is a dict
    ``{"masses": [...], "residues": [...]}`` consumed by the OI LinePlot tag walk;
    when no residues are available it carries only masses (highlight-only)."""
    def _clear_all() -> None:
        state_manager.clear_selection(TAG_MASSES_KEY)
        state_manager.clear_selection(TAG_SPAN_KEY)

    tag_index = state_manager.get_selection(TAG_KEY)
    if tag_index is None:
        _clear_all()
        return

    tags = _lazy(file_manager, experiment_id, "tag_dfs")
    if tags is None:
        _clear_all()
        return

    selected = (
        tags.filter(pl.col("TagIndex") == int(tag_index))
        .select(
            pl.col("mzs")
            .str.strip_chars(",")
            .str.split(",")
            .list.eval(pl.element().cast(pl.Float64, strict=False))
            .alias("tag_masses"),
            pl.col("TagSequence").alias("tag_sequence"),
            # Anchoring + span (legacy TabulatorTagTable.vue:142-173).
            pl.col("StartPos").alias("start_pos"),
            pl.col("EndPos").alias("end_pos"),
            pl.col("Nmass").alias("n_mass"),
        )
        .collect()
    )
    if not selected.height:
        _clear_all()
        return

    raw = selected["tag_masses"][0]
    # Keep STORED order (do not sort) so the reversed-sequence walk aligns for
    # both ascending (C-term) and descending (N-term) tags.
    masses = [m for m in raw if m is not None] if raw is not None else []
    if not masses:
        _clear_all()
        return

    # Residue letter per consecutive-mass gap (len(masses) - 1 gaps): the legacy
    # labels gap i with sequence[len-1-i], i.e. reversed(sequence) over the
    # stored-order gaps. Trim to the number of gaps.
    seq = selected["tag_sequence"][0] or ""
    residues = list(reversed(str(seq)))[: max(len(masses) - 1, 0)]

    # Terminal anchoring (legacy `nTerminal = (Nmass == -1)`): an N-terminal tag is
    # one whose N-mass is the `-1` sentinel. Forwarded into the tag-walk so the
    # LinePlot honors the requested direction.
    n_mass = selected["n_mass"][0]
    n_terminal = (n_mass is not None) and (float(n_mass) == -1.0)

    state_manager.set_selection(
        TAG_MASSES_KEY,
        {"masses": list(masses), "residues": residues, "nTerminal": n_terminal},
    )

    # Tag-span highlight on the SequenceView. StartPos/EndPos are protein-absolute
    # (matching the full-protein residue grid), so they bracket the tag directly.
    start_pos = selected["start_pos"][0]
    end_pos = selected["end_pos"][0]
    if start_pos is not None and end_pos is not None:
        state_manager.set_selection(
            TAG_SPAN_KEY,
            {"start": int(start_pos), "end": int(end_pos), "nTerminal": n_terminal},
        )
    else:
        state_manager.clear_selection(TAG_SPAN_KEY)


def _normalize_mod_ranges(raw) -> list:
    """Normalize the cache ``modifications`` field into the SequenceView
    ``mod_ranges`` shape: a list of ``{start, end, mass_diff, labels}`` dicts.

    The ``modifications`` field is a ``list[struct{start,end,mass_diff,labels}]``
    (sequence_data_store.SCHEMA) carrying ambiguous/spanning modification ranges
    (DISTINCT from per-residue fixed mods). Entries missing start/end are skipped;
    indices are protein-absolute (the SequenceView renders the full protein)."""
    if raw is None:
        return []
    out = []
    for item in raw:
        if not isinstance(item, dict):
            continue
        if item.get("start") is None or item.get("end") is None:
            continue
        md = item.get("mass_diff")
        labels = item.get("labels")
        out.append({
            "start": int(item["start"]),
            "end": int(item["end"]),
            "mass_diff": float(md) if md is not None else 0.0,
            "labels": "" if labels is None else str(labels),
        })
    return out


def _precursor_mass_lookup(file_manager, experiment_id: str) -> dict:
    """``proteoform_index -> PrecursorMass`` from ``protein_dfs`` (or ``{}``).

    The observed precursor mass is not stored in ``sequence_data``; the protein
    frame carries it per proteoform (``index`` == proteoform_index). Used for the
    SequenceView mass header ``precursor_mass`` column."""
    protein = _lazy(file_manager, experiment_id, "protein_dfs")
    if protein is None:
        return {}
    schema = protein.collect_schema().names()
    if "PrecursorMass" not in schema or "index" not in schema:
        return {}
    df = protein.select(["index", "PrecursorMass"]).collect()
    return {
        int(i): float(m)
        for i, m in zip(df["index"].to_list(), df["PrecursorMass"].to_list())
        if i is not None and m is not None
    }


def _build_sequence_frame(
    file_manager, experiment_id: str
) -> Optional[pl.LazyFrame]:
    """Build the SequenceView-ready per-proteoform sequence frame.

    Source: the per-proteoform ``sequence_data`` store (keyed by
    ``proteoform_index``). It carries the FULL protein ``sequence`` list, the
    matching full-length per-residue ``coverage`` / ``maxCoverage``, the
    determined-region bounds ``proteoform_start`` / ``proteoform_end`` (with the
    ``-2`` sentinel = that terminus undetermined / open), the observed
    ``computed_mass`` and the ambiguous ``modifications`` ranges.

    We emit the FULL protein sequence (NOT pre-sliced) so the OpenMS-Insight
    SequenceView Vue renders the truncated N-/C-flanks and undetermined termini
    itself from ``proteoform_start`` / ``proteoform_end``; the full-length
    ``coverage`` stays aligned to the full sequence. Because the full protein is
    rendered, the residue grid index IS the protein-absolute 0-based position, so
    ``sequence_offset`` is always 0 (the residue-click cross-link then emits
    positions that already match tag ``StartPos`` / ``EndPos`` directly).

    Columns emitted: ``proteoform_index`` (filter key), ``sequence`` (str, FULL
    protein), ``precursor_charge`` (=1, neutral/deconvolved peaks), ``coverage``
    (full-length list[f64]), ``maxCoverage`` (f64), ``fixed_modifications``
    (list[str]), ``sequence_offset`` (=0), ``proteoform_start`` / ``proteoform_end``
    (int, sentinel ``-2`` carried through unchanged), ``computed_mass`` (f64),
    ``precursor_mass`` (f64, from protein_dfs), and ``mod_ranges``
    (list[struct{start,end,mass_diff,labels}]).

    ``sequence_data`` is loaded with ``use_polars=True`` and arrives in EITHER of
    two formats handled identically here:

    * **parquet (current ``parseTnT``):** one row per proteoform (schema in
      ``src/render/sequence_data_store.py``) returned as a polars ``LazyFrame``.
    * **pickle dict (legacy ``.pkl.gz`` example caches):** a dict keyed by the
      proteoform index; each value a dict with the same keys."""
    if not file_manager.result_exists(experiment_id, "sequence_data"):
        return None
    store = file_manager.get_results(
        experiment_id, ["sequence_data"], use_polars=True
    )["sequence_data"]

    # Normalise either format into an iterable of per-proteoform row dicts.
    if isinstance(store, pl.LazyFrame):
        rows = store.collect().iter_rows(named=True)
    elif isinstance(store, pl.DataFrame):
        rows = store.iter_rows(named=True)
    elif isinstance(store, dict):
        if not store:
            return None
        rows = (
            {"proteoform_index": pid, **(store[pid] or {})}
            for pid in sorted(store.keys())
        )
    else:
        return None

    precursor_masses = _precursor_mass_lookup(file_manager, experiment_id)

    proteoform_indices: List[int] = []
    sequences: List[str] = []
    coverages: List[list] = []
    max_coverages: List[float] = []
    fixed_mods: List[list] = []
    sequence_offsets: List[int] = []
    proteoform_starts: List[int] = []
    proteoform_ends: List[int] = []
    computed_masses: List[float] = []
    precursor_mass_col: List[float] = []
    mod_ranges_col: List[list] = []
    for entry in rows:
        pid = entry.get("proteoform_index")
        if pid is None:
            continue
        # Emit the FULL protein sequence + full-length coverage (no slicing): the
        # Vue side derives truncation / undetermined termini from
        # proteoform_start/end, so the residue grid index already equals the
        # protein-absolute position (sequence_offset = 0).
        full = list(entry.get("sequence") or [])
        cov = list(entry.get("coverage") or [])
        start = entry.get("proteoform_start")
        end = entry.get("proteoform_end")
        proteoform_indices.append(int(pid))
        sequences.append("".join(str(a) for a in full))
        coverages.append([float(c) for c in cov])
        mc = entry.get("maxCoverage")
        max_coverages.append(float(mc) if mc is not None else 0.0)
        fm = entry.get("fixed_modifications") or []
        fixed_mods.append([str(m) for m in fm])
        # Full protein rendered => residue grid index IS protein-absolute position.
        sequence_offsets.append(0)
        # Carry the determined-region bounds through UNCHANGED, including the
        # `-2` (UNDETERMINED_TERMINUS) sentinel. Absent bound => 0 / last residue
        # default on the Vue side (no truncation); we default to 0 / len-1 here.
        proteoform_starts.append(int(start) if start is not None else 0)
        proteoform_ends.append(
            int(end) if end is not None else (len(full) - 1 if full else 0)
        )
        cm = entry.get("computed_mass")
        computed_masses.append(float(cm) if cm is not None else -1.0)
        precursor_mass_col.append(float(precursor_masses.get(int(pid), 0.0)))
        mod_ranges_col.append(_normalize_mod_ranges(entry.get("modifications")))

    if not proteoform_indices:
        return None

    out = pl.DataFrame({
        "proteoform_index": proteoform_indices,
        "sequence": sequences,
        "precursor_charge": [1] * len(proteoform_indices),
        "coverage": coverages,
        "maxCoverage": max_coverages,
        "fixed_modifications": fixed_mods,
        "sequence_offset": sequence_offsets,
        "proteoform_start": proteoform_starts,
        "proteoform_end": proteoform_ends,
        "computed_mass": computed_masses,
        "precursor_mass": precursor_mass_col,
        "mod_ranges": mod_ranges_col,
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
        # Click / span sources (all routed through the interactivity mapping):
        #  - MASS_KEY: a fragment-table row click sets it to the matched peak's
        #    peak_id (combined-spectrum cross-link).
        #  - AA_KEY: a RESIDUE click sets it to the clicked residue's protein-
        #    absolute position via the "residue_position" sentinel (Tag-Table
        #    range filter). Now that the full protein is rendered, the grid index
        #    already IS the protein-absolute position (sequence_offset == 0).
        #  - TAG_SPAN_KEY: the "tag_span" sentinel does NOT set state on click;
        #    Vue READS this selection value ({start,end,nTerminal}, protein-
        #    absolute) to bracket-highlight the selected tag's span on the
        #    sequence. Published by _resolve_tag_masses.
        interactivity={
            MASS_KEY: "peak_id",
            AA_KEY: "residue_position",
            TAG_SPAN_KEY: "tag_span",
        },
        deconvolved=True,
        compute_fixed_mods=True,
        # TnT path: keep the variable/custom-mod context menu disabled (default).
        disable_variable_modifications=True,
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
        # Charge-state drill-down: per deconv-peak row, the list of constituent
        # signal-peak m/z / charge / intensity (present on combined_spectrum_long
        # from a fresh parse; empty lists for non-signal peaks).
        signal_mz_column="signal_mzs",
        signal_charge_column="signal_charges",
        signal_intensity_column="signal_intensities",
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


def build_component(
    file_manager, experiment_id: str, cache_dir: str, comp_name: str,
    best_per_spectrum: bool = True,
):
    """Instantiate the OpenMS-Insight component for a layout cell, or None."""
    builder = COMPONENT_BUILDERS.get(comp_name)
    if builder is None:
        return None
    if comp_name == "protein_table":
        return _build_protein_table(
            file_manager, experiment_id, cache_dir,
            best_per_spectrum=best_per_spectrum,
        )
    return builder(file_manager, experiment_id, cache_dir)


def _clear_proteoform_dependent_selections(state_manager) -> None:
    """Clear the per-proteoform downstream selections (mirrors legacy
    TabulatorProteinTable.vue:235-237, which resets the selected tag / tag data /
    selected AA on a proteoform change). We also clear the resolved tag masses,
    tag span and selected mass so no stale tag/peak highlight survives the switch."""
    for ident in (TAG_KEY, TAG_MASSES_KEY, AA_KEY, TAG_SPAN_KEY, MASS_KEY):
        state_manager.clear_selection(ident)


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

    # Selection clearing on proteoform change (legacy
    # TabulatorProteinTable.vue:235-237): when the selected proteinIndex differs
    # from the last-seen one for THIS panel, clear the downstream per-proteoform
    # selections (tag / tag masses / selected AA / tag span / selected mass) BEFORE
    # building components so no stale selection leaks across proteoforms.
    last_seen_key = f"{session_key}__last_protein_index"
    current_protein = state_manager.get_selection(PROTEIN_KEY)
    if st.session_state.get(last_seen_key, "__unset__") != current_protein:
        _clear_proteoform_dependent_selections(state_manager)
        st.session_state[last_seen_key] = current_protein

    # Best-per-spectrum toggle (legacy default ON). Per-panel widget key so
    # side-by-side panels toggle independently. Shown only when a protein table is
    # in the layout.
    has_protein_table = any(
        "protein_table" in row for row in layout_info_per_exp
    )
    best_per_spectrum = True
    if has_protein_table:
        best_per_spectrum = st.checkbox(
            "Best per spectrum",
            value=True,
            key=f"tnt_oi_best_per_spectrum_{panel_index}",
            help="Show only the highest-scoring proteoform per spectrum (scan).",
        )

    # Resolve the selected tag (scalar TagIndex set by the Tag Table) into its
    # list of masses BEFORE rendering so the combined-spectrum tagger overlay
    # sees the up-to-date `tagMasses` selection this rerun.
    _resolve_tag_masses(file_manager, experiment_id, state_manager)

    for row_index, row in enumerate(layout_info_per_exp):
        columns = st.columns(len(row))
        for col, (col_index, comp_name) in zip(columns, enumerate(row)):
            with col:
                component = build_component(
                    file_manager, experiment_id, cache_dir, comp_name,
                    best_per_spectrum=best_per_spectrum,
                )
                if component is None:
                    st.warning(f"No data for '{comp_name}'.")
                    continue
                key = f"tnt_oi_{panel_index}_{row_index}_{col_index}_{comp_name}"
                component(key=key, state_manager=state_manager)
