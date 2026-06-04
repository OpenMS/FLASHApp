"""FLASHApp FileManager caches -> OpenMS-Insight-ready tidy parquet.

The oracle render layer shipped *wide, list-column, index-addressed* caches (one
row per scan with array cells; selection by positional ``iloc`` /
``SignalPeaks[massIndex]``). OpenMS-Insight components want **tidy parquet with
stable value IDs** addressed by ``filters`` / ``interactivity``. This module is
the adapter: it reads the existing FileManager caches and writes derived tidy
parquet (via ``file_manager.store_data``, so the derived frames live in the same
SQLite-indexed store and gain a ``result_path`` for ``data_path=``).

It is a **pure post-process** of the ``src/parse/*`` producers' output and does
not touch them.

Public API
----------
``build_insight_caches(file_manager, dataset_id, tool, logger=None,
regenerate=False)`` reads the oracle caches for ``(dataset_id, tool)`` and writes
the tidy parquet the Insight builders (``src/render/render.py``) consume. It is
idempotent + cache-guarded: a target is skipped when its ``name_tag`` already
exists (``file_manager.result_exists``) unless ``regenerate=True``.

Stable IDs minted here (deterministic, dataset-scoped):

* ``scan_id``     -- = oracle scan-table ``index`` (already ``0..N``)
* ``mass_id``     -- per ``(scan, mass)`` global running id
* ``peak_id``     -- per exploded signal/raw peak, global running id
* ``protein_id``  -- = ``protein_df`` ``index``
* ``tag_id``      -- per tag row
* ``feature_id``  -- = ``FeatureGroupIndex``

These become the ``interactivity`` / ``filters`` columns.

See ``migration/specs/PHASE3_PLAN.md`` sections 5.1 + Appendix A for the
per-component cache -> parquet -> filters/interactivity contract.
"""

from __future__ import annotations

import polars as pl

from src.render.scan_resolution import build_proteoform_scan_map
from src.render.sequence_data_store import reconstruct_all


# Insight pushes selections down to the parquet reader; small row groups let the
# predicate pushdown skip non-matching groups for the per-scan / per-protein
# tidy frames (one logical entity may explode to many rows).
TIDY_ROW_GROUP_SIZE = 16384


# --------------------------------------------------------------------------- #
# Generic long-format / explode helpers (all polars-lazy where practical)
# --------------------------------------------------------------------------- #
def _explode_list_cols(
    df: pl.DataFrame, by: list, list_cols: list, id_name: str
) -> pl.DataFrame:
    """Explode parallel list columns to one row per element and mint a global id.

    ``by`` columns are carried (repeated per element); ``list_cols`` are exploded
    together (they must be element-aligned, which the oracle guarantees). A global
    running ``id_name`` (0..N over the whole exploded frame) is added, plus a
    per-group 0-based position ``{id_name}_in_group`` for callers that still need
    the within-scan ordinal (the oracle ``massIndex`` analogue).
    """
    keep = by + list_cols
    src = df.select(keep)
    # Drop rows whose list cell is null/empty BEFORE exploding: polars explodes an
    # empty/null list to a single null row, which would surface a phantom null
    # entry (e.g. a null mass in the Mass Table / a null peak in a spectrum) where
    # the oracle showed nothing for an empty spectrum / zero-mass scan. The
    # ``list_cols`` are element-aligned, so guarding the first is sufficient.
    primary = list_cols[0]
    src = src.filter(
        pl.col(primary).is_not_null() & (pl.col(primary).list.len() > 0)
    )
    exploded = src.explode(list_cols)
    # per-group 0-based position (replacement for the oracle positional index)
    if by:
        exploded = exploded.with_columns(
            pl.int_range(pl.len()).over(by).alias(f"{id_name}_in_group")
        )
    exploded = exploded.with_row_index(id_name)
    return exploded


def _explode_nested_signal_peaks(
    df: pl.DataFrame, scan_id_col: str, col: str, series_label: str
) -> pl.DataFrame:
    """Two-level explode of a ``SignalPeaks`` / ``NoisyPeaks`` nested cell.

    The cell is ``list[mass_idx] -> list[peak] -> [peak_index, mz, intensity,
    charge]`` (all float64; confirmed from ``masstable._compute_peak_cells`` and
    PHASE3_PLAN Appendix B). Returns one row per *point*:
    ``scan_id, mass_in_scan, peak_index, mz, intensity, charge, series`` where
    ``series`` is the supplied label ("Signal" / "Noise").

    Empty / null cells (scans with no masses, masses with no peaks) drop out, so
    the result contains only real points.
    """
    out = (
        df.select([pl.col(scan_id_col).alias("scan_id"), pl.col(col)])
        # level 1: one row per mass within a scan; position == mass_in_scan
        .explode(col)
        .with_columns(pl.int_range(pl.len()).over("scan_id").alias("mass_in_scan"))
        # drop masses whose peak list is null/empty before the inner explode
        .filter(pl.col(col).is_not_null() & (pl.col(col).list.len() > 0))
        # level 2: one row per peak record [peak_index, mz, intensity, charge]
        .explode(col)
        .filter(pl.col(col).is_not_null())
        .with_columns(
            [
                pl.col(col).list.get(0).alias("peak_index"),
                pl.col(col).list.get(1).alias("mz"),
                pl.col(col).list.get(2).alias("intensity"),
                pl.col(col).list.get(3).cast(pl.Int64).alias("charge"),
                pl.lit(series_label).alias("series"),
            ]
        )
        # neutral-ish mass the oracle 3D plot puts on its x-axis ("Mass"):
        # get3DplotInputFromSNRPeaks uses x = peaks[mz] * peaks[charge].
        .with_columns((pl.col("mz") * pl.col("charge")).alias("mass"))
        .drop(col)
    )
    return out


def _comma_split_long(df: pl.DataFrame, by: list, point_cols: dict) -> pl.DataFrame:
    """Explode comma-joined per-trace strings (FLASHQuant) to one row per point.

    ``df`` is the trace-level frame (one row per trace). ``point_cols`` maps a
    source string column (e.g. ``"MZs"``) to the output column (``"mz"``); each
    source cell is a comma-joined list of point values for that trace. All the
    point columns of one trace are element-aligned, so they explode together.
    ``by`` columns (feature_id, charge, isotope, centroid_mz) are repeated.
    """
    src = list(point_cols.keys())
    out = df.select(
        by
        + [
            pl.col(s)
            .cast(pl.Utf8)
            .str.split(",")
            .alias(point_cols[s])
            for s in src
        ]
    )
    out = out.explode([point_cols[s] for s in src])
    out = out.with_columns(
        [pl.col(point_cols[s]).cast(pl.Float64) for s in src]
    ).filter(pl.col(point_cols[src[0]]).is_not_null())
    return out


def _kde_to_long(target_df, decoy_df) -> pl.DataFrame:
    """Concat two ``{x, y}`` KDE frames into one tidy ``{x, y, group}`` frame."""
    frames = []
    for frame, label in ((target_df, "target"), (decoy_df, "decoy")):
        if frame is None:
            continue
        lf = pl.from_pandas(frame) if not isinstance(frame, pl.DataFrame) else frame
        if lf.height == 0:
            # keep schema-consistent empty contribution
            lf = pl.DataFrame({"x": [], "y": []}, schema={"x": pl.Float64, "y": pl.Float64})
        lf = lf.select(
            [pl.col("x").cast(pl.Float64), pl.col("y").cast(pl.Float64)]
        ).with_columns(pl.lit(label).alias("group"))
        frames.append(lf)
    if not frames:
        return pl.DataFrame(
            {"x": [], "y": [], "group": []},
            schema={"x": pl.Float64, "y": pl.Float64, "group": pl.Utf8},
        )
    return pl.concat(frames, how="vertical")


# --------------------------------------------------------------------------- #
# store guard
# --------------------------------------------------------------------------- #
def _store(file_manager, dataset_id, name_tag, frame, regenerate, logger=None,
           row_group_size=None):
    """Store ``frame`` under ``name_tag`` unless already present (cache guard)."""
    if (not regenerate) and file_manager.result_exists(dataset_id, name_tag):
        return False
    file_manager.store_data(dataset_id, name_tag, frame, row_group_size=row_group_size)
    if logger is not None:
        logger.log(f"[schema] wrote {name_tag} for {dataset_id}", level=2)
    return True


def _get(file_manager, dataset_id, name_tag, use_polars=False):
    """Fetch one oracle cache (pandas by default, polars LazyFrame if asked)."""
    return file_manager.get_results(
        dataset_id, [name_tag], use_polars=use_polars
    )[name_tag]


# --------------------------------------------------------------------------- #
# FLASHDeconv builders
# --------------------------------------------------------------------------- #
def _build_scans(file_manager, dataset_id, regenerate, logger):
    """(a) Scan table -> ``scans`` (already tidy; alias index -> scan_id)."""
    if (not regenerate) and file_manager.result_exists(dataset_id, "scans"):
        return
    df = _get(file_manager, dataset_id, "scan_table", use_polars=True)
    scans = df.with_columns(pl.col("index").alias("scan_id"))
    _store(file_manager, dataset_id, "scans", scans, regenerate, logger)


def _build_masses(file_manager, dataset_id, regenerate, logger):
    """(b) Mass table -> ``masses`` (explode list cells to one row per mass)."""
    if (not regenerate) and file_manager.result_exists(dataset_id, "masses"):
        return
    df = _get(file_manager, dataset_id, "mass_table", use_polars=True).collect()
    list_cols = [
        "MonoMass", "SumIntensity", "MinCharges", "MaxCharges",
        "MinIsotopes", "MaxIsotopes", "CosineScore", "SNR", "QScore",
    ]
    masses = _explode_list_cols(
        df.rename({"index": "scan_id"}), ["scan_id"], list_cols, "mass_id"
    ).rename({"mass_id_in_group": "mass_in_scan"})
    _store(file_manager, dataset_id, "masses", masses, regenerate, logger,
           row_group_size=TIDY_ROW_GROUP_SIZE)


def _build_deconv_spectrum(file_manager, dataset_id, regenerate, logger):
    """(c) Deconvolved spectrum -> ``deconv_spectrum`` (one row per peak)."""
    if (not regenerate) and file_manager.result_exists(dataset_id, "deconv_spectrum_tidy"):
        return
    df = _get(file_manager, dataset_id, "deconv_spectrum", use_polars=True).collect()
    tidy = _explode_list_cols(
        df.rename({"index": "scan_id"}),
        ["scan_id"], ["MonoMass", "SumIntensity"], "peak_id",
    ).rename({
        # SequenceView requires a peak-mass column literally named ``mass``; the
        # deconvolved monoisotopic mass IS that neutral mass.
        "MonoMass": "mass",
        # per-scan ordinal == the oracle ``massIndex`` space the 3D S/N plot and
        # the Mass Table share (onPlotClick selects the index into ``MonoMass``).
        "peak_id_in_group": "mass_in_scan",
    })
    _store(file_manager, dataset_id, "deconv_spectrum_tidy", tidy, regenerate, logger,
           row_group_size=TIDY_ROW_GROUP_SIZE)


def _build_anno_spectrum(file_manager, dataset_id, regenerate, logger):
    """(d.1) Annotated spectrum -> ``anno_spectrum`` (raw m/z, is_signal flag).

    Explode ``MonoMass_Anno`` / ``SumIntensity_Anno`` (raw m/z arrays). ``is_signal``
    marks peaks whose positional index appears in any ``SignalPeaks`` record's
    ``peak_index`` for that scan -> the LinePlot ``highlight_column``.
    """
    if (not regenerate) and file_manager.result_exists(dataset_id, "anno_spectrum_tidy"):
        return
    df = _get(file_manager, dataset_id, "combined_spectrum", use_polars=True).collect()
    df = df.rename({"index": "scan_id"})

    # set of signal peak_index values per scan, from the nested SignalPeaks cell
    sig = _explode_nested_signal_peaks(df, "scan_id", "SignalPeaks", "Signal")
    sig_idx = (
        sig.select(["scan_id", pl.col("peak_index").cast(pl.Int64)])
        .unique()
        .with_columns(pl.lit(1).alias("is_signal"))
    )

    tidy = _explode_list_cols(
        df, ["scan_id"], ["MonoMass_Anno", "SumIntensity_Anno"], "peak_id"
    ).drop("peak_id_in_group")
    # positional index within scan -> match against SignalPeaks peak_index
    tidy = (
        tidy.with_columns(
            pl.int_range(pl.len()).over("scan_id").cast(pl.Int64).alias("peak_index")
        )
        .join(sig_idx, on=["scan_id", "peak_index"], how="left")
        .with_columns(pl.col("is_signal").fill_null(0).cast(pl.Int64))
        .rename({"MonoMass_Anno": "mz", "SumIntensity_Anno": "intensity"})
        .select(["scan_id", "peak_id", "mz", "intensity", "is_signal"])
    )
    _store(file_manager, dataset_id, "anno_spectrum_tidy", tidy, regenerate, logger,
           row_group_size=TIDY_ROW_GROUP_SIZE)


def _build_anno_highlight_link(file_manager, dataset_id, regenerate, logger):
    """(d.3) Annotated-spectrum highlight linkage -> ``anno_highlight_link``.

    Selection-driven highlighting: when a deconvolved *mass* is selected, the
    annotated spectrum should highlight that mass's SIGNAL peaks (and expose each
    peak's charge). This frame is the value-based map from a deconvolved mass to
    the annotated raw peaks that are its signal peaks, keyed by the SAME
    ``peak_id`` as ``anno_spectrum_tidy`` so a viewer can ``filter`` it by the
    selected ``(scan, mass)`` and read off the ``peak_id`` set to highlight.

    Columns EXACTLY: ``scan_id, peak_id, mass_in_scan, charge`` where

    * ``peak_id``      -- the ``anno_spectrum_tidy`` peak_id of the annotated raw peak
    * ``mass_in_scan`` -- the within-scan deconvolved-mass ordinal the peak is a
      signal peak for (same ordinal space as ``masses`` / ``deconv_spectrum_tidy``
      / ``precursor_signals`` -- the outer ``SignalPeaks`` index, which the oracle
      ``combined_spectrum`` join guarantees is aligned to ``MonoMass``)
    * ``charge``       -- that signal peak's charge (``SignalPeaks`` tuple[3])

    The nested ``SignalPeaks`` cell is ``list[mass_idx] -> list[peak] ->
    [annotated_peak_index, mz, intensity, charge]``. ``annotated_peak_index``
    (tuple[0]) is the positional index of the peak within the (sorted) raw
    annotated spectrum (``MonoMass_Anno``) -- the SAME positional index
    ``_build_anno_spectrum`` matches to set ``is_signal``. We join the exploded
    signal points on ``(scan_id, that positional index)`` against the
    positionally-indexed ``anno_spectrum_tidy`` to attach the stable ``peak_id``.

    1:many: a single annotated raw peak CAN be a signal peak for MULTIPLE
    deconvolved masses (the same observed m/z explained by different charge states
    of different masses), so ``(scan_id, peak_id)`` is NOT unique here -- the frame
    is one row per ``(peak, mass)`` pair (verified against the real
    ``masstable._compute_peak_cells`` algorithm; see tests).
    """
    if (not regenerate) and file_manager.result_exists(dataset_id, "anno_highlight_link"):
        return
    # need the anno tidy frame for the stable peak_id <-> positional index map
    if not file_manager.result_exists(dataset_id, "anno_spectrum_tidy"):
        return
    df = _get(file_manager, dataset_id, "combined_spectrum", use_polars=True).collect()
    df = df.rename({"index": "scan_id"})

    # one row per signal point: scan_id, mass_in_scan, peak_index(=positional), charge
    sig = _explode_nested_signal_peaks(df, "scan_id", "SignalPeaks", "Signal")
    sig = sig.select(
        [
            "scan_id",
            "mass_in_scan",
            pl.col("peak_index").cast(pl.Int64),
            pl.col("charge").cast(pl.Int64),
        ]
    )

    # rebuild the same positional index -> peak_id map the anno tidy frame uses
    # (peak_id is assigned by exploding MonoMass_Anno per scan in scan order, so the
    # within-scan positional index is the join key against SignalPeaks' peak_index).
    anno = pl.read_parquet(file_manager.result_path(dataset_id, "anno_spectrum_tidy"))
    # peak_id is the global running explode index (monotonic in scan-then-position
    # order); sort by it so the per-scan positional index is reconstructed
    # deterministically regardless of parquet row-group read order.
    pos_map = anno.select(["scan_id", "peak_id"]).sort("peak_id").with_columns(
        pl.int_range(pl.len()).over("scan_id").cast(pl.Int64).alias("peak_index")
    )

    link = (
        sig.join(pos_map, on=["scan_id", "peak_index"], how="inner")
        .select(["scan_id", "peak_id", "mass_in_scan", "charge"])
        .sort(["scan_id", "mass_in_scan", "peak_id"])
    )
    _store(file_manager, dataset_id, "anno_highlight_link", link, regenerate, logger,
           row_group_size=TIDY_ROW_GROUP_SIZE)


def _build_combined_tagger(file_manager, dataset_id, regenerate, logger):
    """(d.2) Augmented spectrum -> ``combined_tagger`` (per-scan list columns).

    ``LinePlot.tagger`` does its own explode, so this writes one row per scan
    with the list columns it consumes:
    ``scan_id, MonoMass, SumIntensity, SignalPeaks, Mzs, MzIntensities``.
    """
    if (not regenerate) and file_manager.result_exists(dataset_id, "combined_tagger"):
        return
    df = _get(file_manager, dataset_id, "combined_spectrum", use_polars=True)
    tagger = df.select(
        [
            pl.col("index").alias("scan_id"),
            pl.col("MonoMass"),
            pl.col("SumIntensity"),
            pl.col("SignalPeaks"),
            pl.col("MonoMass_Anno").alias("Mzs"),
            pl.col("SumIntensity_Anno").alias("MzIntensities"),
        ]
    )
    _store(file_manager, dataset_id, "combined_tagger", tagger, regenerate, logger,
           row_group_size=TIDY_ROW_GROUP_SIZE)


def _build_precursor_signals(file_manager, dataset_id, regenerate, logger):
    """(e) 3D S/N plot -> ``precursor_signals`` (fully exploded Signal+Noise points)."""
    if (not regenerate) and file_manager.result_exists(dataset_id, "precursor_signals"):
        return
    df = _get(file_manager, dataset_id, "threedim_SN_plot", use_polars=True).collect()
    df = df.rename({"index": "scan_id"})
    sig = _explode_nested_signal_peaks(df, "scan_id", "SignalPeaks", "Signal")
    noi = _explode_nested_signal_peaks(df, "scan_id", "NoisyPeaks", "Noise")
    both = pl.concat([sig, noi], how="vertical").with_row_index("peak_id")
    out = both.select(
        [
            "scan_id", "mass_in_scan", "peak_id",
            "mass", "mz", "charge", "intensity", "series",
        ]
    )
    _store(file_manager, dataset_id, "precursor_signals", out, regenerate, logger,
           row_group_size=TIDY_ROW_GROUP_SIZE)


def _build_qscore_density(file_manager, dataset_id, regenerate, logger,
                          target_tag, decoy_tag, out_tag):
    """(g) Score distribution -> tidy long ``{x, y, group}``."""
    if (not regenerate) and file_manager.result_exists(dataset_id, out_tag):
        return
    if not file_manager.result_exists(dataset_id, target_tag):
        return
    target = _get(file_manager, dataset_id, target_tag)
    decoy = (
        _get(file_manager, dataset_id, decoy_tag)
        if file_manager.result_exists(dataset_id, decoy_tag)
        else None
    )
    long = _kde_to_long(target, decoy)
    _store(file_manager, dataset_id, out_tag, long, regenerate, logger)


def _build_seq_deconv(file_manager, dataset_id, regenerate, logger):
    """(j, deconv) Sequence view -> ``seq_deconv`` (one row per scan, global seq).

    The global input sequence lives in the ``('sequence','sequence')`` cache.
    SequenceView enumerates + matches fragments itself, so we only need
    ``scan_id, sequence, precursor_charge`` per scan; peaks come from the
    deconv-spectrum long frame (neutral masses).
    """
    if (not regenerate) and file_manager.result_exists(dataset_id, "seq_deconv"):
        return
    if not file_manager.result_exists("sequence", "sequence"):
        return
    seq = file_manager.get_results("sequence", "sequence")["sequence"]
    sequence = seq["input_sequence"]
    scans = _get(file_manager, dataset_id, "scan_table", use_polars=True)
    # precursor charge is not tracked per scan in the oracle deconv cache; use the
    # nearest integer of PrecursorMass/MonoMass is unavailable here, so default
    # charge 1 (neutral-mass matching is charge-agnostic for deconvolved=True).
    seq_df = scans.select(
        [
            pl.col("index").alias("scan_id"),
            pl.lit(sequence).alias("sequence"),
            pl.lit(1).cast(pl.Int64).alias("precursor_charge"),
        ]
    )
    _store(file_manager, dataset_id, "seq_deconv", seq_df, regenerate, logger,
           row_group_size=TIDY_ROW_GROUP_SIZE)


# --------------------------------------------------------------------------- #
# FLASHTnT builders
# --------------------------------------------------------------------------- #
def _build_proteins(file_manager, dataset_id, regenerate, logger):
    """(h) Protein table -> ``proteins`` (already tidy; index -> protein_id).

    Also denormalize ``scan_id`` (the proteoform's representative deconv-scan row
    index) onto each protein row. This is the value-based form of the oracle's
    ``proteoform_scan_map[proteinIndex]['deconv_index']``: a protein-row click can
    then set BOTH the ``protein`` selection and the ``scan`` selection, so all the
    scan-keyed panels (augmented spectrum, sequence-view peaks, tag table) follow
    the selected proteoform to its scan -- exactly as the oracle's render-time
    scan resolution did. Proteoforms whose scan is absent get ``scan_id = -1``.

    Also mint ``is_best_per_scan`` (1/0): the oracle ProteinTable defaults to
    "best per spectrum" = the single highest-``Score`` proteoform per ``Scan``
    (ties -> first occurrence). Exactly one row per ``Scan`` gets 1. A later step
    adds the viewer toggle + filter on this flag.
    """
    if (not regenerate) and file_manager.result_exists(dataset_id, "proteins"):
        return
    df = _get(file_manager, dataset_id, "protein_dfs")  # pandas
    scan_pd = _get(file_manager, dataset_id, "scan_table")  # pandas
    scan_map = build_proteoform_scan_map(
        df[["index", "Scan"]], scan_pd[["index", "Scan"]]
    )
    scan_to_deconv = {pid: v["deconv_index"] for pid, v in scan_map.items()}
    pdf = pl.from_pandas(df)
    proteins = pdf.with_columns(
        pl.col("index").cast(pl.Int64).alias("protein_id"),
    ).with_columns(
        pl.col("protein_id")
        .map_elements(lambda p: scan_to_deconv.get(int(p), -1), return_dtype=pl.Int64)
        .alias("scan_id"),
    ).with_columns(
        # round-8 finding 3-tables-002: the oracle ProteinTable defaults to "best
        # per spectrum" = the single highest-Score proteoform per Scan (ties ->
        # first-seen, matching the oracle's keep-first ``>`` semantics). Flag that
        # representative row 1, else 0. ``rank("ordinal", descending=True)`` gives a
        # strict 1..N ranking with NO ties, so EXACTLY one row per Scan == 1; the
        # ordinal tiebreak follows row order (first occurrence wins on equal Score).
        # A later step adds the viewer toggle + filter; we only mint the flag.
        # Null/NaN/non-numeric-Scan proteoforms are PASSED THROUGH (flagged best):
        # the oracle filterBestPerSpectrum keeps every row whose Scan is
        # `typeof !== 'number' || isNaN(scan)` rather than collapsing them into one
        # .over(Scan) group (round-10 finding 3-best-002). A missing Scan from
        # protein.tsv reads as float NaN (not a polars null), so is_null() alone
        # would miss it -- cast to f64 (non-numeric -> null) then treat null|NaN as
        # missing (dtype-safe: is_nan errors on an int column without the cast).
        (
            (pl.col("Score").rank("ordinal", descending=True).over("Scan") == 1)
            | pl.col("Scan").cast(pl.Float64, strict=False).is_nan().fill_null(True)
        )
        .cast(pl.Int64)
        .alias("is_best_per_scan"),
    )
    _store(file_manager, dataset_id, "proteins", proteins, regenerate, logger)


def _build_tags(file_manager, dataset_id, regenerate, logger):
    """(i) Tag table -> ``tags`` with a denormalized ``scan_id`` column.

    Tags are scan (spectrum) data. The oracle resolved the selected proteoform ->
    its scan via ``proteoform_scan_map`` and filtered the tag table by ``Scan``,
    so EVERY tag on that scan showed for ANY proteoform sharing the scan. We keep
    that semantics value-based: each tag carries the ``scan_id`` (deconv-row index)
    of its ``Scan``, and the builder filters ``{"scan": "scan_id"}`` -- driven by
    the protein-row click that also sets the ``scan`` selection (see
    ``_build_proteins``). We deliberately do NOT bake a per-tag ``protein_id``:
    that collapsed multi-proteoform-per-scan to one proteoform (last-wins) and hid
    the other proteoforms' tags. Tags whose scan is absent get ``scan_id = -1``.
    """
    if (not regenerate) and file_manager.result_exists(dataset_id, "tags"):
        return
    tag_pd = _get(file_manager, dataset_id, "tag_dfs")  # pandas
    protein_pd = _get(file_manager, dataset_id, "protein_dfs")  # pandas
    scan_pd = _get(file_manager, dataset_id, "scan_table")  # pandas

    # Scan number -> deconv-row index (scan_id), via the proteoform scan map.
    scan_map = build_proteoform_scan_map(
        protein_pd[["index", "Scan"]], scan_pd[["index", "Scan"]]
    )
    scan_to_deconv = {v["scan"]: v["deconv_index"] for v in scan_map.values()}

    tdf = pl.from_pandas(tag_pd).with_row_index("tag_id")
    tdf = tdf.with_columns(
        pl.col("Scan")
        .map_elements(lambda s: scan_to_deconv.get(int(s), -1)
                      if s is not None else -1, return_dtype=pl.Int64)
        .alias("scan_id"),
    )
    _store(file_manager, dataset_id, "tags", tdf, regenerate, logger,
           row_group_size=TIDY_ROW_GROUP_SIZE)


def _build_seq_tnt(file_manager, dataset_id, regenerate, logger):
    """(j, tnt) Sequence view -> ``seq_tnt`` (one row per proteoform).

    Coverage / proteoform terminals come straight from the oracle
    ``sequence_data`` store entry. SequenceView matches fragments itself from
    ``sequence`` + ``annotation_config``; the precomputed theoretical-fragment
    list-of-lists is no longer needed.
    """
    if (not regenerate) and file_manager.result_exists(dataset_id, "seq_tnt"):
        return
    if not file_manager.result_exists(dataset_id, "sequence_data"):
        return
    seq_ds = file_manager.get_results(
        dataset_id, ["sequence_data"], use_pyarrow=True
    )["sequence_data"]
    entries = reconstruct_all(seq_ds)  # {proteoform_index: entry}

    rows = []
    for pid in sorted(entries):
        e = entries[pid]
        rows.append(
            {
                "protein_id": int(pid),
                "sequence": "".join(e["sequence"]),
                "precursor_charge": 1,
                "coverage": [float(c) for c in (e.get("coverage") or [])],
                "proteoform_start": int(e.get("proteoform_start", -1)),
                "proteoform_end": int(e.get("proteoform_end", -1)),
            }
        )
    if not rows:
        return
    seq_df = pl.DataFrame(
        rows,
        schema={
            "protein_id": pl.Int64,
            "sequence": pl.Utf8,
            "precursor_charge": pl.Int64,
            "coverage": pl.List(pl.Float64),
            "proteoform_start": pl.Int64,
            "proteoform_end": pl.Int64,
        },
    )
    _store(file_manager, dataset_id, "seq_tnt", seq_df, regenerate, logger,
           row_group_size=TIDY_ROW_GROUP_SIZE)


# --------------------------------------------------------------------------- #
# FLASHQuant builders
# --------------------------------------------------------------------------- #
_QUANT_SCALAR_RENAME = {
    "FeatureGroupIndex": "feature_id",
    "StartRetentionTime(FWHM)": "StartRT",
    "EndRetentionTime(FWHM)": "EndRT",
    "HighestApexRetentionTime": "ApexRT",
    "AllAreaUnderTheCurve": "AllAUC",
}


def _build_quant(file_manager, dataset_id, regenerate, logger):
    """(k) FLASHQuant -> ``quant_features`` (tidy scalars) + ``quant_traces`` (long).

    The oracle ``quant_dfs`` is one row per FeatureGroup with scalar columns plus
    list columns (``Charges/IsotopeIndices/CentroidMzs``) and comma-joined
    per-trace strings (``RTs/MZs/Intensities``). We split into:

    * ``quant_features`` -- one row per feature (scalars), ``feature_id`` minted.
    * ``quant_traces``   -- one row per trace *point* (comma-split + explode);
      each point carries ``trace_in_feature``, a stable per-feature running id of
      its parent trace so the 3D can break the polyline per-trace (the oracle's
      -1000 z sentinel) even when two traces share ``(charge, isotope)``.
    """
    need_feat = regenerate or not file_manager.result_exists(dataset_id, "quant_features")
    need_traces = regenerate or not file_manager.result_exists(dataset_id, "quant_traces")
    if not (need_feat or need_traces):
        return
    df = _get(file_manager, dataset_id, "quant_dfs")  # pandas
    pdf = pl.from_pandas(df)

    # ---- feature scalars ----
    if need_feat:
        scalar_cols = [c for c in pdf.columns if c not in (
            "Charges", "IsotopeIndices", "CentroidMzs", "RTs", "MZs", "Intensities",
        )]
        feats = pdf.select(scalar_cols).rename(
            {k: v for k, v in _QUANT_SCALAR_RENAME.items() if k in scalar_cols}
        )
        feats = feats.with_columns(pl.col("feature_id").cast(pl.Int64))
        _store(file_manager, dataset_id, "quant_features", feats, regenerate, logger)

    # ---- trace points (one row per trace, then comma-split to one row per point) ----
    if need_traces:
        # explode the per-trace list columns (Charges/IsotopeIndices/CentroidMzs and
        # the comma-joined RTs/MZs/Intensities strings move together, one per trace)
        trace_lists = ["Charges", "IsotopeIndices", "CentroidMzs", "RTs", "MZs", "Intensities"]
        per_trace = (
            pdf.select(
                [pl.col("FeatureGroupIndex").cast(pl.Int64).alias("feature_id")]
                + [pl.col(c) for c in trace_lists]
            )
            .explode(trace_lists)
            # Stable per-feature running trace id (round-8 finding 3-quant-005): the
            # 3D wraps EVERY trace in a -1000 z sentinel, so the polyline must break
            # per-trace. (charge, isotope) is NOT unique -- two traces of one feature
            # can share it -- so mint a distinct id per exploded trace row and carry
            # it through to every point so a trace can be drawn as one isolated line.
            .with_columns(
                pl.int_range(pl.len()).over("feature_id").alias("trace_in_feature")
            )
            .rename(
                {
                    "Charges": "charge",
                    "IsotopeIndices": "isotope",
                    "CentroidMzs": "centroid_mz",
                }
            )
            .with_columns(
                [
                    pl.col("charge").cast(pl.Int64),
                    pl.col("isotope").cast(pl.Int64),
                    pl.col("centroid_mz").cast(pl.Float64),
                ]
            )
        )
        traces = _comma_split_long(
            per_trace,
            ["feature_id", "charge", "isotope", "centroid_mz", "trace_in_feature"],
            {"RTs": "rt", "MZs": "mz", "Intensities": "intensity"},
        )
        _store(file_manager, dataset_id, "quant_traces", traces, regenerate, logger,
               row_group_size=TIDY_ROW_GROUP_SIZE)


# --------------------------------------------------------------------------- #
# public entry point
# --------------------------------------------------------------------------- #
def build_insight_caches(file_manager, dataset_id, tool, logger=None,
                         regenerate=False) -> None:
    """Read the oracle caches for ``(dataset_id, tool)`` and write Insight tidy parquet.

    Idempotent + cache-guarded: a target is skipped when its ``name_tag`` already
    exists unless ``regenerate=True``. ``tool`` selects the panel set:

    * ``"flashdeconv"`` -- scans, masses, deconv/anno/tagger spectra, the
      annotated-spectrum highlight linkage, 3D S/N, qscore density, (optional)
      global sequence view. Heatmaps reuse the existing full-resolution
      ``ms*_{deconv,raw}_heatmap`` caches as-is.
    * ``"flashtnt"`` -- everything deconv has, plus proteins, tags, per-proteoform
      sequence view, and the id-FDR density.
    * ``"flashquant"`` -- quant feature scalars + exploded trace points.
    """
    tool = (tool or "").lower()

    if tool == "flashquant":
        _build_quant(file_manager, dataset_id, regenerate, logger)
        return

    # ---- shared deconv-style panels (flashdeconv + flashtnt) ----
    _build_scans(file_manager, dataset_id, regenerate, logger)
    _build_masses(file_manager, dataset_id, regenerate, logger)
    _build_deconv_spectrum(file_manager, dataset_id, regenerate, logger)
    _build_anno_spectrum(file_manager, dataset_id, regenerate, logger)
    _build_anno_highlight_link(file_manager, dataset_id, regenerate, logger)
    _build_combined_tagger(file_manager, dataset_id, regenerate, logger)
    _build_precursor_signals(file_manager, dataset_id, regenerate, logger)

    if tool == "flashdeconv":
        _build_qscore_density(
            file_manager, dataset_id, regenerate, logger,
            "density_target", "density_decoy", "qscore_density",
        )
        _build_seq_deconv(file_manager, dataset_id, regenerate, logger)
    elif tool == "flashtnt":
        _build_qscore_density(
            file_manager, dataset_id, regenerate, logger,
            "density_id_target", "density_id_decoy", "qscore_density_id",
        )
        _build_proteins(file_manager, dataset_id, regenerate, logger)
        _build_tags(file_manager, dataset_id, regenerate, logger)
        _build_seq_tnt(file_manager, dataset_id, regenerate, logger)
