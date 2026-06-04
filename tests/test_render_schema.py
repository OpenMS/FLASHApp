"""Construct-smoke for ``src.render.schema.build_insight_caches``.

Builds synthetic FLASHApp FileManager caches (matching the ``src/parse/*`` output
schemas), runs ``build_insight_caches`` for each tool, and asserts the Insight-ready
tidy parquet is produced with the stable-ID columns and the right explode shapes.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

import pandas as pd

from src.workflow.FileManager import FileManager
from src.render.schema import (
    build_insight_caches,
    _explode_list_cols,
    _explode_nested_signal_peaks,
    _comma_split_long,
    _kde_to_long,
    _build_proteins,
    _build_seq_tnt,
)
from tests.conftest import make_deconv_caches, make_tnt_caches, make_quant_caches, \
    make_sequence_cache


def _fm(workspace):
    return FileManager(workspace, Path(workspace, "cache"))


# --------------------------------------------------------------------------- #
# helper-level unit checks (the explode/comma-split/kde primitives)
# --------------------------------------------------------------------------- #
def test_explode_list_cols_mints_global_and_group_ids():
    df = pl.DataFrame({"scan_id": [0, 1], "MonoMass": [[100.0, 200.0], [300.0]],
                       "SumIntensity": [[1.0, 2.0], [3.0]]})
    out = _explode_list_cols(df, ["scan_id"], ["MonoMass", "SumIntensity"], "peak_id")
    assert out.height == 3
    assert out["peak_id"].to_list() == [0, 1, 2]
    assert out["peak_id_in_group"].to_list() == [0, 1, 0]  # per-scan ordinal


def test_explode_list_cols_drops_empty_and_null_cells():
    # a scan with an empty mass list (zero-mass scan) and one with null must NOT
    # surface a phantom null row (the oracle showed nothing for an empty spectrum).
    df = pl.DataFrame(
        {"scan_id": [0, 1, 2], "MonoMass": [[100.0, 200.0], [], None],
         "SumIntensity": [[1.0, 2.0], [], None]},
        schema={"scan_id": pl.Int64, "MonoMass": pl.List(pl.Float64),
                "SumIntensity": pl.List(pl.Float64)})
    out = _explode_list_cols(df, ["scan_id"], ["MonoMass", "SumIntensity"], "peak_id")
    assert out.height == 2  # only scan 0's two real masses
    assert out["scan_id"].to_list() == [0, 0]
    assert out["MonoMass"].null_count() == 0


def test_explode_nested_signal_peaks_two_levels():
    sp = pl.DataFrame(
        {"scan_id": [0, 1],
         "SignalPeaks": [
             [[[0.0, 75.0, 3.0, 12.0], [1.0, 75.1, 1.0, 12.0]], [[3.0, 125.0, 4.0, 5.0]]],
             [[[5.0, 100.0, 1.0, 1.0]]]]},
        schema={"scan_id": pl.Int64, "SignalPeaks": pl.List(pl.List(pl.List(pl.Float64)))})
    out = _explode_nested_signal_peaks(sp, "scan_id", "SignalPeaks", "Signal")
    assert out.height == 4
    assert out["mass_in_scan"].to_list() == [0, 0, 1, 0]
    assert out["charge"].to_list() == [12, 12, 5, 1]
    assert set(out["series"].unique().to_list()) == {"Signal"}
    # oracle 3D x = mz * charge
    assert out["mass"].to_list() == [75.0 * 12, 75.1 * 12, 125.0 * 5, 100.0 * 1]


def test_explode_nested_handles_empty_cells():
    sp = pl.DataFrame(
        {"scan_id": [0], "SignalPeaks": [[[]]]},
        schema={"scan_id": pl.Int64, "SignalPeaks": pl.List(pl.List(pl.List(pl.Float64)))})
    out = _explode_nested_signal_peaks(sp, "scan_id", "SignalPeaks", "Noise")
    assert out.height == 0


def test_comma_split_long_explodes_points():
    tr = pl.DataFrame({"feature_id": [0], "charge": [2], "isotope": [0],
                       "centroid_mz": [500.0], "RTs": ["1.0,2.0,3.0"],
                       "MZs": ["500.1,500.2,500.3"], "Intensities": ["10,20,30"]})
    out = _comma_split_long(tr, ["feature_id", "charge", "isotope", "centroid_mz"],
                            {"RTs": "rt", "MZs": "mz", "Intensities": "intensity"})
    assert out.height == 3
    assert out["rt"].to_list() == [1.0, 2.0, 3.0]
    assert out["intensity"].to_list() == [10.0, 20.0, 30.0]


def test_kde_to_long_concats_with_group_and_handles_missing_decoy():
    import pandas as pd
    t = pd.DataFrame({"x": [0.1, 0.2], "y": [1.0, 2.0]})
    d = pd.DataFrame({"x": [0.3], "y": [0.5]})
    out = _kde_to_long(t, d)
    assert out.height == 3
    assert set(out["group"].unique().to_list()) == {"target", "decoy"}
    # decoy absent -> only target rows
    assert set(_kde_to_long(t, None)["group"].unique().to_list()) == {"target"}


# --------------------------------------------------------------------------- #
# FLASHDeconv tidy parquet
# --------------------------------------------------------------------------- #
def test_build_insight_caches_flashdeconv(temp_workspace):
    fm = _fm(temp_workspace)
    ds = make_deconv_caches(fm)
    make_sequence_cache(fm)

    build_insight_caches(fm, ds, "flashdeconv")

    expected = ["scans", "masses", "deconv_spectrum_tidy", "anno_spectrum_tidy",
                "combined_tagger", "precursor_signals", "qscore_density", "seq_deconv"]
    for tag in expected:
        assert fm.result_exists(ds, tag), f"missing tidy cache: {tag}"

    masses = pl.read_parquet(fm.result_path(ds, "masses"))
    assert {"scan_id", "mass_id", "mass_in_scan"}.issubset(masses.columns)
    assert masses["mass_id"].n_unique() == masses.height  # stable unique id
    assert masses.height == 3  # 2 + 1 masses exploded

    # deconv spectrum: SequenceView needs a 'mass' column; the per-scan ordinal is
    # exposed as 'mass_in_scan' (the oracle massIndex space shared with the 3D).
    deconv = pl.read_parquet(fm.result_path(ds, "deconv_spectrum_tidy"))
    assert {"scan_id", "peak_id", "mass", "mass_in_scan"}.issubset(deconv.columns)
    assert deconv.filter(pl.col("scan_id") == 0)["mass_in_scan"].to_list() == [0, 1]

    ps = pl.read_parquet(fm.result_path(ds, "precursor_signals"))
    assert {"scan_id", "mass_in_scan", "peak_id", "mass", "mz", "charge",
            "intensity", "series"}.issubset(ps.columns)
    assert ps["peak_id"].n_unique() == ps.height
    assert set(ps["series"].unique().to_list()) <= {"Signal", "Noise"}
    # 3D x-axis is the oracle "Mass" = mz * charge (get3DplotInputFromSNRPeaks),
    # not raw m/z.
    assert ps.select(
        (pl.col("mass") - pl.col("mz") * pl.col("charge")).abs().max()
    ).item() < 1e-9

    anno = pl.read_parquet(fm.result_path(ds, "anno_spectrum_tidy"))
    assert {"scan_id", "peak_id", "mz", "intensity", "is_signal"}.issubset(anno.columns)
    # is_signal = membership in SignalPeaks.peak_index. scan 0 SignalPeaks cover
    # peak_index {0,1,3} (3 signal); scan 1 covers {0} (1 signal) -> 4 total.
    assert int(anno["is_signal"].sum()) == 4

    seq = pl.read_parquet(fm.result_path(ds, "seq_deconv"))
    assert {"scan_id", "sequence", "precursor_charge"}.issubset(seq.columns)
    assert seq["sequence"].unique().to_list() == ["PEPTIDEK"]


def test_anno_highlight_link(temp_workspace):
    """Annotated-spectrum highlight linkage: maps each annotated SIGNAL peak to the
    deconvolved mass (mass_in_scan) + charge it is a signal peak for, keyed by the
    SAME peak_id as anno_spectrum_tidy. Verifies columns, peak_id consistency,
    a known peak's (mass_in_scan, charge), and the 1:many relationship."""
    fm = _fm(temp_workspace)
    ds = make_deconv_caches(fm)
    build_insight_caches(fm, ds, "flashdeconv")

    assert fm.result_exists(ds, "anno_highlight_link")
    link = pl.read_parquet(fm.result_path(ds, "anno_highlight_link"))
    anno = pl.read_parquet(fm.result_path(ds, "anno_spectrum_tidy"))

    # EXACT columns.
    assert link.columns == ["scan_id", "peak_id", "mass_in_scan", "charge"]

    # anno_spectrum_tidy peak_id is a stable per-row id (unique within the frame).
    assert anno["peak_id"].n_unique() == anno.height

    # every link peak_id exists in anno_spectrum_tidy with the same scan_id (the
    # linkage is keyed by the anno peak_id), and only SIGNAL peaks are linked.
    joined = link.join(anno, on=["scan_id", "peak_id"], how="left")
    assert joined["mz"].null_count() == 0           # all link peak_ids resolve
    assert int(joined["is_signal"].min()) == 1      # linked peaks are all signal

    # Known signal peak. Synthetic combined_spectrum scan 0:
    #   anno peaks (sorted): idx0 m/z 75.0, idx1 75.1, idx2 125.0, idx3 99.0
    #   SignalPeaks mass0 -> anno idx0(z12), idx1(z12); mass1 -> idx3(z5), idx0(z5)
    # so anno idx0 (peak_id of scan0/pos0) links to mass_in_scan 0 (z12) AND 1 (z5).
    pid0 = anno.filter(pl.col("scan_id") == 0).sort("peak_id")["peak_id"].to_list()[0]
    rows0 = link.filter((pl.col("scan_id") == 0) & (pl.col("peak_id") == pid0)).sort(
        "mass_in_scan"
    )
    assert rows0["mass_in_scan"].to_list() == [0, 1]
    assert rows0["charge"].to_list() == [12, 5]

    # CRITICAL 1:1 vs 1:many finding: a single annotated raw peak CAN belong to
    # MULTIPLE deconvolved masses, so (scan_id, peak_id) is NOT unique -> the frame
    # is 1:many (one row per (peak, mass) pair). Assert the dup pair is present.
    dup = (
        link.group_by(["scan_id", "peak_id"]).len().filter(pl.col("len") > 1)
    )
    assert dup.height >= 1, "expected at least one annotated peak mapping to >1 mass"
    # and that the link allows it (the dup we constructed: scan 0, peak pos 0).
    assert (
        link.filter((pl.col("scan_id") == 0) & (pl.col("peak_id") == pid0)).height == 2
    )


def test_build_insight_caches_idempotent(temp_workspace):
    fm = _fm(temp_workspace)
    ds = make_deconv_caches(fm)
    build_insight_caches(fm, ds, "flashdeconv")
    before = pl.read_parquet(fm.result_path(ds, "masses")).height
    # second call must not error and must leave the cache untouched (guarded)
    build_insight_caches(fm, ds, "flashdeconv")
    after = pl.read_parquet(fm.result_path(ds, "masses")).height
    assert before == after


# --------------------------------------------------------------------------- #
# FLASHTnT tidy parquet
# --------------------------------------------------------------------------- #
def test_build_insight_caches_flashtnt(temp_workspace):
    fm = _fm(temp_workspace)
    ds = make_tnt_caches(fm)

    build_insight_caches(fm, ds, "flashtnt")

    for tag in ["proteins", "tags", "seq_tnt", "qscore_density_id"]:
        assert fm.result_exists(ds, tag), f"missing tidy cache: {tag}"

    proteins = pl.read_parquet(fm.result_path(ds, "proteins"))
    assert {"protein_id", "scan_id", "is_best_per_scan"}.issubset(proteins.columns)
    assert proteins["protein_id"].to_list() == [0, 1]
    # protein row carries its scan (deconv-row index): Scan 10 -> 0, Scan 20 -> 1,
    # so a protein-row click can resolve protein -> scan (value-based scan map).
    assert proteins["scan_id"].to_list() == [0, 1]
    # round-8 finding 3-tables-002: exactly one is_best_per_scan==1 per Scan. Here
    # each Scan (10, 20) has a single proteoform, so both rows are best.
    assert proteins["is_best_per_scan"].to_list() == [1, 1]
    best_per_scan = proteins.filter(pl.col("is_best_per_scan") == 1)
    assert best_per_scan["Scan"].n_unique() == best_per_scan.height

    tags = pl.read_parquet(fm.result_path(ds, "tags"))
    # tags are scan-keyed (NOT collapsed to a per-scan protein_id): each tag carries
    # the deconv-row index of its Scan, and the tag table follows protein->scan.
    assert {"tag_id", "scan_id"}.issubset(tags.columns)
    assert "protein_id" not in tags.columns
    # Scan 10 -> scan_id 0, Scan 20 -> scan_id 1 (from scan_table index)
    m = {r["Scan"]: r["scan_id"] for r in tags.select(["Scan", "scan_id"]).to_dicts()}
    assert m == {10: 0, 20: 1}

    seqt = pl.read_parquet(fm.result_path(ds, "seq_tnt"))
    assert {"protein_id", "sequence", "coverage", "proteoform_start",
            "proteoform_end"}.issubset(seqt.columns)
    assert sorted(seqt["sequence"].to_list()) == ["ACDEFGHK", "PEPTIDEK"]


def _make_truncated_proteoform_seq_cache(fm, ds="exp_pf"):
    """Write a ``sequence_data`` cache mirroring the oracle ``parseTnT`` output for
    a TRUNCATED proteoform (round-17 3-seqview-009).

    Full protein ``MKPEPTIDEK``; the determined proteoform is ``PEPTIDEK``
    (1-based StartPosition 3, EndPosition 10). The oracle stores the FULL protein
    in ``sequence`` but computes the fragment grid on the SLICED sub-sequence
    ``str(sequence)[start_index:end_index+1]`` and stores 0-based
    ``proteoform_start``/``proteoform_end`` (StartPosition-1 / EndPosition-1).
    """
    import numpy as np
    import pyarrow.parquet as pq
    from src.render.sequence import getFragmentDataFromSeq
    from src.render.sequence_data_store import build_table, ROW_GROUP_SIZE

    full = "MKPEPTIDEK"
    # Oracle src/parse/tnt.py slice derivation for StartPosition=3, EndPosition=10.
    start_position, end_position = 3, 10
    start_index = 0 if start_position <= 0 else start_position - 1
    end_index = len(full) - 1 if end_position <= 0 else end_position - 1
    cov = np.array([1.0] * len(full))
    # Oracle: getFragmentDataFromSeq on the SLICED sub-sequence.
    entry = getFragmentDataFromSeq(
        full[start_index:end_index + 1], list(cov / cov.max()), cov.max(), []
    )
    entry["sequence"] = list(full)                 # FULL protein in the grid
    entry["proteoform_start"] = start_position - 1  # 0-based -> 2
    entry["proteoform_end"] = end_position - 1      # 0-based -> 9
    entry["computed_mass"] = 900.0
    entry["theoretical_mass"] = 1100.0
    entry["modifications"] = []
    tbl = build_table({0: entry})
    with fm.parquet_sink(ds, "sequence_data") as p:
        pq.write_table(tbl, p, row_group_size=ROW_GROUP_SIZE)
    return ds, full, start_index, end_index


def test_seq_tnt_truncated_proteoform_carries_full_seq_and_terminals(temp_workspace):
    """``seq_tnt`` keeps the FULL protein + the 0-based proteoform terminals.

    The migrated ``_build_seq_tnt`` must surface the FULL ``sequence`` (the display
    grid) plus the reported 0-based ``proteoform_start``/``proteoform_end`` so the
    Insight SequenceView can slice the fragment grid + offset the mapping
    (3-seqview-009). It must NOT slice the stored ``sequence`` itself.
    """
    fm = _fm(temp_workspace)
    ds, full, _, _ = _make_truncated_proteoform_seq_cache(fm)

    # _build_seq_tnt only consumes the sequence_data cache; call it directly so we
    # do not need the full deconv-style cache set for this proteoform-region check.
    _build_seq_tnt(fm, ds, regenerate=True, logger=None)
    seqt = pl.read_parquet(fm.result_path(ds, "seq_tnt"))

    row = seqt.filter(pl.col("protein_id") == 0).to_dicts()[0]
    assert row["sequence"] == full          # full protein, NOT the sub-region
    assert row["proteoform_start"] == 2     # StartPosition(3) - 1
    assert row["proteoform_end"] == 9       # EndPosition(10) - 1


def test_seq_tnt_truncated_proteoform_sequenceview_matches_oracle(temp_workspace):
    """End-to-end: the SequenceView wired from ``seq_tnt`` computes the fragment
    grid on the PROTEOFORM SUB-region, numerically matching the oracle.

    Reproduces the oracle FLASHApp ``getFragmentDataFromSeq`` on the SLICED
    sub-sequence (3-seqview-009): the migrated Insight SequenceView slices
    ``sequence[proteoform_start..proteoform_end]`` and the resulting grid +
    offset match the oracle exactly (b1 == 97.05 for PEPTIDEK, not 131.04 for the
    full MKPEPTIDEK).
    """
    from openms_insight.components.sequenceview import (
        SequenceView,
        calculate_fragment_masses_pyopenms,
    )

    fm = _fm(temp_workspace)
    ds, full, start_index, end_index = _make_truncated_proteoform_seq_cache(fm)
    _build_seq_tnt(fm, ds, regenerate=True, logger=None)

    # Wire the SequenceView exactly as src/render/render.py does for flashtnt
    # (proteoform terminal columns configured).
    sv = SequenceView(
        cache_id="pf_e2e",
        sequence_data_path=fm.result_path(ds, "seq_tnt"),
        cache_path=str(Path(temp_workspace, "insight_cache")),
        filters={"protein": "protein_id"},
        proteoform_start_column="proteoform_start",
        proteoform_end_column="proteoform_end",
        deconvolved=True,
    )
    seq = sv._prepare_vue_data({"protein": 0})["sequenceData"]

    # Grid shows the full protein; fragments are on the sub-region with the offset.
    assert len(seq["sequence"]) == len(full)
    assert seq["proteoform_fragments"] is True
    assert seq["fragment_grid_offset"] == start_index  # 2

    # Numerically identical to the oracle sub-region grid.
    sub = full[start_index:end_index + 1]
    assert sub == "PEPTIDEK"
    oracle_sub_grid = calculate_fragment_masses_pyopenms(sub)
    for ion in ("a", "b", "c", "x", "y", "z"):
        assert seq[f"fragment_masses_{ion}"] == oracle_sub_grid[f"fragment_masses_{ion}"]
    # The finding's concrete example: b1 of the proteoform region.
    assert seq["fragment_masses_b"][0][0] == __import__("pytest").approx(
        97.0527642233, abs=1e-6
    )


def test_proteins_is_best_per_scan(temp_workspace):
    """round-8 finding 3-tables-002: is_best_per_scan == 1 for the single
    highest-Score proteoform per Scan, with ties broken by first occurrence
    (oracle keep-first ``>``). Build a cache directly with a multi-proteoform Scan
    AND a Score tie so exactly one row per Scan is flagged."""
    fm = _fm(temp_workspace)
    ds = "exp1"

    # Two scans. Scan 10 has THREE proteoforms incl. a Score tie (5.0 == 5.0);
    # Scan 20 has two. The deconv scan_table maps Scan 10 -> deconv 0, Scan 20 -> 1.
    fm.store_data(ds, "scan_table", pd.DataFrame({
        "index": [0, 1], "Scan": [10, 20]}))
    fm.store_data(ds, "protein_dfs", pd.DataFrame({
        "index":  [0,    1,    2,    3,    4],
        "Scan":   [10,   10,   10,   20,   20],
        # Scan 10: max is 7.0 (proteoform 1). The 5.0 tie (0 and 2) must NOT both win.
        "Score":  [5.0,  7.0,  5.0,  3.0,  9.0],
        "accession": ["a", "b", "c", "d", "e"]}))

    _build_proteins(fm, ds, regenerate=True, logger=None)
    proteins = pl.read_parquet(fm.result_path(ds, "proteins")).sort("protein_id")

    assert "is_best_per_scan" in proteins.columns
    # exactly one best per Scan
    best = proteins.filter(pl.col("is_best_per_scan") == 1)
    assert best.height == proteins["Scan"].n_unique() == 2
    assert best["Scan"].n_unique() == best.height  # one per scan, no dup
    # the right rows: Scan 10 -> proteoform 1 (Score 7.0), Scan 20 -> proteoform 4 (9.0)
    assert set(best["protein_id"].to_list()) == {1, 4}
    # the 5.0 tie on Scan 10 produced exactly ZERO winners (max was 7.0), and even a
    # tie AT the max is broken keep-first (ordinal rank): verify per-scan sum == 1.
    by_scan = proteins.group_by("Scan").agg(pl.col("is_best_per_scan").sum())
    assert sorted(r["is_best_per_scan"] for r in by_scan.to_dicts()) == [1, 1]


def test_proteins_is_best_per_scan_tie_keeps_first(temp_workspace):
    """A Score tie AT the per-Scan maximum is broken keep-first (oracle ``>``):
    the FIRST-occurring max-Score row wins, not the later one."""
    fm = _fm(temp_workspace)
    ds = "exp1"
    fm.store_data(ds, "scan_table", pd.DataFrame({"index": [0], "Scan": [10]}))
    # both proteoforms on Scan 10 tie at the max Score 8.0; first (index 0) wins.
    fm.store_data(ds, "protein_dfs", pd.DataFrame({
        "index": [0, 1], "Scan": [10, 10], "Score": [8.0, 8.0],
        "accession": ["first", "second"]}))
    _build_proteins(fm, ds, regenerate=True, logger=None)
    proteins = pl.read_parquet(fm.result_path(ds, "proteins")).sort("protein_id")
    assert proteins["is_best_per_scan"].to_list() == [1, 0]


def test_proteins_is_best_per_scan_passthrough_missing_scan(temp_workspace):
    """round-10 finding 3-best-002: proteoforms with a MISSING Scan (NaN/null) are
    PASSED THROUGH (every one flagged best), matching the oracle filterBestPerSpectrum
    which keeps each row whose Scan is non-numeric/NaN -- NOT collapsed into one
    .over(Scan) group (which would flag only one). A missing Scan from protein.tsv
    arrives as float NaN, so the flag must catch NaN, not just polars null."""
    fm = _fm(temp_workspace)
    ds = "exp1"
    fm.store_data(ds, "scan_table", pd.DataFrame({"index": [0], "Scan": [10]}))
    # Scan 10 (one proteoform) + THREE missing-Scan (NaN) proteoforms.
    fm.store_data(ds, "protein_dfs", pd.DataFrame({
        "index": [0, 1, 2, 3],
        "Scan": [10, None, None, None],  # -> float64 [10.0, NaN, NaN, NaN]
        "Score": [5.0, 1.0, 9.0, 3.0],
        "accession": ["a", "b", "c", "d"]}))
    _build_proteins(fm, ds, regenerate=True, logger=None)
    proteins = pl.read_parquet(fm.result_path(ds, "proteins")).sort("protein_id")
    # Scan 10 -> best (alone); ALL three missing-Scan rows -> best (passthrough).
    assert proteins["is_best_per_scan"].to_list() == [1, 1, 1, 1]


def test_proteins_is_best_per_scan_nan_score_loses(temp_workspace):
    """round-11 finding 3-best-003: a NaN/missing Score must NOT win best-per-spectrum
    (oracle toScore maps NaN/non-numeric -> -Infinity, sorting it last). On a Scan
    with one real Score and one missing (NaN) Score, the REAL-Score proteoform is
    flagged best -- NOT the NaN one (which polars rank(descending) would otherwise
    rank largest)."""
    fm = _fm(temp_workspace)
    ds = "exp1"
    fm.store_data(ds, "scan_table", pd.DataFrame({"index": [0], "Scan": [10]}))
    fm.store_data(ds, "protein_dfs", pd.DataFrame({
        "index": [0, 1],
        "Scan": [10, 10],
        "Score": [5.0, None],  # -> float64 [5.0, NaN]
        "accession": ["real", "noscore"]}))
    _build_proteins(fm, ds, regenerate=True, logger=None)
    proteins = pl.read_parquet(fm.result_path(ds, "proteins")).sort("protein_id")
    # the real-Score row (5.0) wins; the NaN-Score row does NOT.
    assert proteins["is_best_per_scan"].to_list() == [1, 0]


# --------------------------------------------------------------------------- #
# FLASHQuant tidy parquet
# --------------------------------------------------------------------------- #
def test_build_insight_caches_flashquant(temp_workspace):
    fm = _fm(temp_workspace)
    ds = make_quant_caches(fm)

    build_insight_caches(fm, ds, "flashquant")

    for tag in ["quant_features", "quant_traces"]:
        assert fm.result_exists(ds, tag), f"missing tidy cache: {tag}"

    feats = pl.read_parquet(fm.result_path(ds, "quant_features"))
    assert "feature_id" in feats.columns
    assert {"StartRT", "EndRT", "ApexRT", "AllAUC"}.issubset(feats.columns)
    assert feats["feature_id"].to_list() == [0, 1, 12]

    traces = pl.read_parquet(fm.result_path(ds, "quant_traces"))
    assert {"feature_id", "charge", "isotope", "centroid_mz", "rt", "mz",
            "intensity", "trace_in_feature"}.issubset(traces.columns)
    # feature 0: 3+2 points, feature 1: 2 points, feature 12: 2+3 points
    per = {r["feature_id"]: r["len"]
           for r in traces.group_by("feature_id").len().to_dicts()}
    assert per == {0: 5, 1: 2, 12: 5}

    # round-8 finding 3-quant-005: trace_in_feature is a stable per-feature running
    # trace id, distinct PER TRACE -- even for two traces that share (charge,
    # isotope). Each feature's trace ids run 0..(#traces-1).
    assert traces.filter(pl.col("feature_id") == 0)["trace_in_feature"] \
        .unique().sort().to_list() == [0, 1]
    assert traces.filter(pl.col("feature_id") == 1)["trace_in_feature"] \
        .unique().to_list() == [0]
    # feature 12 / charge 13 / isotope 11 appears as TWO distinct traces: the dup
    # (charge, isotope) must NOT collapse -> two distinct trace_in_feature values.
    dup = traces.filter(
        (pl.col("feature_id") == 12) & (pl.col("charge") == 13)
        & (pl.col("isotope") == 11)
    )
    assert dup["trace_in_feature"].unique().sort().to_list() == [0, 1]
    # within one trace, every point shares the SAME trace_in_feature (one id/trace)
    per_trace_pts = dup.group_by("trace_in_feature").len().sort("trace_in_feature")
    assert per_trace_pts["len"].to_list() == [2, 3]
