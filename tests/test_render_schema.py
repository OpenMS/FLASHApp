"""Construct-smoke for ``src.render.schema.build_insight_caches``.

Builds synthetic FLASHApp FileManager caches (matching the ``src/parse/*`` output
schemas), runs ``build_insight_caches`` for each tool, and asserts the Insight-ready
tidy parquet is produced with the stable-ID columns and the right explode shapes.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from src.workflow.FileManager import FileManager
from src.render.schema import (
    build_insight_caches,
    _explode_list_cols,
    _explode_nested_signal_peaks,
    _comma_split_long,
    _kde_to_long,
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
    assert {"protein_id", "scan_id"}.issubset(proteins.columns)
    assert proteins["protein_id"].to_list() == [0, 1]
    # protein row carries its scan (deconv-row index): Scan 10 -> 0, Scan 20 -> 1,
    # so a protein-row click can resolve protein -> scan (value-based scan map).
    assert proteins["scan_id"].to_list() == [0, 1]

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
    assert feats["feature_id"].to_list() == [0, 1]

    traces = pl.read_parquet(fm.result_path(ds, "quant_traces"))
    assert {"feature_id", "charge", "isotope", "centroid_mz", "rt", "mz",
            "intensity"}.issubset(traces.columns)
    # feature 0: 3+2 points, feature 1: 2 points -> 7 total
    per = {r["feature_id"]: r["len"]
           for r in traces.group_by("feature_id").len().to_dicts()}
    assert per == {0: 5, 1: 2}
