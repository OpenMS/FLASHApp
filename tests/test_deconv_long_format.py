"""
Tests for the Stage B long-format producers in src/parse/deconv.py.

FLASHApp's legacy render path filters per-scan spectra/masses by ROW INDEX
(``iloc[scanIndex]``) and stores arrays-per-scan. The OpenMS-Insight components
filter by COLUMN VALUE and expect LONG format (one row per peak/mass). These
tests validate the additive long-format producers:

  - row-count fidelity: exploded rows == legacy per-column max-length expansion
    (the TabulatorMassTable.vue ``forEach``-per-column semantics);
  - index-filter parity: ``filter(index == k)`` reproduces, position by position,
    the legacy ``iloc[k]`` array contents (with shorter columns padded to None);
  - ragged-scan handling: when ``mz_array`` (full spectrum) is longer than the
    per-mass charge/score arrays, trailing rows carry mass/intensity with null
    charge/score cells; ``is_signal`` is False past the SignalPeaks axis;
  - id columns: ``peak_id`` / ``mass_id`` are 0-based positions within each scan,
    and the deconv peak axis and mass-table mass axis are aligned.

The producers are pure polars (no Streamlit / pyopenms), so they are unit
testable without booting the app. ``pyopenms`` is stubbed at import time only
because src/parse/deconv.py imports src/parse/masstable.py at module load.
"""

import os
import sys
import types

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Stub pyopenms so importing src.parse.deconv works without the native dep
# (the long-format producers do not use it).
if "pyopenms" not in sys.modules:
    _m = types.ModuleType("pyopenms")
    for _a in ("MSExperiment", "MzMLFile", "SpectrumLookup", "Constants"):
        setattr(_m, _a, type(_a, (), {"PROTON_MASS_U": 1.0, "C13C12_MASSDIFF_U": 1.0}))
    sys.modules["pyopenms"] = _m

import polars as pl

from src.parse.deconv import (
    anno_spectrum_long,
    combined_spectrum_long,
    deconv_spectrum_long,
    mass_table_long,
    threedim_SN_plot,
)


def _deconv():
    # scan 0 is RAGGED: mz_array length 4 > per-mass arrays length 3.
    # scan 1 is EMPTY. scans 2,3 have 2 and 1 masses.
    return pl.DataFrame(
        {
            "mz_array": [[1000.1, 2000.2, 3000.3, 4000.4], [], [500.5, 600.6], [777.7]],
            "intensity_array": [[10.0, 20.0, 30.0, 40.0], [], [5.0, 6.0], [7.0]],
            "MinCharges": [[1, 2, 3], [], [1, 2], [4]],
            "MaxCharges": [[5, 6, 7], [], [3, 4], [8]],
            "MinIsotopes": [[0, 1, 2], [], [0, 1], [3]],
            "MaxIsotopes": [[4, 5, 6], [], [2, 3], [7]],
            "cos": [[0.9, 0.8, 0.7], [], [0.95, 0.85], [0.6]],
            "snr": [[11.0, 12.0, 13.0], [], [14.0, 15.0], [16.0]],
            "qscore": [[0.99, 0.98, 0.97], [], [0.96, 0.95], [0.94]],
            "SignalPeaks": [
                [[[0.0, 1000.1, 10.0, 1.0]], [], [[2.0, 3000.3, 30.0, 3.0]]],
                [],
                [[[0.0, 500.5, 5.0, 1.0]], []],
                [[]],
            ],
        }
    ).with_row_index("index")


def _anno():
    return pl.DataFrame(
        {
            "mz_array": [[101.1, 102.2], [201.1], [], [401.1, 402.2, 403.3]],
            "intensity_array": [[1.0, 2.0], [3.0], [], [4.0, 5.0, 6.0]],
        }
    ).with_row_index("index")


def _max_len_expansion(row, cols):
    arrays = {c: list(row[c]) for c in cols}
    n = max((len(a) for a in arrays.values()), default=0)
    return [{c: (arrays[c][i] if i < len(arrays[c]) else None) for c in arrays} for i in range(n)]


def test_deconv_spectrum_long_schema_and_rowcount():
    df = deconv_spectrum_long(_deconv().lazy()).collect()
    assert df.columns == [
        "index", "peak_id", "MonoMass", "SumIntensity",
        "signal_mzs", "signal_charges", "signal_intensities",
    ]
    # 4 + 0 + 2 + 1 = 7
    assert df.height == 7


def test_anno_spectrum_long_index_filter_parity():
    anno = _anno()
    df = anno_spectrum_long(anno.lazy()).collect()
    assert df.columns == ["index", "peak_id", "MonoMass_Anno", "SumIntensity_Anno"]
    apd = anno.to_pandas()
    for k in range(len(apd)):
        sub = df.filter(pl.col("index") == k).sort("peak_id")
        want_mass = list(apd.iloc[k]["mz_array"])
        assert sub["MonoMass_Anno"].to_list() == want_mass
        assert sub["peak_id"].to_list() == list(range(len(want_mass)))


def test_mass_table_long_ragged_padding():
    deconv = _deconv()
    df = mass_table_long(deconv.lazy()).collect()
    expected_cols = [
        "index", "mass_id", "MonoMass", "SumIntensity",
        "MinCharges", "MaxCharges", "MinIsotopes", "MaxIsotopes",
        "CosineScore", "SNR", "QScore",
    ]
    assert df.columns == expected_cols
    # Scan 0 has 4 mass/intensity positions but only 3 charge positions →
    # row at mass_id 3 carries MonoMass=4000.4 with null MinCharges.
    scan0 = df.filter(pl.col("index") == 0).sort("mass_id")
    assert scan0.height == 4
    last = scan0.row(3, named=True)
    assert last["MonoMass"] == 4000.4
    assert last["MinCharges"] is None
    assert last["QScore"] is None
    # First three rows have full charge/score data.
    assert scan0.row(0, named=True)["MinCharges"] == 1


def test_mass_table_long_rowcount_and_empty_scan():
    df = mass_table_long(_deconv().lazy()).collect()
    # max-length per scan: 4 + 0 + 2 + 1 = 7
    assert df.height == 7
    # Empty scan contributes no rows.
    assert df.filter(pl.col("index") == 1).height == 0


def test_combined_spectrum_long_is_signal():
    deconv = _deconv()
    df = combined_spectrum_long(deconv.lazy()).collect()
    assert df.columns == [
        "index", "peak_id", "MonoMass", "SumIntensity", "is_signal",
        "signal_mzs", "signal_charges", "signal_intensities",
    ]
    dpd = deconv.to_pandas()
    for r in df.iter_rows(named=True):
        sp = dpd[dpd["index"] == r["index"]].iloc[0]["SignalPeaks"]
        pid = r["peak_id"]
        want = (pid < len(sp)) and (len(sp[pid]) > 0)
        assert bool(r["is_signal"]) == want
    # Ragged past-end position (scan 0, peak_id 3) must be is_signal False.
    row3 = df.filter((pl.col("index") == 0) & (pl.col("peak_id") == 3)).row(0, named=True)
    assert row3["is_signal"] is False


# Signal-peak quadruple layout in the nested SignalPeaks column:
# [peak_index, mz, intensity, charge]. (Verified against the example caches:
# SignalPeaks has dtype List(List(List(Float64))) — scan -> mass -> peak -> tuple.)
_SP_MZ, _SP_INT, _SP_CH = 1, 2, 3


def _check_signal_arrays(df):
    """Shared assertions for the per-mass signal_* list columns on a long frame."""
    # Column dtypes: lists of f64 / i64 / f64.
    assert df.schema["signal_mzs"] == pl.List(pl.Float64)
    assert df.schema["signal_intensities"] == pl.List(pl.Float64)
    assert df.schema["signal_charges"] == pl.List(pl.Int64)

    deconv = _deconv()
    dpd = deconv.to_pandas()
    for r in df.iter_rows(named=True):
        sp = dpd[dpd["index"] == r["index"]].iloc[0]["SignalPeaks"]
        pid = r["peak_id"]
        peaks = list(sp[pid]) if pid < len(sp) else []

        mzs = r["signal_mzs"]
        chs = r["signal_charges"]
        ints = r["signal_intensities"]

        # Never null — past-end / no-signal positions are empty lists.
        assert mzs is not None and chs is not None and ints is not None
        # The three arrays are mutually aligned (one entry per signal peak).
        assert len(mzs) == len(chs) == len(ints) == len(peaks)
        # Contents match the per-mass signal peaks at this position.
        assert mzs == [p[_SP_MZ] for p in peaks]
        assert ints == [p[_SP_INT] for p in peaks]
        assert chs == [int(p[_SP_CH]) for p in peaks]

        # Alignment with is_signal (combined frame only) / non-emptiness.
        if "is_signal" in r:
            assert bool(r["is_signal"]) == (len(mzs) > 0)


def test_deconv_spectrum_long_signal_arrays():
    df = deconv_spectrum_long(_deconv().lazy()).collect()
    _check_signal_arrays(df)
    # Concrete spot check: scan 0, peak 0 has one signal peak (mz 1000.1, ch 1).
    r0 = df.filter((pl.col("index") == 0) & (pl.col("peak_id") == 0)).row(0, named=True)
    assert r0["signal_mzs"] == [1000.1]
    assert r0["signal_charges"] == [1]
    assert r0["signal_intensities"] == [10.0]
    # Non-signal peak (scan 0, peak 1) and ragged past-end peak (scan 0, peak 3)
    # both carry empty lists.
    for pid in (1, 3):
        rr = df.filter((pl.col("index") == 0) & (pl.col("peak_id") == pid)).row(0, named=True)
        assert rr["signal_mzs"] == []
        assert rr["signal_charges"] == []
        assert rr["signal_intensities"] == []


def test_combined_spectrum_long_signal_arrays():
    df = combined_spectrum_long(_deconv().lazy()).collect()
    _check_signal_arrays(df)
    # Signal-flagged peaks have non-empty, equal-length signal_* lists; non-signal
    # peaks have empty lists across all three.
    for r in df.iter_rows(named=True):
        if r["is_signal"]:
            assert len(r["signal_mzs"]) > 0
            assert len(r["signal_mzs"]) == len(r["signal_charges"]) == len(r["signal_intensities"])
        else:
            assert r["signal_mzs"] == []
            assert r["signal_charges"] == []
            assert r["signal_intensities"] == []


def _deconv_3d():
    # MS1 precursor scan (scan 100) with two masses; MS2 fragment scan (scan 101)
    # isolated from precursor mass 2000.2 in scan 100.
    return pl.DataFrame(
        {
            "Scan": [100, 101],
            "PrecursorScan": [0.0, 100.0],
            "PrecursorMass": [0.0, 2000.2],
            "mz_array": [[1000.1, 2000.2], [3000.3]],
            "intensity_array": [[10.0, 20.0], [30.0]],
            "SignalPeaks": [
                [[[0.0, 1000.1, 10.0, 1.0]], [[1.0, 2000.2, 20.0, 2.0]]],
                [[[0.0, 3000.3, 30.0, 1.0]]],
            ],
            "NoisyPeaks": [
                [[], []],
                [[]],
            ],
        }
    ).with_row_index("index")


def test_threedim_SN_plot_precursor_lookup_columns():
    df = threedim_SN_plot(_deconv_3d().lazy()).collect()
    assert df.columns == [
        "index", "Scan", "PrecursorScan", "PrecursorMass",
        "MonoMass", "SignalPeaks", "NoisyPeaks",
    ]
    # Precursor-lookup key dtypes.
    assert df.schema["Scan"] == pl.Int64
    assert df.schema["PrecursorScan"] == pl.Float64
    assert df.schema["PrecursorMass"] == pl.Float64
    # MonoMass is the per-mass array (== mz_array).
    assert df.schema["MonoMass"] == pl.List(pl.Float64)

    # The MS2 fragment scan's precursor resolves to a mass in its precursor scan:
    # find the precursor-scan row (Scan == PrecursorScan) and the MonoMass index
    # matching PrecursorMass — the position the Scatter3D uses for SignalPeaks.
    ms2 = df.filter(pl.col("Scan") == 101).row(0, named=True)
    assert ms2["PrecursorScan"] == 100.0
    assert ms2["PrecursorMass"] == 2000.2
    prec = df.filter(pl.col("Scan") == int(ms2["PrecursorScan"])).row(0, named=True)
    mass_index = prec["MonoMass"].index(ms2["PrecursorMass"])
    assert mass_index == 1
    # That per-mass position carries the matching signal peaks in the precursor scan.
    assert len(prec["SignalPeaks"][mass_index]) > 0


def test_peak_id_and_mass_id_share_mass_axis():
    deconv = _deconv()
    ds = deconv_spectrum_long(deconv.lazy()).collect()
    mt = mass_table_long(deconv.lazy()).collect()
    join = ds.join(
        mt.select(["index", "mass_id", "MonoMass"]),
        left_on=["index", "peak_id"],
        right_on=["index", "mass_id"],
        how="inner",
        suffix="_mt",
    )
    assert (join["MonoMass"] == join["MonoMass_mt"]).all()
