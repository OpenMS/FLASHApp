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
    assert df.columns == ["index", "peak_id", "MonoMass", "SumIntensity"]
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
    assert df.columns == ["index", "peak_id", "MonoMass", "SumIntensity", "is_signal"]
    dpd = deconv.to_pandas()
    for r in df.iter_rows(named=True):
        sp = dpd[dpd["index"] == r["index"]].iloc[0]["SignalPeaks"]
        pid = r["peak_id"]
        want = (pid < len(sp)) and (len(sp[pid]) > 0)
        assert bool(r["is_signal"]) == want
    # Ragged past-end position (scan 0, peak_id 3) must be is_signal False.
    row3 = df.filter((pl.col("index") == 0) & (pl.col("peak_id") == 3)).row(0, named=True)
    assert row3["is_signal"] is False


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
