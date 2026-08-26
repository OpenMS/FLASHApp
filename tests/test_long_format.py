"""
Tests for the long-format adapters used by the OpenMS-Insight migration.

FLASHApp stores spectra as arrays-per-scan and the old grid filtered by row
index. OpenMS-Insight filters by column value over long-format frames (one row
per peak). These adapters bridge the two; the migration's highest-risk change
(index->value & array explosion) lives here, so the transforms are unit-tested
directly. They are pure polars (no Streamlit), so testable without the app.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import polars as pl

from src.parse.long_format import (
    density_series_long,
    explode_combined_spectrum_long,
    explode_signal_peaks_long,
    explode_spectrum_long,
)


def _per_scan_spectrum():
    # Two scans, arrays-per-scan (MonoMass / SumIntensity), like deconv_spectrum.
    return pl.LazyFrame(
        {
            "index": [0, 1],
            "MonoMass": [[100.0, 200.0, 300.0], [400.0, 500.0]],
            "SumIntensity": [[10.0, 0.0, 30.0], [40.0, 50.0]],
        }
    )


class TestExplodeSpectrumLong:
    def test_basic_explosion(self):
        out = explode_spectrum_long(_per_scan_spectrum()).collect()
        # 3 + 2 = 5 peaks
        assert out.height == 5
        assert out.columns == ["index", "mass", "intensity", "mass_id"]

    def test_mass_id_per_scan(self):
        out = explode_spectrum_long(_per_scan_spectrum()).collect()
        scan0 = out.filter(pl.col("index") == 0).sort("mass_id")
        assert scan0["mass_id"].to_list() == [0, 1, 2]
        assert scan0["mass"].to_list() == [100.0, 200.0, 300.0]
        scan1 = out.filter(pl.col("index") == 1).sort("mass_id")
        assert scan1["mass_id"].to_list() == [0, 1]

    def test_value_filter_reproduces_iloc(self):
        """Filtering index==1 yields exactly scan 1's peaks (the iloc[1] rows)."""
        out = explode_spectrum_long(_per_scan_spectrum()).collect()
        scan1 = out.filter(pl.col("index") == 1)
        assert scan1.height == 2
        assert sorted(scan1["mass"].to_list()) == [400.0, 500.0]

    def test_mass_id_assigned_before_intensity_filter(self):
        """mass_id matches original array position even when peaks are dropped."""
        out = explode_spectrum_long(
            _per_scan_spectrum(), drop_nonpositive_intensity=True
        ).collect()
        scan0 = out.filter(pl.col("index") == 0).sort("mass_id")
        # The zero-intensity middle peak (mass_id 1) is dropped, leaving 0 and 2
        assert scan0["mass_id"].to_list() == [0, 2]
        assert scan0["mass"].to_list() == [100.0, 300.0]

    def test_row_count_preserved_without_filter(self):
        """Total exploded rows == sum of array lengths (no peaks lost)."""
        src = _per_scan_spectrum().collect()
        expected = sum(len(a) for a in src["MonoMass"].to_list())
        out = explode_spectrum_long(src.lazy()).collect()
        assert out.height == expected


class TestExplodeCombinedSpectrum:
    def test_two_series(self):
        per_scan = pl.LazyFrame(
            {
                "index": [0],
                "MonoMass": [[100.0, 200.0]],
                "SumIntensity": [[10.0, 20.0]],
                "MonoMass_Anno": [[101.0, 201.0, 301.0]],
                "SumIntensity_Anno": [[5.0, 15.0, 25.0]],
            }
        )
        deconv, anno = explode_combined_spectrum_long(per_scan)
        d = deconv.collect()
        a = anno.collect()
        assert d.height == 2
        assert a.height == 3
        assert d["mass"].to_list() == [100.0, 200.0]
        assert a["mass"].to_list() == [101.0, 201.0, 301.0]


class TestExplodeSignalPeaks:
    def test_signal_noise_explosion(self):
        # One scan, two masses; each mass has signal peaks and noisy peaks.
        # Peak format: [_, mz, intensity, charge]
        # Peaks arrive as uniform-float lists (pyarrow float columns).
        per_scan = pl.LazyFrame(
            {
                "index": [5],
                "SignalPeaks": [
                    [
                        [[0.0, 500.0, 1000.0, 2.0], [0.0, 510.0, 800.0, 2.0]],  # mass 0
                        [[0.0, 700.0, 1500.0, 3.0]],                            # mass 1
                    ]
                ],
                "NoisyPeaks": [
                    [
                        [[0.0, 505.0, 50.0, 2.0]],  # mass 0
                        [],                          # mass 1 (no noise)
                    ]
                ],
            }
        )
        out = explode_signal_peaks_long(per_scan).collect()
        # 3 signal + 1 noise = 4 peaks
        assert out.height == 4
        assert set(out.columns) == {
            "index",
            "mass_id",
            "mz",
            "charge",
            "intensity",
            "kind",
        }
        signal = out.filter(pl.col("kind") == "signal")
        noise = out.filter(pl.col("kind") == "noise")
        assert signal.height == 3
        assert noise.height == 1
        # massIndex isolation: mass 0 has 2 signal + 1 noise
        mass0 = out.filter(pl.col("mass_id") == 0)
        assert mass0.height == 3
        # mass 1 has 1 signal, 0 noise
        mass1 = out.filter(pl.col("mass_id") == 1)
        assert mass1.height == 1
        assert mass1["kind"].to_list() == ["signal"]

    def test_charge_and_intensity_extracted(self):
        per_scan = pl.LazyFrame(
            {
                "index": [0],
                "SignalPeaks": [[[[0.0, 500.0, 1000.0, 4.0]]]],
                "NoisyPeaks": [[[]]],
            }
        )
        out = explode_signal_peaks_long(per_scan).collect()
        assert out.height == 1
        row = out.row(0, named=True)
        assert row["mz"] == 500.0
        assert row["intensity"] == 1000.0
        assert row["charge"] == 4.0


class TestDensitySeriesLong:
    def test_stacks_target_and_decoy(self):
        target = pl.DataFrame({"x": [0.1, 0.2], "y": [1.0, 2.0]})
        decoy = pl.DataFrame({"x": [0.3], "y": [0.5]})
        out = density_series_long(target, decoy)
        assert out.columns == ["series", "x", "y"]
        assert out.filter(pl.col("series") == "Target").height == 2
        assert out.filter(pl.col("series") == "Decoy").height == 1

    def test_empty_decoy(self):
        target = pl.DataFrame({"x": [0.1], "y": [1.0]})
        empty = pl.DataFrame(schema={"x": pl.Float64, "y": pl.Float64})
        out = density_series_long(target, empty)
        assert out.filter(pl.col("series") == "Decoy").height == 0
        assert out.filter(pl.col("series") == "Target").height == 1

    def test_no_data(self):
        out = density_series_long(None, None)
        assert out.height == 0
        assert out.columns == ["series", "x", "y"]
