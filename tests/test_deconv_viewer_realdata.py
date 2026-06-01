"""End-to-end verification of the Phase-1 FLASHDeconv OpenMS-Insight viewer
against the bundled real example workspace.

These construct every FLASHDeconv component through ``build_component`` from the
actual ``example_fd`` parquet caches and exercise the cross-link filtering
(scan-table click → spectra / mass-table / 3D), verifying the index→value
migration reproduces the original row counts exactly.

Everything except the final Streamlit ``render()`` call is covered (rendering
needs a browser); the data path — load → long-format explode → component
filter — is fully verified. Skipped automatically when example data or
openms_insight is absent.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

pl = pytest.importorskip("polars")
pytest.importorskip("openms_insight")

from pathlib import Path  # noqa: E402

_FD = Path(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
) / "example-data" / "workspaces" / "default" / "flashdeconv" / "cache" / "files" / "example_fd"

pytestmark = pytest.mark.skipif(
    not _FD.is_dir(), reason="example workspace data not available"
)


class _FakeFM:
    """Serves the real example_fd parquet files like FileManager.get_results."""

    def __init__(self, cache_path):
        self.cache_path = str(cache_path)

    def get_results(
        self, dataset_id, names, use_polars=False, use_pyarrow=False, partial=False
    ):
        out = {}
        for n in names:
            p = _FD / f"{n}.pq"
            out[n] = pl.scan_parquet(p) if use_polars else p
        return out

    def result_exists(self, a, b):
        return False  # no submitted sequence in this fixture


ALL_DECONV_COMPONENTS = [
    "ms1_deconv_heat_map",
    "ms2_deconv_heat_map",
    "ms1_raw_heatmap",
    "ms2_raw_heatmap",
    "scan_table",
    "mass_table",
    "deconv_spectrum",
    "anno_spectrum",
    "3D_SN_plot",
    "fdr_plot",
]


@pytest.fixture
def fake_fm(tmp_path):
    return _FakeFM(tmp_path)


@pytest.mark.parametrize("comp", ALL_DECONV_COMPONENTS)
def test_every_component_builds_from_real_cache(comp, fake_fm, monkeypatch):
    import streamlit as st

    monkeypatch.setattr(st, "session_state", {}, raising=False)
    from src.render_oi.deconv_viewer import build_component

    builder = build_component(comp, "example_fd", fake_fm, None, key_prefix="p0")
    assert callable(builder), f"{comp} did not produce a render callable"


def test_scan_click_cross_link_row_counts(fake_fm, monkeypatch, tmp_path):
    """Scan-table click filters spectra / mass-table / 3D to exact row counts."""
    import streamlit as st

    monkeypatch.setattr(st, "session_state", {}, raising=False)

    from openms_insight import LinePlot, Scatter3D, Table

    from src.parse.long_format import (
        explode_signal_peaks_long,
        explode_spectrum_long,
    )
    from src.render_oi.deconv_viewer import _explode_mass_table

    mt_long = _explode_mass_table(pl.scan_parquet(_FD / "mass_table.pq")).collect()
    busiest = mt_long.group_by("index").len().sort("len", descending=True)
    scan = int(busiest["index"][0])
    n_mass = int(busiest["len"][0])

    # Mass table filtered by scanIndex
    mass_tbl = Table(
        cache_id="x_mass",
        data=_explode_mass_table(pl.scan_parquet(_FD / "mass_table.pq")),
        filters={"scanIndex": "index"},
        interactivity={"massIndex": "mass_id"},
        index_field="mass_id",
        cache_path=str(tmp_path),
    )
    mvd = mass_tbl._prepare_vue_data({"scanIndex": scan})
    assert mvd["_pagination"]["total_rows"] == n_mass

    # Deconvolved spectrum filtered by scanIndex
    dec = explode_spectrum_long(pl.scan_parquet(_FD / "deconv_spectrum.pq"))
    lp = LinePlot(
        cache_id="x_dec",
        data=dec,
        filters={"scanIndex": "index"},
        x_column="mass",
        y_column="intensity",
        cache_path=str(tmp_path),
    )
    exp_peaks = (
        explode_spectrum_long(pl.scan_parquet(_FD / "deconv_spectrum.pq"))
        .collect()
        .filter(pl.col("index") == scan)
        .height
    )
    assert len(lp._prepare_vue_data({"scanIndex": scan})["plotData"]) == exp_peaks

    # 3D S/N: all masses for scan, isolate mass 0
    sn = explode_signal_peaks_long(pl.scan_parquet(_FD / "threedim_SN_plot.pq"))
    snc = sn.collect()
    s3 = Scatter3D(
        cache_id="x_3d",
        data=sn,
        filters={"scanIndex": "index"},
        optional_filters={"massIndex": "mass_id"},
        cache_path=str(tmp_path),
    )
    all_peaks = len(s3._prepare_vue_data({"scanIndex": scan})["scatter3dData"])
    assert all_peaks == snc.filter(pl.col("index") == scan).height
    mass0 = len(
        s3._prepare_vue_data({"scanIndex": scan, "massIndex": 0})["scatter3dData"]
    )
    assert mass0 == snc.filter(
        (pl.col("index") == scan) & (pl.col("mass_id") == 0)
    ).height
