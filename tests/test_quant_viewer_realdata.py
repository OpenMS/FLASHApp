"""End-to-end verification of the Phase-3 FLASHQuant OpenMS-Insight viewer
against the bundled real quant workspace.

Verifies the long-format trace explosion preserves point counts and the
FeatureView filters to exact per-feature-group points.

Skipped automatically when example data or openms_insight is absent.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

pl = pytest.importorskip("polars")
pytest.importorskip("openms_insight")

from pathlib import Path  # noqa: E402

_FQ = (
    Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    / "example-data" / "workspaces" / "default" / "flashquant"
    / "cache" / "files" / "example"
)

pytestmark = pytest.mark.skipif(
    not (_FQ / "quant_dfs.pq").exists(),
    reason="example FLASHQuant workspace data not available",
)


class _FakeFM:
    def __init__(self, cache_path):
        self.cache_path = str(cache_path)

    def get_results(
        self, dataset_id, names, use_polars=False, use_pyarrow=False, partial=False
    ):
        return {
            n: (pl.scan_parquet(_FQ / f"{n}.pq") if use_polars else _FQ / f"{n}.pq")
            for n in names
        }

    def result_exists(self, a, b):
        return True


@pytest.fixture
def fake_fm(tmp_path):
    return _FakeFM(tmp_path)


def test_quant_long_format_preserves_points():
    """Exploding the per-group arrays preserves every trace point."""
    from src.parse.long_format import explode_quant_traces_long

    q = pl.read_parquet(_FQ / "quant_dfs.pq")
    long = explode_quant_traces_long(q.lazy()).collect()

    # One unique feature group per source row
    assert long["feature_group"].n_unique() == q.height

    # Group 0: total points == sum of split lengths over its traces
    g0 = q.row(0, named=True)
    expected_pts = sum(len(s.split(",")) for s in g0["RTs"])
    got_pts = long.filter(pl.col("feature_group") == 0).height
    assert got_pts == expected_pts


def test_quant_view_builds(fake_fm, monkeypatch):
    import streamlit as st

    monkeypatch.setattr(st, "session_state", {}, raising=False)
    from src.render_oi.quant_viewer import build_quant_components

    render = build_quant_components("example", fake_fm, None, key_prefix="p0")
    assert callable(render)


def test_feature_view_filters_to_group(fake_fm, monkeypatch, tmp_path):
    import streamlit as st

    monkeypatch.setattr(st, "session_state", {}, raising=False)
    from openms_insight import FeatureView

    from src.parse.long_format import explode_quant_traces_long

    traces = explode_quant_traces_long(pl.scan_parquet(_FQ / "quant_dfs.pq"))
    tc = traces.collect()
    fv = FeatureView(
        cache_id="fv_x",
        data=traces,
        filters={"featureGroup": "feature_group"},
        charge_column="charge",
        mz_column="mz",
        rt_column="rt",
        intensity_column="intensity",
        isotope_column="isotope",
        cache_path=str(tmp_path),
    )
    for fg in (0, 5, 100):
        vd = fv._prepare_vue_data({"featureGroup": fg})
        data_key = next(k for k in vd if not k.startswith("_"))
        assert len(vd[data_key]) == tc.filter(pl.col("feature_group") == fg).height
