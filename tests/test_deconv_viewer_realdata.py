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


def _capture_built_components(monkeypatch):
    """Patch OI component __call__ so building a render closure and invoking it
    records the constructed component instance + its render kwargs, without a
    browser. Returns a dict {component_type: (instance, call_kwargs)}.
    """
    import openms_insight as oi

    captured = {}

    def make_spy(cls, name):
        orig = cls.__call__

        def spy(self, *args, **kwargs):  # noqa: ANN001
            captured[name] = (self, kwargs)
            return None

        monkeypatch.setattr(cls, "__call__", spy, raising=False)
        return orig

    make_spy(oi.Table, "table")
    make_spy(oi.LinePlot, "lineplot")
    make_spy(oi.Scatter3D, "scatter3d")
    return captured


def test_scan_table_column_definitions(fake_fm, monkeypatch):
    """scan_table passes legacy column_definitions covering every real field."""
    import streamlit as st

    monkeypatch.setattr(st, "session_state", {}, raising=False)
    from openms_insight import StateManager

    from src.render_oi.deconv_viewer import build_component

    captured = _capture_built_components(monkeypatch)
    sm = StateManager(session_key="oi_scan_coldef")
    build_component("scan_table", "example_fd", fake_fm, sm, key_prefix="p0")()

    tbl, _ = captured["table"]
    defs = tbl._column_definitions
    assert defs, "scan_table must pass explicit column_definitions"
    by_field = {c["field"]: c for c in defs}
    # Every real scan_table.pq column is covered (no column dropped).
    real_cols = set(
        pl.scan_parquet(_FD / "scan_table.pq").collect_schema().names()
    )
    assert real_cols <= set(by_field), (real_cols, set(by_field))
    # Legacy titles and descriptive tooltips.
    assert by_field["MSLevel"]["title"] == "MS Level"
    assert by_field["RT"]["title"] == "Retention time"
    assert by_field["PrecursorMass"]["title"] == "Precursor Mass"
    assert by_field["#Masses"]["title"] == "#Masses"
    assert all(isinstance(c["headerTooltip"], str) for c in defs)
    # Numeric formatting on the float columns (RT / PrecursorMass).
    assert by_field["RT"]["formatter"] == "money"
    assert by_field["PrecursorMass"]["formatter"] == "money"


def test_mass_table_column_definitions(fake_fm, monkeypatch):
    """mass_table passes legacy column_definitions covering every exploded field."""
    import streamlit as st

    monkeypatch.setattr(st, "session_state", {}, raising=False)
    from openms_insight import StateManager

    from src.render_oi.deconv_viewer import _explode_mass_table, build_component

    captured = _capture_built_components(monkeypatch)
    sm = StateManager(session_key="oi_mass_coldef")
    build_component("mass_table", "example_fd", fake_fm, sm, key_prefix="p0")()

    tbl, _ = captured["table"]
    defs = tbl._column_definitions
    assert defs, "mass_table must pass explicit column_definitions"
    by_field = {c["field"]: c for c in defs}
    # Every exploded per-mass field is shown EXCEPT the scan "index" — that
    # column is the cross-link filter target (filters={scanIndex: index}) and was
    # not a visible column in the legacy Mass Table; it stays available for
    # filtering (Table._get_columns_to_select adds filter columns) but is not
    # displayed. The displayed row index is mass_id ("Index"), matching legacy.
    exploded_cols = set(
        _explode_mass_table(pl.scan_parquet(_FD / "mass_table.pq"))
        .collect_schema()
        .names()
    )
    assert (exploded_cols - {"index"}) <= set(by_field), (exploded_cols, set(by_field))
    assert "index" not in by_field  # scan index is a filter column, not displayed
    assert "mass_id" in by_field  # displayed per-mass index ("Index")
    assert tbl._filters == {"scanIndex": "index"}  # filter column still present
    assert by_field["MonoMass"]["title"] == "Monoisotopic mass"
    assert by_field["SumIntensity"]["title"] == "Sum intensity"
    assert by_field["MinCharges"]["title"] == "Min charge"
    assert by_field["MaxCharges"]["title"] == "Max charge"
    assert by_field["CosineScore"]["title"] == "Cosine score"
    assert by_field["QScore"]["title"] == "QScore"
    # Numeric formatting where sensible (masses / intensity / scores).
    for f in ("MonoMass", "SumIntensity", "CosineScore", "SNR", "QScore"):
        assert by_field[f]["formatter"] == "money", f


def test_deconv_spectrum_no_always_on_annotation(fake_fm, monkeypatch):
    """The deconvolved (mass-axis) spectrum carries NO always-on per-peak label:
    the legacy draws annotation boxes only for the SELECTED mass and gates z=N
    charge labels to the m/z axis. It relies on the massIndex cross-link for the
    gold selected-peak highlight, and its rows carry mass_id (the highlight key,
    aligned with the mass table)."""
    import streamlit as st

    monkeypatch.setattr(st, "session_state", {}, raising=False)
    from openms_insight import StateManager

    from src.render_oi.deconv_viewer import MASS, build_component

    captured = _capture_built_components(monkeypatch)
    sm = StateManager(session_key="oi_deconv_ann")
    build_component("deconv_spectrum", "example_fd", fake_fm, sm, key_prefix="p0")()

    lp, _ = captured["lineplot"]
    # No always-on labels (the regression was a z=N label on every stick).
    assert lp._annotation_column is None
    # The mass-table <-> spectrum highlight cross-link is wired via massIndex.
    assert lp._interactivity == {MASS: "mass_id"}

    # mass_id is present (the highlight lookup key) and rows are non-empty.
    from src.parse.long_format import explode_spectrum_long

    df = explode_spectrum_long(pl.scan_parquet(_FD / "deconv_spectrum.pq")).collect()
    assert "mass_id" in df.columns and df.height > 0


def test_anno_spectrum_has_no_charge_annotation(fake_fm, monkeypatch):
    """anno_spectrum has no charge data, so it carries no annotation column."""
    import streamlit as st

    monkeypatch.setattr(st, "session_state", {}, raising=False)
    from openms_insight import StateManager

    from src.render_oi.deconv_viewer import build_component

    captured = _capture_built_components(monkeypatch)
    sm = StateManager(session_key="oi_anno_ann")
    build_component("anno_spectrum", "example_fd", fake_fm, sm, key_prefix="p0")()

    lp, _ = captured["lineplot"]
    assert lp._annotation_column is None


def test_3d_plot_title_reflects_selection(fake_fm, monkeypatch):
    """3D_SN_plot title is 'Precursor signals' with no mass selected and
    'Mass signals' once a massIndex is selected."""
    import streamlit as st

    monkeypatch.setattr(st, "session_state", {}, raising=False)
    from openms_insight import StateManager

    from src.render_oi.deconv_viewer import build_component

    captured = _capture_built_components(monkeypatch)
    sm = StateManager(session_key="oi_3d_title")
    render = build_component("3D_SN_plot", "example_fd", fake_fm, sm, key_prefix="p0")

    sm.set_selection("scanIndex", 0)
    render()
    assert captured["scatter3d"][0]._title == "Precursor signals"

    sm.set_selection("massIndex", 0)
    render()
    assert captured["scatter3d"][0]._title == "Mass signals"

    sm.set_selection("massIndex", None)
    render()
    assert captured["scatter3d"][0]._title == "Precursor signals"


def test_scan_change_resets_mass(monkeypatch):
    """A scan-table click resets massIndex to the new scan's first peak (legacy
    updateSelectedMass(0)); a heatmap click that changes scan AND mass together
    keeps the heatmap's mass instead of being clobbered."""
    import streamlit as st

    monkeypatch.setattr(st, "session_state", {}, raising=False)
    from openms_insight import StateManager

    from src.render_oi.deconv_viewer import MASS, SCAN, _reset_mass_on_scan_change

    sm = StateManager(session_key="oi_reset_test")
    sm.set_selection(SCAN, 1)
    sm.set_selection(MASS, 5)
    _reset_mass_on_scan_change(sm)  # first call only records the markers
    assert sm.get_selection(MASS) == 5

    # Scan-table click: scan changes, mass unchanged -> reset to first peak.
    sm.set_selection(SCAN, 2)
    _reset_mass_on_scan_change(sm)
    assert sm.get_selection(MASS) == 0

    # Heatmap click: scan AND mass change together -> keep the heatmap's mass.
    sm.set_selection(SCAN, 3)
    sm.set_selection(MASS, 9)
    _reset_mass_on_scan_change(sm)
    assert sm.get_selection(MASS) == 9
