"""End-to-end verification of the Phase-2 FLASHTnT OpenMS-Insight viewer
against the bundled real antibody workspace.

Builds every FLASHTnT component through ``build_component_tnt`` from the real
caches and verifies the proteoform→scan resolution plus the cross-link filters
(protein click → tag table / sequence view / combined spectrum).

Skipped automatically when example data or openms_insight is absent.
"""

import gzip
import os
import pickle
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

pl = pytest.importorskip("polars")
pytest.importorskip("openms_insight")

from pathlib import Path  # noqa: E402

_TNT = (
    Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    / "example-data" / "workspaces" / "default" / "flashtnt" / "cache" / "files"
)


def _first_dataset():
    if not _TNT.is_dir():
        return None
    for d in sorted(_TNT.iterdir()):
        if (d / "protein_dfs.pq").exists():
            return d
    return None


_DS = _first_dataset()

pytestmark = pytest.mark.skipif(
    _DS is None, reason="example FLASHTnT workspace data not available"
)


class _FakeFM:
    def __init__(self, cache_path):
        self.cache_path = str(cache_path)

    def get_results(
        self, dataset_id, names, use_polars=False, use_pyarrow=False, partial=False
    ):
        out = {}
        for n in names:
            pq = _DS / f"{n}.pq"
            pkl = _DS / f"{n}.pkl.gz"
            if pq.exists():
                out[n] = pl.scan_parquet(pq) if use_polars else pq
            elif pkl.exists():
                out[n] = pkl
            elif not partial:
                raise KeyError(n)
        return out

    def result_exists(self, a, b):
        return True


ALL_TNT_COMPONENTS = [
    "protein_table",
    "tag_table",
    "combined_spectrum",
    "sequence_view",
    "id_fdr_plot",
    "ms1_raw_heatmap",
    "ms1_deconv_heat_map",
]


@pytest.fixture
def fake_fm(tmp_path):
    return _FakeFM(tmp_path)


@pytest.mark.parametrize("comp", ALL_TNT_COMPONENTS)
def test_every_tnt_component_builds(comp, fake_fm, monkeypatch):
    import streamlit as st

    monkeypatch.setattr(st, "session_state", {}, raising=False)
    from src.render_oi.tnt_viewer import build_component_tnt

    builder = build_component_tnt(comp, "ds", fake_fm, None, key_prefix="p0")
    assert callable(builder), f"{comp} did not produce a render callable"


def _captured_table(builder, monkeypatch):
    """Run a build_component_tnt render callable and capture the Table instance it
    constructs (the protein table builds its Table inside the render closure to
    honour the runtime "Best per spectrum" checkbox)."""
    from openms_insight import Table

    captured = {}

    def _spy(self, *a, **k):
        captured["table"] = self
        return None

    monkeypatch.setattr(Table, "__call__", _spy, raising=False)
    builder()
    return captured["table"]


def test_tag_table_column_definitions_and_sort(fake_fm, monkeypatch):
    """Tag table: legacy titles incl. restored Nmass/Cmass, initial_sort by Score
    desc, go_to fields, plus the proteinIndex filter / tagIndex interactivity."""
    import streamlit as st

    monkeypatch.setattr(st, "session_state", {}, raising=False)
    from openms_insight import StateManager

    from src.render_oi.tnt_viewer import build_component_tnt

    sm = StateManager(session_key="oi_tnt_test_tag")
    builder = build_component_tnt("tag_table", "ds", fake_fm, sm, key_prefix="p0")
    tbl = _captured_table(builder, monkeypatch)
    args = tbl._get_component_args()

    titles = {c["title"]: c["field"] for c in args["columnDefinitions"]}
    # Restored columns that the pre-parity build dropped.
    assert titles.get("N mass") == "Nmass"
    assert titles.get("C mass") == "Cmass"
    # Every column definition carries a header tooltip.
    assert all("headerTooltip" in c for c in args["columnDefinitions"])
    assert args["initialSort"] == [{"column": "Score", "dir": "desc"}]
    assert "Scan" in args["goToFields"] and "TagSequence" in args["goToFields"]
    assert tbl._filters == {"proteinIndex": "ProteinIndex"}
    assert args["interactivity"] == {"tagIndex": "TagIndex"}
    # The private dash-sentinel marker never reaches the frontend.
    assert all("_dashSentinel" not in c for c in args["columnDefinitions"])


def test_protein_table_column_definitions_and_sort(fake_fm, monkeypatch):
    """Protein table: legacy + parity columns, initial_sort by Score desc, go_to
    fields, and the proteinIndex interactivity preserved."""
    import streamlit as st

    monkeypatch.setattr(st, "session_state", {}, raising=False)
    monkeypatch.setattr(st, "checkbox", lambda *a, **k: True, raising=False)
    from openms_insight import StateManager

    from src.render_oi.tnt_viewer import build_component_tnt

    sm = StateManager(session_key="oi_tnt_test_prot")
    builder = build_component_tnt("protein_table", "ds", fake_fm, sm, key_prefix="p0")
    tbl = _captured_table(builder, monkeypatch)
    args = tbl._get_component_args()

    fields = {c["field"] for c in args["columnDefinitions"]}
    # Legacy columns the pre-parity build had dropped, plus the requested Coverage.
    for restored in (
        "Scan",
        "accession",
        "description",
        "length",
        "ProteoformMass",
        "Coverage(%)",
        "MatchingFragments",
        "ModCount",
        "TagCount",
        "Score",
        "ProteoformLevelQvalue",
    ):
        assert restored in fields, f"protein column {restored} missing"
    assert all("headerTooltip" in c for c in args["columnDefinitions"])
    assert args["initialSort"] == [{"column": "Score", "dir": "desc"}]
    assert args["goToFields"] == ["Scan", "accession"]
    assert args["interactivity"] == {"proteinIndex": "index"}
    assert all("_dashSentinel" not in c for c in args["columnDefinitions"])


def test_protein_best_per_spectrum_reduces_rows(fake_fm, monkeypatch):
    """The "Best per spectrum" toggle collapses the protein table to one (top
    Score) row per Scan; default-on shows fewer rows than the unfiltered table,
    and the reduced count equals the number of distinct scans."""
    import streamlit as st

    monkeypatch.setattr(st, "session_state", {}, raising=False)
    from openms_insight import StateManager

    from src.render_oi.tnt_viewer import build_component_tnt

    def _rows(best_only):
        monkeypatch.setattr(st, "checkbox", lambda *a, **k: best_only, raising=False)
        sm = StateManager(session_key=f"oi_tnt_test_bps_{best_only}")
        builder = build_component_tnt(
            "protein_table", "ds", fake_fm, sm, key_prefix="p0"
        )
        tbl = _captured_table(builder, monkeypatch)
        return tbl._prepare_vue_data({})["_pagination"]["total_rows"]

    all_rows = _rows(False)
    best_rows = _rows(True)

    n_scans = (
        pl.scan_parquet(_DS / "protein_dfs.pq").select("Scan").collect()["Scan"].n_unique()
    )
    assert best_rows == n_scans
    assert best_rows < all_rows, "Best per spectrum should reduce the row count"


def test_protein_dash_sentinel_nulled(fake_fm, monkeypatch):
    """The -1.0 'unmatched' sentinel for the dash columns (ProteoformMass /
    ProteoformLevelQvalue) is nulled at the data layer so the cell renders blank
    instead of -1, while the column stays numeric."""
    import streamlit as st

    monkeypatch.setattr(st, "session_state", {}, raising=False)
    monkeypatch.setattr(st, "checkbox", lambda *a, **k: False, raising=False)
    from openms_insight import StateManager

    from src.render_oi.tnt_viewer import build_component_tnt

    raw = pl.scan_parquet(_DS / "protein_dfs.pq").collect()
    if "ProteoformMass" not in raw.columns:
        pytest.skip("ProteoformMass not present in this dataset")
    neg1 = int((raw["ProteoformMass"] == -1).sum())
    if neg1 == 0:
        pytest.skip("no -1 sentinel present to null in this dataset")

    sm = StateManager(session_key="oi_tnt_test_dash")
    builder = build_component_tnt("protein_table", "ds", fake_fm, sm, key_prefix="p0")
    tbl = _captured_table(builder, monkeypatch)
    df = tbl._prepare_vue_data({})["tableData"]
    # No -1 sentinel survives; the rows that had it are now null.
    assert (df["ProteoformMass"] == -1).sum() == 0
    assert df["ProteoformMass"].isna().sum() >= neg1


def test_proteoform_scan_resolution(fake_fm, monkeypatch):
    import streamlit as st

    monkeypatch.setattr(st, "session_state", {}, raising=False)
    from src.render_oi.tnt_viewer import _build_proteoform_scan_map

    scan_map = _build_proteoform_scan_map(fake_fm, "ds")
    # Every proteoform that resolves carries both scan and deconv_index.
    assert scan_map, "scan map should not be empty"
    for entry in scan_map.values():
        assert "scan" in entry and "deconv_index" in entry


def test_protein_click_cross_links(fake_fm, monkeypatch, tmp_path):
    """Protein click → tag table filters; sequence view resolves coverage."""
    import streamlit as st

    monkeypatch.setattr(st, "session_state", {}, raising=False)

    from openms_insight import SequenceView, Table

    from src.render_oi.tnt_viewer import _build_proteoform_scan_map, _sequence_table

    scan_map = _build_proteoform_scan_map(fake_fm, "ds")
    pid = sorted(scan_map.keys())[0]

    # Tag table filtered by proteinIndex
    tags = pl.scan_parquet(_DS / "tag_dfs.pq")
    cols = [
        c
        for c in ["TagIndex", "TagSequence", "ProteinIndex", "StartPos", "EndPos"]
        if c in tags.collect_schema().names()
    ]
    tt = Table(
        cache_id="tt_x",
        data=tags.select(cols),
        filters={"proteinIndex": "ProteinIndex"},
        index_field="TagIndex",
        cache_path=str(tmp_path),
    )
    got = tt._prepare_vue_data({"proteinIndex": pid})["_pagination"]["total_rows"]
    expected = tags.collect().filter(pl.col("ProteinIndex") == pid).height
    assert got == expected

    # Sequence view filtered by proteinIndex carries coverage
    seq_tbl = _sequence_table(fake_fm, "ds")
    if seq_tbl is not None:
        sv = SequenceView(
            cache_id="sv_x",
            sequence_data=seq_tbl,
            filters={"proteinIndex": "proteoform_index"},
            deconvolved=True,
            coverage_column="coverage",
            max_coverage_column="max_coverage",
            cache_path=str(tmp_path),
        )
        sd = sv._prepare_vue_data({"proteinIndex": pid})["sequenceData"]
        assert len(sd["sequence"]) > 0


def test_combined_spectrum_overlay(fake_fm, monkeypatch, tmp_path):
    """Augmented spectrum: deconv primary + annotated overlay, both present."""
    import streamlit as st

    monkeypatch.setattr(st, "session_state", {}, raising=False)

    from openms_insight import LinePlot

    from src.parse.long_format import explode_combined_spectrum_long

    per_scan = pl.scan_parquet(_DS / "combined_spectrum.pq")
    deconv_long, anno_long = explode_combined_spectrum_long(per_scan)
    dl, al = deconv_long.collect(), anno_long.collect()

    lp = LinePlot(
        cache_id="cs_x",
        data=deconv_long,
        overlay_data=anno_long,
        filters={"deconvIndex": "index"},
        x_column="mass",
        y_column="intensity",
        overlay_x_column="mass",
        overlay_y_column="intensity",
        cache_path=str(tmp_path),
    )
    vd = lp._prepare_vue_data({"deconvIndex": 0})
    assert len(vd["plotData"]) == dl.filter(pl.col("index") == 0).height
    assert len(vd.get("plotDataOverlay", [])) == al.filter(pl.col("index") == 0).height
    assert lp._get_component_args().get("has_overlay") is True


def test_tagger_overlay_data(fake_fm):
    """Tagger plumbing: per-peak signal_peaks triplets + tagData from a tag row."""
    from src.render_oi.tnt_viewer import (
        _deconv_signal_peaks_long,
        _load_polars,
        _tag_data,
    )

    per_scan = _load_polars(fake_fm, "ds", "combined_spectrum")
    dl = _deconv_signal_peaks_long(per_scan).collect()
    assert {"index", "mass", "intensity", "signal_peaks"} <= set(dl.columns)
    sp = dl.row(0, named=True)["signal_peaks"]
    assert sp and len(sp[0]) == 3  # [mz, intensity, charge] (binIdx dropped)

    tid = pl.scan_parquet(_DS / "tag_dfs.pq").collect().row(0, named=True)["TagIndex"]
    td = _tag_data(fake_fm, "ds", tid)
    assert set(td) == {
        "masses",
        "sequence",
        "nTerminal",
        "startPos",
        "endPos",
        "selectedAA",
    }
    assert isinstance(td["masses"], list) and isinstance(td["sequence"], str)
