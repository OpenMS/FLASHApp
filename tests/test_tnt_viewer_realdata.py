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
