"""Shared pytest fixtures for the FLASHApp render/schema construct-smoke tests.

Mirrors the OpenMS-Insight ``tests/conftest.py`` ``mock_streamlit`` fixture (patch
``st.session_state`` with a dict so components run without a Streamlit server) and
adds light mocks for the Streamlit *layout* primitives that the frozen
``render_linked_grid`` touches (``st.columns`` / ``st.container`` / ``st.warning``).

These cannot run as a Streamlit ``AppTest`` because OpenMS-Insight's subprocess
(spawn) preprocessing is incompatible with AppTest's runtime; instead the smoke
constructs synthetic FileManager caches, runs ``build_insight_caches`` +
``make_builders``, and exercises each component's ``_prepare_vue_data`` /
``_get_component_args`` over its on-disk ``data_path=`` cache.
"""

from __future__ import annotations

import sys
import tempfile
import shutil
from contextlib import contextmanager
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import polars as pl
import pyarrow.parquet as pq
import pytest

# Ensure the FLASHApp repo root is importable (``src`` package).
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))


class MockSessionState(dict):
    """Mock Streamlit session_state that behaves like a dict (attr + item access)."""

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError as exc:  # pragma: no cover - defensive
            raise AttributeError(name) from exc

    def __setattr__(self, name, value):
        self[name] = value


class _MockColumn:
    """Stand-in for a Streamlit column/container: context manager + no-op widgets."""

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def warning(self, *a, **k):
        return None

    def info(self, *a, **k):
        return None

    def container(self, *a, **k):
        return _MockColumn()


@pytest.fixture
def mock_streamlit():
    """Patch ``st.session_state`` + the layout primitives ``render_linked_grid`` uses."""
    session = MockSessionState()

    def _columns(spec, *a, **k):
        n = spec if isinstance(spec, int) else len(spec)
        return [_MockColumn() for _ in range(n)]

    @contextmanager
    def _container(*a, **k):
        yield _MockColumn()

    with patch("streamlit.session_state", session), \
         patch("streamlit.columns", _columns), \
         patch("streamlit.container", lambda *a, **k: _MockColumn()), \
         patch("streamlit.divider", lambda *a, **k: None), \
         patch("streamlit.warning", lambda *a, **k: None):
        yield session


@pytest.fixture
def temp_workspace():
    """A throwaway FLASHApp workspace directory (with its own cache)."""
    tmp = tempfile.mkdtemp(prefix="flashapp_render_test_")
    yield Path(tmp)
    shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------- #
# Synthetic oracle-cache builders (matching the src/parse/* output schemas)
# --------------------------------------------------------------------------- #
def _sp_schema():
    return {
        "index": pl.Int64,
        "MonoMass": pl.List(pl.Float64),
        "SumIntensity": pl.List(pl.Float64),
        "SignalPeaks": pl.List(pl.List(pl.List(pl.Float64))),
        "MonoMass_Anno": pl.List(pl.Float64),
        "SumIntensity_Anno": pl.List(pl.Float64),
    }


def _sn_schema():
    return {
        "index": pl.Int64,
        "PrecursorScan": pl.Float64,
        "SignalPeaks": pl.List(pl.List(pl.List(pl.Float64))),
        "NoisyPeaks": pl.List(pl.List(pl.List(pl.Float64))),
    }


def make_deconv_caches(fm, ds="exp1"):
    """Write a tiny set of FLASHDeconv-style oracle caches (deconv + raw heatmaps)."""
    fm.store_data(ds, "scan_table", pl.DataFrame({
        "index": [0, 1], "Scan": [10, 20], "MSLevel": [1, 2],
        "RT": [1.0, 2.0], "PrecursorMass": [1000.0, 2000.0], "#Masses": [2, 1]}))
    fm.store_data(ds, "mass_table", pl.DataFrame({
        "index": [0, 1],
        "MonoMass": [[100.0, 200.0], [300.0]],
        "SumIntensity": [[1.0, 2.0], [3.0]],
        "MinCharges": [[1, 2], [3]], "MaxCharges": [[2, 3], [4]],
        "MinIsotopes": [[0, 0], [0]], "MaxIsotopes": [[1, 1], [1]],
        "CosineScore": [[0.9, 0.8], [0.7]], "SNR": [[5.0, 4.0], [3.0]],
        "QScore": [[0.99, 0.98], [0.97]]}))
    fm.store_data(ds, "deconv_spectrum", pl.DataFrame({
        "index": [0, 1], "MonoMass": [[100.0, 200.0], [300.0]],
        "SumIntensity": [[1.0, 2.0], [3.0]]}))
    fm.store_data(ds, "combined_spectrum", pl.DataFrame({
        "index": [0, 1],
        "MonoMass": [[100.0, 200.0], [300.0]],
        "SumIntensity": [[1.0, 2.0], [3.0]],
        "SignalPeaks": [
            [[[0.0, 75.0, 3.0, 12.0], [1.0, 75.1, 1.0, 12.0]], [[3.0, 125.0, 4.0, 5.0]]],
            [[[0.0, 150.0, 2.0, 2.0]]],
        ],
        "MonoMass_Anno": [[75.0, 75.1, 125.0, 99.0], [150.0]],
        "SumIntensity_Anno": [[3.0, 1.0, 4.0, 0.5], [2.0]],
    }, schema=_sp_schema()))
    fm.store_data(ds, "threedim_SN_plot", pl.DataFrame({
        "index": [0, 1], "PrecursorScan": [0.0, 0.0],
        "SignalPeaks": [
            [[[0.0, 75.0, 3.0, 12.0], [1.0, 75.1, 1.0, 12.0]], [[3.0, 125.0, 4.0, 5.0]]],
            [[[0.0, 150.0, 2.0, 2.0]]],
        ],
        "NoisyPeaks": [[[[2.0, 80.0, 0.5, 12.0]], []], [[]]],
    }, schema=_sn_schema()))
    # full-resolution heatmaps (already tidy: rt, mass, intensity)
    for tag in ("ms1_deconv_heatmap", "ms2_deconv_heatmap",
                "ms1_raw_heatmap", "ms2_raw_heatmap"):
        fm.store_data(ds, tag, pl.DataFrame({
            "rt": [1.0, 1.0, 2.0, 2.0],
            "mass": [100.0, 200.0, 300.0, 400.0],
            "intensity": [10.0, 20.0, 30.0, 40.0]}))
    fm.store_data(ds, "density_target", pd.DataFrame({"x": [0.1, 0.2], "y": [1.0, 2.0]}))
    fm.store_data(ds, "density_decoy", pd.DataFrame({"x": [0.3, 0.4], "y": [0.5, 0.6]}))
    return ds


def make_sequence_cache(fm):
    """Write the global deconv sequence cache ('sequence','sequence')."""
    fm.store_data("sequence", "sequence", {
        "input_sequence": "PEPTIDEK",
        "fixed_mod_cysteine": False,
        "fixed_mod_methionine": False,
    })


def make_tnt_caches(fm, ds="exp1"):
    """Write FLASHTnT-style oracle caches (proteins, tags, sequence_data, settings)."""
    from src.render.sequence import getFragmentDataFromSeq
    from src.render.sequence_data_store import build_table, ROW_GROUP_SIZE

    make_deconv_caches(fm, ds)  # tnt reuses the deconv-style spectra

    # Mirror the real protein.tsv columns that survive parse/tnt.py's rename
    # (ProteoformIndex->index, ProteinAccession->accession, etc. + added length),
    # including the curated-display fields the oracle ProteinTable shows
    # (MatchingFragments / ModCount / TagCount / Score) so the ported
    # column_definitions + initialSort(Score desc) exercise real columns.
    protein_df = pd.DataFrame({
        "index": [0, 1], "accession": ["P1", "DECOY_P2"],
        "description": ["d1", "d2"], "sequence": ["PEPTIDEK", "ACDEFGHK"],
        "length": [8, 8], "ProteoformMass": [900.4, 800.3],
        "MatchingFragments": [12, 8], "Coverage(%)": [55.0, 40.0],
        "ModCount": [0, 1], "TagCount": [2, 1], "Score": [5.0, 6.0],
        "ProteoformLevelQvalue": [0.01, 0.5], "Scan": [10, 20]})
    fm.store_data(ds, "protein_dfs", protein_df)

    # Mirror the real tags.tsv columns that survive parse/tnt.py's rename
    # (DeNovoScore->Score, Masses->mzs, StartPosition->StartPos + added EndPos),
    # including Nmass / Cmass / DeltaMass the oracle TagTable shows (Nmass/Cmass use
    # the -1->"-" placeholder). -1 in Nmass/Cmass exercises that formatter's data.
    tag_df = pd.DataFrame({
        "Scan": [10, 10, 20], "TagSequence": ["PEP", "TID", "ACD"],
        "StartPos": [0, 3, 0], "EndPos": [2, 5, 2], "Length": [3, 3, 3],
        "Score": [5.0, 4.0, 6.0], "mzs": ["1,2,3", "4,5,6", "7,8,9"],
        "Nmass": [-1.0, 100.5, 200.5], "Cmass": [300.5, -1.0, 400.5],
        "DeltaMass": [0.1, 0.2, 0.3], "ProteinIndex": [0, 0, 1]})
    fm.store_data(ds, "tag_dfs", tag_df, row_group_size=128)

    seqdata = {}
    for pid, seq in [(0, "PEPTIDEK"), (1, "ACDEFGHK")]:
        cov = np.array([1.0] * len(seq))
        entry = getFragmentDataFromSeq(seq, list(cov / cov.max()), cov.max(), [])
        entry["sequence"] = list(seq)
        entry["proteoform_start"] = -1
        entry["proteoform_end"] = -1
        entry["computed_mass"] = 900.0
        entry["theoretical_mass"] = 900.0
        entry["modifications"] = []
        seqdata[pid] = entry
    tbl = build_table(seqdata)
    with fm.parquet_sink(ds, "sequence_data") as p:
        pq.write_table(tbl, p, row_group_size=ROW_GROUP_SIZE)

    fm.store_data(ds, "settings", {"tolerance": 10.0, "ion_types": ["b", "y"]})
    fm.store_data(ds, "density_id_target", pd.DataFrame({"x": [0.1, 0.2], "y": [1.0, 2.0]}))
    fm.store_data(ds, "density_id_decoy", pd.DataFrame(columns=["x", "y"]))
    return ds


def make_quant_caches(fm, ds="exp1"):
    """Write a FLASHQuant-style oracle quant_dfs cache."""
    quant = pd.DataFrame({
        "FeatureGroupIndex": [0, 1],
        "MonoisotopicMass": [1000.0, 2000.0], "AverageMass": [1000.5, 2000.5],
        "StartRetentionTime(FWHM)": [1.0, 3.0], "EndRetentionTime(FWHM)": [2.0, 4.0],
        "HighestApexRetentionTime": [1.5, 3.5], "FeatureGroupQuantity": [100.0, 200.0],
        "AllAreaUnderTheCurve": [150.0, 250.0], "MinCharge": [1, 2], "MaxCharge": [3, 4],
        "MostAbundantFeatureCharge": [2, 3], "IsotopeCosineScore": [0.99, 0.98],
        "Charges": [np.array([2, 3]), np.array([4])],
        "IsotopeIndices": [np.array([0, 1]), np.array([0])],
        "CentroidMzs": [np.array([500.1, 500.2]), np.array([501.0])],
        "RTs": [["1.0,1.5,2.0", "1.1,1.6"], ["3.0,3.5"]],
        "MZs": [["500.10,500.12,500.14", "500.20,500.22"], ["501.00,501.05"]],
        "Intensities": [["10,20,15", "5,8"], ["30,25"]],
    })
    fm.store_data(ds, "quant_dfs", quant)
    return ds
