"""Regression test for the FLASHTnT Sequence View data loader.

A fresh FLASHTnT run stores ``sequence_data`` as a *parquet dataset*. The real
``FileManager.get_results`` returns a pandas DataFrame for such a column *unless*
``use_pyarrow=True`` is requested. ``_sequence_table`` must request pyarrow so
``reconstruct_all`` can read the dataset; otherwise it raises, gets swallowed, and
the Sequence View renders "Component unavailable: sequence_view".

This test uses a FileManager fake that mirrors the real loader semantics
(pandas by default, pyarrow Dataset on ``use_pyarrow=True``), so it fails on the
pre-fix code and passes on the fix — independent of the bundled example data.
"""

import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq
import pytest

pl = pytest.importorskip("polars")
pytest.importorskip("openms_insight")


def _write_sequence_data(path):
    fb = [[[float(i)] for i in range(8)], [[float(i)] for i in range(8)]]
    table = pa.table(
        {
            "proteoform_index": pa.array([0, 1], type=pa.int64()),
            "sequence": pa.array(
                [list("PEPTIDER"), list("ACDEFGHK")], type=pa.list_(pa.string())
            ),
            "coverage": pa.array(
                [[0.0, 1.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0], [0.0] * 8],
                type=pa.list_(pa.float64()),
            ),
            "maxCoverage": pa.array([2.0, 0.0], type=pa.float64()),
            "theoretical_mass": pa.array([1000.0, 2000.0], type=pa.float64()),
            # proteoform 1 carries FLASHTnT's -1.0 "unmatched" sentinel.
            "computed_mass": pa.array([1001.5, -1.0], type=pa.float64()),
            "fragment_masses_b": pa.array(
                fb, type=pa.list_(pa.list_(pa.float64()))
            ),
        }
    )
    pq.write_table(table, path)


class _RealisticFM:
    """Mirrors FileManager.get_results for a stored_data ``.pq`` column: pandas by
    default, a pyarrow Dataset when ``use_pyarrow=True`` (see
    src/workflow/FileManager.py:get_results)."""

    def __init__(self, pq_path):
        self._pq = pq_path

    def get_results(
        self, dataset_id, names, partial=False, use_pyarrow=False, use_polars=False
    ):
        out = {}
        for n in names:
            if n != "sequence_data":
                continue
            if use_pyarrow:
                out[n] = ds.dataset(self._pq, format="parquet")
            elif use_polars:
                out[n] = pl.scan_parquet(self._pq)
            else:
                out[n] = pd.read_parquet(self._pq)  # the trap the bug fell into
        return out


def test_sequence_table_reads_parquet_dataset(tmp_path):
    from src.render_oi.tnt_viewer import _sequence_table

    pqp = tmp_path / "sequence_data.pq"
    _write_sequence_data(pqp)
    fm = _RealisticFM(pqp)

    # Sanity: a default get_results would hand back a pandas DataFrame (the form
    # that broke reconstruct_all). The loader must avoid this by using pyarrow.
    assert isinstance(
        fm.get_results("ds", ["sequence_data"])["sequence_data"], pd.DataFrame
    )

    lf = _sequence_table(fm, "ds")
    assert lf is not None, "parquet sequence_data must load (regression: use_pyarrow)"
    df = lf.collect()
    assert df.height == 2
    assert set(df["proteoform_index"].to_list()) == {0, 1}

    row0 = df.filter(pl.col("proteoform_index") == 0).row(0, named=True)
    assert row0["sequence"] == "PEPTIDER"
    assert row0["coverage"][:3] == [0.0, 1.0, 2.0]
    assert row0["max_coverage"] == 2.0

    # Header masses + precomputed fragments flow through; the -1.0 sentinel on
    # proteoform 1 maps to a null observed mass so the header omits it.
    assert row0["theoretical_mass"] == 1000.0
    assert row0["observed_mass"] == 1001.5
    assert row0["fragment_masses_b"][0] == [0.0]
    row1 = df.filter(pl.col("proteoform_index") == 1).row(0, named=True)
    assert row1["observed_mass"] is None
