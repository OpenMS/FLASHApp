"""
Tests for downsample_heatmap, the heatmap point-reduction helper.

The deconvolved-heatmap panels are built per MS level, and a level with no
peaks (e.g. an MS1-only run that still shows an MS2 panel) yields an empty
frame. That empty frame used to reach scipy's binned_statistic_2d, which
raises "zero-size array to reduction operation minimum which has no identity"
while computing bin edges. downsample_heatmap must short-circuit on empty
input and return an empty, schema-preserving frame instead of crashing.

The helper is pure polars/numpy/scipy (no Streamlit), so it is unit-testable
without booting the app.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import polars as pl

from src.render.compression import downsample_heatmap


HEATMAP_SCHEMA = {"mass": pl.Float64, "rt": pl.Float64, "intensity": pl.Float64}


def test_empty_dataframe_returns_empty_with_schema():
    result = downsample_heatmap(pl.DataFrame(schema=HEATMAP_SCHEMA)).collect()
    assert result.is_empty()
    assert set(result.columns) >= {"mass", "rt", "intensity"}


def test_empty_lazyframe_returns_empty_with_schema():
    result = downsample_heatmap(pl.LazyFrame(schema=HEATMAP_SCHEMA)).collect()
    assert result.is_empty()
    assert set(result.columns) >= {"mass", "rt", "intensity"}


def test_nonempty_input_passes_through_binning():
    df = pl.DataFrame(
        {
            "mass": [100.0, 200.0, 300.0],
            "rt": [1.0, 2.0, 3.0],
            "intensity": [10.0, 20.0, 30.0],
        }
    )
    result = downsample_heatmap(df).collect()
    assert set(result.columns) >= {"mass", "rt", "intensity"}
    assert result.height <= df.height
