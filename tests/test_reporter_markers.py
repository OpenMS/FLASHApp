"""The reporter ion plot marks every channel when openms-insight has x markers."""

import polars as pl
import pytest

pytest.importorskip("openms_insight")

from src.render.render import REPORTER_TOLERANCE, reporter_plot_args  # noqa: E402

CHANNELS = ["126", "127N", "127C"]
MZS = [126.127726, 127.124761, 0.0]
FRAME = pl.DataFrame({
    "scan_idx": [0] * 5,
    "mz": [125.6, 126.127726, 127.124761, 127.131081, 128.6],
    "intensity": [0.0, 5.0, 0.0, 3.0, 0.0],
    "label": ["", "126", "127N", "127C", ""],
})


def test_channels_as_markers():
    data_args, render = reporter_plot_args(CHANNELS, MZS, markers=True)
    # A channel without a reporter m/z is not marked
    assert render == {
        "x_markers": {"126": 126.127726, "127N": 127.124761},
        "x_marker_tolerance": REPORTER_TOLERANCE,
    }
    args = data_args(FRAME)
    assert "annotation_column" not in args
    # The channel without signal is left to its marker; the edge points stay
    assert args["data"]["label"].to_list() == ["", "126", "127C", ""]


def test_older_insight_labels_the_sticks():
    data_args, render = reporter_plot_args(CHANNELS, MZS, markers=False)
    assert render == {}
    args = data_args(FRAME)
    assert args["annotation_column"] == "label"
    assert args["data"].height == FRAME.height
