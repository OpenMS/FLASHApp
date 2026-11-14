"""Draw a small table as an openms-insight ``SeriesPlot`` (lines, bars or boxes).

The plots of a comparison's quality control, of its comparison steps and of
the GO enrichment are computed in the session from tables already in the
cache, so they are cached here by their content and settings: a rerun reuses
the plot, a new table or setting builds it again.
"""

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl

# Rows kept per box: evenly spaced values of the sorted ones, which draw the same box
BOX_POINTS = 201


def thin(frame: pl.DataFrame, by: list[str], value: str, points: int = BOX_POINTS) -> pl.DataFrame:
    """At most ``points`` evenly spaced values of ``value`` per group of ``by``."""
    frame = frame.filter(pl.col(value).is_not_null() & pl.col(value).is_finite())
    parts = []
    for _, part in frame.group_by(by, maintain_order=True):
        if part.height > points:
            rows = np.unique(np.round(np.linspace(0, part.height - 1, points)).astype(int))
            part = part.sort(value)[rows.tolist()]
        parts.append(part)
    return pl.concat(parts) if parts else frame


def _fingerprint(frame: pl.DataFrame, settings: dict) -> str:
    digest = hashlib.sha1(repr(sorted(settings.items())).encode())
    digest.update(repr((frame.columns, frame.shape)).encode())
    if frame.height:
        digest.update(str(int(frame.hash_rows(seed=0).sum())).encode())
    return digest.hexdigest()[:16]


def show_series(data, name: str, cache_path, key: str, height=None, state_manager=None, **plot_kwargs):
    """Draw ``data`` (polars or pandas) as a ``SeriesPlot``; ``plot_kwargs`` go to it."""
    from openms_insight import SeriesPlot

    frame = pl.from_pandas(data) if isinstance(data, pd.DataFrame) else data
    cache_id = f'{name}-{_fingerprint(frame, plot_kwargs)}'
    plot = SeriesPlot(
        cache_id=cache_id, cache_path=str(Path(cache_path, 'insight', 'series')), data=frame.lazy(),
        **plot_kwargs,
    )
    plot(key=key, state_manager=state_manager, height=height or 400)
    return plot
