"""The comparison plots return tables for an openms-insight SeriesPlot."""
import polars as pl

from src.render.series import thin
from src.viz.plumbing import fingerprint


def test_thin_keeps_at_most_points_per_group_and_the_extremes():
    frame = pl.DataFrame({'g': ['a'] * 1000 + ['b'] * 5, 'v': [float(i) for i in range(1000)] + [1.0, 2.0, None, float('nan'), 3.0]})
    out = thin(frame, ['g'], 'v', points=11)
    a = out.filter(pl.col('g') == 'a')['v']
    assert a.len() == 11 and a.min() == 0 and a.max() == 999
    assert out.filter(pl.col('g') == 'b')['v'].to_list() == [1.0, 2.0, 3.0]


def test_fingerprint_changes_with_data_and_settings():
    frame = pl.DataFrame({'x': [1, 2], 'y': [3, 4]})
    assert fingerprint(frame, ('mode', 'bars')) == fingerprint(frame.clone(), ('mode', 'bars'))
    assert fingerprint(frame, ('mode', 'bars')) != fingerprint(frame, ('mode', 'lines'))
    assert fingerprint(frame) != fingerprint(frame.with_columns(y=pl.lit(5)))
