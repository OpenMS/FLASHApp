"""The comparison and GO plots return tables for an openms-insight SeriesPlot."""
import pandas as pd
import polars as pl

from src.render import go_view
from src.render.series import _fingerprint, thin


def test_thin_keeps_at_most_points_per_group_and_the_extremes():
    frame = pl.DataFrame({'g': ['a'] * 1000 + ['b'] * 5, 'v': [float(i) for i in range(1000)] + [1.0, 2.0, None, float('nan'), 3.0]})
    out = thin(frame, ['g'], 'v', points=11)
    a = out.filter(pl.col('g') == 'a')['v']
    assert a.len() == 11 and a.min() == 0 and a.max() == 999
    assert out.filter(pl.col('g') == 'b')['v'].to_list() == [1.0, 2.0, 3.0]


def test_fingerprint_changes_with_data_and_settings():
    frame = pl.DataFrame({'x': [1, 2], 'y': [3, 4]})
    assert _fingerprint(frame, {'mode': 'bars'}) == _fingerprint(frame.clone(), {'mode': 'bars'})
    assert _fingerprint(frame, {'mode': 'bars'}) != _fingerprint(frame, {'mode': 'lines'})
    assert _fingerprint(frame, {}) != _fingerprint(frame.with_columns(y=pl.lit(5)), {})


def test_go_bar_chart_puts_the_most_significant_term_first():
    table = pd.DataFrame({
        'GO_Term': ['late', 'early'], 'p_adj': [0.04, 1e-5], 'Count': [2, 5],
        'GeneRatio': ['2/10', '5/10'], 'FoldEnrichment': [1.5, 3.0],
    })
    spec = go_view._bar_chart(table, 'BP', 10)
    assert spec['data']['term'].to_list() == ['early', 'late']
    assert spec['data']['score'][0] == 5.0
    assert spec['mode'] == 'bars' and spec['orientation'] == 'h'


def test_go_count_chart_reports_the_share_of_annotated_proteins():
    table = pl.DataFrame({'term': ['a', 'b'], 'count': [4, 1]})
    spec = go_view._count_chart(table, 'CC', 1, 8)
    assert spec['data']['proteins'].to_list() == [4]
    assert '50.0% of annotated' in spec['data']['hover'][0]
