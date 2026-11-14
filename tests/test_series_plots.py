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


def test_go_bar_chart_puts_the_most_significant_term_first_per_category():
    tables = {
        'BP': pd.DataFrame({
            'GO_Term': ['late', 'early'], 'p_adj': [0.04, 1e-5], 'Count': [2, 5],
            'GeneRatio': ['2/10', '5/10'], 'FoldEnrichment': [1.5, 3.0],
        }),
        'MF': pd.DataFrame({
            'GO_Term': ['binding'], 'p_adj': [0.01], 'Count': [3],
            'GeneRatio': ['3/10'], 'FoldEnrichment': [2.0],
        }),
        'CC': None,
    }
    spec = go_view._bar_chart(tables)
    data = spec['data']
    assert data.filter(pl.col('category') == 'BP')['term'].to_list() == ['early', 'late']
    assert data.filter(pl.col('category') == 'MF')['term'].to_list() == ['binding']
    assert data['score'][0] == 5.0
    assert spec['mode'] == 'bars' and spec['orientation'] == 'h'
    # The chart picks the category and the number of terms itself
    assert spec['view_column'] == 'category' and spec['view_identifier'] == go_view.CATEGORY_IDENTIFIER
    assert spec['top_choices'] == go_view.TERM_COUNTS and spec['top_default'] == go_view.TOP_TERMS
    assert (spec['x_column'], spec['y_column']) == ('term', 'score')


def test_go_count_chart_reports_the_share_of_annotated_proteins():
    tables = {'CC': pl.DataFrame({'term': ['a', 'b'], 'count': [4, 1]}), 'BP': None}
    spec = go_view._count_chart(tables, 8)
    assert spec['data']['proteins'].to_list() == [4, 1]
    assert spec['data']['category'].to_list() == ['CC', 'CC']
    assert '50.0% of annotated' in spec['data']['hover'][0]


def test_go_term_table_holds_every_category():
    frames = {'BP': pl.DataFrame({'term': ['a']}), 'CC': pl.DataFrame({'term': ['b', 'c']})}
    table = go_view._all_categories(frames)
    assert table['category'].to_list() == ['BP', 'CC', 'CC']
    assert table['id'].to_list() == [0, 1, 2]


def test_row_height_lines_up_the_views_of_a_row():
    from src.render import quant_cells

    # The tallest view of the row sets it; views that size themselves do not count
    assert quant_cells.row_height(['quant_feature_map', 'quant_feature_intensity']) == 480
    assert quant_cells.row_height(['quant_correlation', 'quant_cv']) == 480
    assert quant_cells.row_height(['quant_totals', 'quant_go']) == quant_cells.DEFAULT_HEIGHT
    # A view alone in its row keeps its own height
    assert quant_cells.row_height(['quant_pvalues']) is None
    # The normalization view grows with the samples, and its partner with it
    assert quant_cells.row_height(['quant_totals', 'quant_normalization'], samples=20) == 44 * 20 + 140


def test_supports_falls_back_for_a_series_plot_without_bars():
    from src.render.series import supports

    class Old:
        def __init__(self, cache_id, data=None, x_column='x', y_column='y', mode='lines+markers', title=None):
            pass

    class New(Old):
        def __init__(self, cache_id, data=None, x_column='x', y_column='y', mode='lines+markers', title=None,
                     orientation='v', view_column=None, **kwargs):
            pass

    assert supports(Old, {'mode': 'lines', 'title': 't'})
    assert not supports(Old, {'mode': 'bars'})
    assert not supports(Old, {'mode': 'lines', 'view_column': 'c'})
    assert supports(New, {'mode': 'bars', 'orientation': 'h', 'view_column': 'c'})


def test_figure_draws_grouped_and_horizontal_bars():
    from src.render.series import figure

    frame = pl.DataFrame({'x': ['a', 'b', 'a', 'b'], 'y': [1, 2, 3, 4], 's': ['p', 'p', 'q', 'q']})
    fig = figure(frame, series_column='s', series_colors={'p': '#ff0000'}, mode='bars')
    assert [t.type for t in fig.data] == ['bar', 'bar'] and fig.layout.barmode == 'group'
    assert fig.data[0].marker.color == '#ff0000'
    fig = figure(frame.filter(pl.col('s') == 'p'), mode='bars', orientation='h', x_label='GO term')
    assert fig.data[0].orientation == 'h' and list(fig.data[0].y) == ['a', 'b']
    assert fig.layout.yaxis.title.text == 'GO term'
    assert figure(frame, series_column='s', mode='box').data[0].type == 'box'
