"""The 3D MS1 map is drawn by the MS1 heatmap now: saved layouts and views
that name it show the heatmap instead, once."""

import json

from src.render import topdown_layout
from src.render.layout_editor.catalog import DECONV_COMPONENTS


def test_3d_map_becomes_the_heatmap():
    assert topdown_layout.merge_components([['ms1_3d_map', 'scan_table']]) == [['ms1_raw_heatmap', 'scan_table']]


def test_3d_map_beside_the_heatmap_goes():
    rows = [['ms1_raw_heatmap', 'scan_table'], ['ms1_3d_map']]
    assert topdown_layout.merge_components(rows) == [['ms1_raw_heatmap', 'scan_table']]


def test_other_rows_untouched():
    rows = [['scan_table'], ['mass_table', 'deconv_spectrum']]
    assert topdown_layout.merge_components(rows) == rows


def test_catalog_no_longer_offers_the_3d_map():
    assert 'ms1_3d_map' not in {c['name'] for c in DECONV_COMPONENTS}


def test_saved_views_are_merged(tmp_path):
    path = topdown_layout.views_path(tmp_path)
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps({
        'layout': [[['ms1_3d_map'], ['scan_table']], [['ms1_raw_heatmap', 'ms1_3d_map']]],
        'current': 'custom-1',
        'views': [{'name': 'deconvolution', 'in_bar': True},
                  {'name': 'custom-1', 'label': 'Mine', 'in_bar': True, 'rows': [['ms1_3d_map', 'mass_table']]}],
    }))
    topdown_layout.merge_saved_views(tmp_path)
    data = json.loads(path.read_text())
    assert data['layout'] == [[['ms1_raw_heatmap'], ['scan_table']], [['ms1_raw_heatmap']]]
    assert data['views'][1]['rows'] == [['ms1_raw_heatmap', 'mass_table']]
    assert data['current'] == 'custom-1'


def test_views_without_it_are_not_rewritten(tmp_path):
    path = topdown_layout.views_path(tmp_path)
    path.parent.mkdir(parents=True)
    text = '{"views": [], "layout": [[["scan_table"]]]}'
    path.write_text(text)
    topdown_layout.merge_saved_views(tmp_path)
    assert path.read_text() == text
