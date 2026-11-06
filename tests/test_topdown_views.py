"""The Viewer's views: the ready-made ones fit the catalog, a former layout is adopted."""

from src.layout_manager import load_views, save_views, select_view, validate_layout
from src.render import topdown_layout


class _Results:
    """The part of FileManager the views read: a saved layout and the sequence."""

    def __init__(self, layout=None, sequence=False):
        self.layout = layout
        self.sequence = sequence

    def result_exists(self, dataset_id, name):
        if name == 'sequence':
            return self.sequence
        return dataset_id == topdown_layout.LAYOUT_DATASET and self.layout is not None

    def get_results(self, dataset_id, name):
        return {'layout': self.layout}


def test_presets_are_valid_views():
    for sequence in (False, True):
        fm = _Results(sequence=sequence)
        components = topdown_layout.components(fm)
        for preset in topdown_layout.presets(fm):
            assert validate_layout([preset['layout']], components) == '', preset['name']


def test_default_view_follows_the_run():
    assert topdown_layout.default_view(False, None) == 'deconvolution'
    assert topdown_layout.default_view(True, 'Native') == 'native'
    assert topdown_layout.default_view(True, 'Lysate') == 'lysate'


def test_no_view_until_one_is_picked(tmp_path):
    fm = _Results()
    assert topdown_layout.view_rows(fm, tmp_path, True, 'Native') is None
    assert topdown_layout.experiments_shown(tmp_path) is None
    path = topdown_layout.views_path(tmp_path)
    assert select_view(path, 'lysate', topdown_layout.components(fm), topdown_layout.presets(fm),
                       topdown_layout.deconvolution_default(fm))
    assert topdown_layout.view_rows(fm, tmp_path, True, 'Native') == topdown_layout.TNT_DEFAULT_LAYOUTS['Lysate']
    # A labelled run gets its reporter ions in a ready-made view
    rows = topdown_layout.view_rows(fm, tmp_path, True, 'Native', labelled=True)
    assert list(topdown_layout.REPORTER_ROW) in rows
    assert topdown_layout.experiments_shown(tmp_path) == (1, False)


def test_former_layout_is_adopted(tmp_path):
    own = [[['scan_table', 'mass_table']], [['scan_table']]]
    fm = _Results(layout={'layout': own, 'side_by_side': True})
    topdown_layout.adopt_saved_layout(fm, tmp_path)
    state = load_views(topdown_layout.views_path(tmp_path), topdown_layout.components(fm),
                       topdown_layout.presets(fm), topdown_layout.deconvolution_default(fm))
    mine = next(v for v in state['views'] if v['name'] == state['current'])
    assert mine['label'] == 'My layout' and mine['rows'] == [['scan_table', 'mass_table']]
    assert topdown_layout.experiments_shown(tmp_path) == (2, True)
    # Adopted once: a later save is not overwritten
    state['current'] = 'deconvolution'
    save_views(topdown_layout.views_path(tmp_path), state)
    topdown_layout.adopt_saved_layout(fm, tmp_path)
    assert topdown_layout.view_rows(fm, tmp_path, False, None) == topdown_layout.DECONV_DEFAULT_LAYOUT
