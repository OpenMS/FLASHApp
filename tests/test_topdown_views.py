"""The Viewer's views: the ready-made ones fit the catalog, a former layout is adopted."""

from src.layout_manager import load_views, save_views, select_view, validate_layout
from src.render import topdown_layout


class _Results:
    """The part of FileManager the views read: saved layouts and the sequence."""

    def __init__(self, layout=None, sequence=False, **legacy):
        # legacy: deconvolution= / identification= layouts the former Viewers saved
        self.layouts = {topdown_layout.LAYOUT_DATASET: layout,
                        topdown_layout.LEGACY_LAYOUT_DATASET[False]: legacy.get('deconvolution'),
                        topdown_layout.LEGACY_LAYOUT_DATASET[True]: legacy.get('identification')}
        self.sequence = sequence

    def result_exists(self, dataset_id, name):
        if name == 'sequence':
            return self.sequence
        return self.layouts.get(dataset_id) is not None

    def get_results(self, dataset_id, name):
        return {'layout': self.layouts[dataset_id]}


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


def _views(fm, tmp_path):
    return load_views(topdown_layout.views_path(tmp_path), topdown_layout.components(fm),
                      topdown_layout.presets(fm), topdown_layout.deconvolution_default(fm))


def test_layouts_of_the_former_viewers_become_views(tmp_path):
    deconv = [['ms1_deconv_heat_map'], ['scan_table', 'deconv_spectrum']]
    ident = [['protein_table', 'tag_table'], ['sequence_view']]
    fm = _Results(deconvolution={'layout': [deconv], 'side_by_side': False},
                  identification={'layout': [ident], 'side_by_side': False})
    topdown_layout.adopt_saved_layout(fm, tmp_path)
    views = {v['label']: v for v in _views(fm, tmp_path)['views']}
    assert set(views) == {'Deconvolution', 'Native', 'Lysate', 'Label-free', 'TMT, iTRAQ',
                          'My deconvolution layout', 'My identification layout'}
    assert views['My deconvolution layout']['rows'] == deconv
    assert views['My identification layout']['rows'] == ident
    # Nothing counts as picked: each run keeps the default of its kind, shown as its view
    assert topdown_layout.view_rows(fm, tmp_path, False, None) is None
    assert topdown_layout.default_view(False, None, fm) == 'custom-deconvolution'
    assert topdown_layout.default_view(True, 'Lysate', fm) == 'custom-identification'
    assert topdown_layout.default_rows(fm, True, 'Lysate') == ident


def test_a_former_layout_equal_to_a_ready_made_view_is_that_view(tmp_path):
    fm = _Results(identification={'layout': [topdown_layout.TNT_DEFAULT_LAYOUTS['Lysate']], 'side_by_side': False})
    topdown_layout.adopt_saved_layout(fm, tmp_path)
    assert not topdown_layout.views_path(tmp_path).exists()
    assert topdown_layout.default_view(True, 'Native', fm) == 'lysate'


def test_every_former_layout_is_kept_next_to_an_edit_layout_one(tmp_path):
    deconv = [['ms1_deconv_heat_map'], ['scan_table']]
    own = [[['scan_table', 'mass_table']], [['mass_table']]]
    fm = _Results(layout={'layout': own, 'side_by_side': True},
                  deconvolution={'layout': [deconv], 'side_by_side': False})
    topdown_layout.adopt_saved_layout(fm, tmp_path)
    state = _views(fm, tmp_path)
    labels = [v['label'] for v in state['views']]
    assert 'My layout' in labels and 'My deconvolution layout' in labels
    assert next(v for v in state['views'] if v['name'] == state['current'])['label'] == 'My layout'
    # The experiments of the Edit layout one are kept as they were
    assert topdown_layout.experiments_shown(tmp_path) == (2, True)


def test_every_view_is_listed_at_once(tmp_path):
    fm = _Results()
    labels = [v['label'] for v in _views(fm, tmp_path)['views']]
    assert labels == ['Deconvolution', 'Native', 'Lysate', 'Label-free', 'TMT, iTRAQ']
