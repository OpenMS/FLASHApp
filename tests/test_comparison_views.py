"""The comparison views: Label-free and TMT, iTRAQ, the default by the
kind of comparison, and the layout saved before the view bar."""

import json

from src.layout_manager import load_views, select_view, validate_layout
from src.render import quant_view, topdown_layout
from src.render.layout_editor.catalog import QUANT_COMPONENTS


class _Results:
    """The part of FileManager the comparison views read."""

    def __init__(self, layout=None, kinds=None):
        self.layout = layout
        self.kinds = kinds or {}

    def result_exists(self, dataset_id, name):
        return dataset_id == quant_view.LAYOUT_DATASET and self.layout is not None

    def get_results(self, dataset_id, names):
        if dataset_id == quant_view.LAYOUT_DATASET:
            return {'layout': {'layout': [self.layout], 'side_by_side': False}}
        return {'quant_manifest': {'kind': self.kinds.get(dataset_id, 'label-free')}}


def test_presets_are_valid_views():
    for preset in quant_view.presets():
        assert validate_layout([preset['layout']], QUANT_COMPONENTS) == '', preset['name']
    assert [p['label'] for p in quant_view.presets()] == ['Label-free', 'TMT, iTRAQ']


def test_default_view_follows_the_comparison(tmp_path):
    fm = _Results(kinds={'comparison_tmt': 'isobaric'})
    assert quant_view.default_view(fm, 'comparison_lfq') == quant_view.LABEL_FREE_VIEW
    assert quant_view.default_view(fm, 'comparison_tmt') == quant_view.ISOBARIC_VIEW
    by_name = {p['name']: p['layout'] for p in quant_view.presets()}
    assert quant_view.view_rows(fm, tmp_path, 'comparison_tmt') == by_name[quant_view.ISOBARIC_VIEW]
    assert quant_view.view_rows(fm, tmp_path, 'comparison_lfq') == by_name[quant_view.LABEL_FREE_VIEW]


def test_a_picked_view_is_shown_for_every_comparison(tmp_path):
    fm = _Results(kinds={'comparison_tmt': 'isobaric'})
    path = topdown_layout.views_path(tmp_path)
    select_view(path, quant_view.LABEL_FREE_VIEW, topdown_layout.components(fm), topdown_layout.presets(fm),
                topdown_layout.deconvolution_default(fm))
    assert quant_view.view_rows(fm, tmp_path, 'comparison_tmt') == quant_view.presets()[0]['layout']


def test_former_layout_becomes_my_comparison_layout(tmp_path):
    own = [['quant_table'], ['quant_volcano']]
    fm = _Results(layout=own, kinds={'comparison_tmt': 'isobaric'})
    quant_view.adopt_saved_layout(fm, tmp_path)
    path = quant_view.views_path(tmp_path)
    state = load_views(path, QUANT_COMPONENTS, quant_view.presets(), quant_view.presets()[0]['layout'])
    labels = [v['label'] for v in state['views']]
    assert labels == ['Label-free', 'TMT, iTRAQ', 'My comparison layout']
    # Adopted once
    saved = path.read_text()
    quant_view.adopt_saved_layout(fm, tmp_path)
    assert path.read_text() == saved
    # ... and listed with the views of a run
    quant_view.merge_views(tmp_path)
    assert not path.exists()
    state = load_views(topdown_layout.views_path(tmp_path), topdown_layout.components(fm),
                       topdown_layout.presets(fm), topdown_layout.deconvolution_default(fm))
    assert [v['label'] for v in state['views']] == [
        'Deconvolution', 'Native', 'Lysate', 'Label-free', 'TMT, iTRAQ', 'My comparison layout']
    # It stays the default of every comparison until a view is picked
    assert quant_view.view_rows(fm, tmp_path, 'comparison_tmt') == own
    assert quant_view.view_rows(fm, tmp_path, 'comparison_lfq') == own


def test_comparison_views_join_the_viewers_views(tmp_path):
    fm = _Results()
    run_views = topdown_layout.views_path(tmp_path)
    run_views.parent.mkdir(parents=True)
    run_views.write_text(json.dumps({'views': [{'name': 'deconvolution', 'in_bar': True},
                                               {'name': 'custom-1', 'label': 'Mine', 'rows': [['scan_table']]}]}))
    quant_view.views_path(tmp_path).write_text(json.dumps({'views': [
        {'name': 'isobaric', 'in_bar': False},
        {'name': 'custom-1', 'label': 'My volcano', 'rows': [['quant_volcano']]}]}))
    quant_view.merge_views(tmp_path)
    state = load_views(run_views, topdown_layout.components(fm), topdown_layout.presets(fm),
                       topdown_layout.deconvolution_default(fm))
    by_label = {v['label']: v for v in state['views']}
    assert by_label['My volcano']['rows'] == [['quant_volcano']] and by_label['My volcano']['name'] != 'custom-1'
    assert by_label['Mine']['rows'] == [['scan_table']]
    assert by_label['TMT, iTRAQ']['in_bar'] is False
    assert set(by_label) == {'Deconvolution', 'Native', 'Lysate', 'Label-free', 'TMT, iTRAQ', 'Mine', 'My volcano'}
    assert quant_view.views_path(tmp_path).with_name('comparison.json.merged').exists()


def test_former_layout_equal_to_a_preset_adds_no_view(tmp_path):
    fm = _Results(layout=quant_view.presets()[1]['layout'])
    quant_view.adopt_saved_layout(fm, tmp_path)
    assert not quant_view.views_path(tmp_path).exists()
    assert quant_view.default_view(fm, 'comparison_lfq') == quant_view.ISOBARIC_VIEW
