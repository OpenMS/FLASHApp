"""The comparison views: Label-free and TMT, iTRAQ, the default by the
kind of comparison, and the layout saved before the view bar."""

import json

from src.layout_manager import load_views, select_view, validate_layout
from src.render import quant_view
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
    path = quant_view.views_path(tmp_path)
    select_view(path, quant_view.LABEL_FREE_VIEW, QUANT_COMPONENTS, quant_view.presets(),
                quant_view.presets()[0]['layout'])
    assert quant_view.view_rows(fm, tmp_path, 'comparison_tmt') == quant_view.presets()[0]['layout']


def test_former_layout_becomes_my_comparison_layout(tmp_path):
    own = [['quant_table'], ['quant_volcano']]
    fm = _Results(layout=own, kinds={'comparison_tmt': 'isobaric'})
    quant_view.adopt_saved_layout(fm, tmp_path)
    path = quant_view.views_path(tmp_path)
    state = load_views(path, QUANT_COMPONENTS, quant_view.presets(), quant_view.presets()[0]['layout'])
    labels = [v['label'] for v in state['views']]
    assert labels == ['Label-free', 'TMT, iTRAQ', 'My comparison layout']
    # It stays the default of every comparison until a view is picked
    assert quant_view.view_rows(fm, tmp_path, 'comparison_tmt') == own
    assert quant_view.view_rows(fm, tmp_path, 'comparison_lfq') == own
    # Adopted once
    path.write_text(json.dumps({'views': []}))
    quant_view.adopt_saved_layout(fm, tmp_path)
    assert json.loads(path.read_text()) == {'views': []}


def test_former_layout_equal_to_a_preset_adds_no_view(tmp_path):
    fm = _Results(layout=quant_view.presets()[1]['layout'])
    quant_view.adopt_saved_layout(fm, tmp_path)
    assert not quant_view.views_path(tmp_path).exists()
    assert quant_view.default_view(fm, 'comparison_lfq') == quant_view.ISOBARIC_VIEW
