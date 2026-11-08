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


LABELS = ['Identification', 'Quantification', 'QC: label-free', 'QC: TMT, iTRAQ']


def test_presets_are_valid_views():
    for preset in quant_view.presets():
        assert validate_layout([preset['layout']], QUANT_COMPONENTS) == '', preset['name']
    assert [p['label'] for p in quant_view.presets()] == LABELS


def test_every_component_is_in_a_ready_made_view():
    shown = {n for p in quant_view.presets() for row in p['layout'] for n in row}
    assert shown == {c['name'] for c in QUANT_COMPONENTS}


def test_quantification_is_the_default_of_every_comparison(tmp_path):
    fm = _Results(kinds={'comparison_tmt': 'isobaric'})
    assert quant_view.default_view(fm, 'comparison_lfq') == quant_view.QUANTIFICATION_VIEW
    assert quant_view.default_view(fm, 'comparison_tmt') == quant_view.QUANTIFICATION_VIEW
    by_name = {p['name']: p['layout'] for p in quant_view.presets()}
    assert quant_view.view_rows(fm, tmp_path, 'comparison_tmt') == by_name[quant_view.QUANTIFICATION_VIEW]


def test_a_picked_view_is_shown_for_every_comparison(tmp_path):
    fm = _Results(kinds={'comparison_tmt': 'isobaric'})
    path = quant_view.views_path(tmp_path)
    by_name = {p['name']: p['layout'] for p in quant_view.presets()}
    select_view(path, quant_view.ISOBARIC_QC_VIEW, QUANT_COMPONENTS, quant_view.presets(),
                quant_view.presets()[0]['layout'])
    assert quant_view.view_rows(fm, tmp_path, 'comparison_lfq') == by_name[quant_view.ISOBARIC_QC_VIEW]


def test_former_ready_made_views_become_quantification(tmp_path):
    path = quant_view.views_path(tmp_path)
    path.parent.mkdir(parents=True)
    own = [['quant_table']]
    path.write_text(json.dumps({
        'current': 'isobaric', 'layout': [[['quant_summary']]], 'side_by_side': False,
        'views': [{'name': 'isobaric', 'in_bar': True}, {'name': 'label-free', 'in_bar': False},
                  {'name': 'custom-mine', 'label': 'Mine', 'rows': own, 'in_bar': True}],
    }))
    rows = quant_view.view_rows(_Results(), tmp_path, 'comparison_lfq')
    assert rows == quant_view.QUANT_QUANTIFICATION_LAYOUT
    data = json.loads(path.read_text())
    assert data['current'] == quant_view.QUANTIFICATION_VIEW
    assert [v['name'] for v in data['views']] == [quant_view.QUANTIFICATION_VIEW, 'custom-mine']
    state = load_views(path, QUANT_COMPONENTS, quant_view.presets(), quant_view.presets()[0]['layout'])
    labels = [v['label'] for v in state['views']]
    assert labels[0] == 'Quantification' and 'Mine' in labels and len(labels) == 5


def test_former_layout_becomes_my_comparison_layout(tmp_path):
    own = [['quant_table'], ['quant_volcano']]
    fm = _Results(layout=own, kinds={'comparison_tmt': 'isobaric'})
    quant_view.adopt_saved_layout(fm, tmp_path)
    path = quant_view.views_path(tmp_path)
    state = load_views(path, QUANT_COMPONENTS, quant_view.presets(), quant_view.presets()[0]['layout'])
    labels = [v['label'] for v in state['views']]
    assert labels == LABELS + ['My comparison layout']
    # It stays the default of every comparison until a view is picked
    assert quant_view.view_rows(fm, tmp_path, 'comparison_tmt') == own
    assert quant_view.view_rows(fm, tmp_path, 'comparison_lfq') == own
    # Adopted once
    path.write_text(json.dumps({'views': []}))
    quant_view.adopt_saved_layout(fm, tmp_path)
    assert json.loads(path.read_text()) == {'views': []}


def test_former_layout_equal_to_a_preset_adds_no_view(tmp_path):
    fm = _Results(layout=quant_view.presets()[3]['layout'])
    quant_view.adopt_saved_layout(fm, tmp_path)
    assert not quant_view.views_path(tmp_path).exists()
    assert quant_view.default_view(fm, 'comparison_lfq') == quant_view.ISOBARIC_QC_VIEW


def test_diagnostics_are_four_components():
    names = {c['name'] for c in QUANT_COMPONENTS}
    assert 'quant_diagnostics' not in names
    split = {n for row in quant_view.QUANT_DIAGNOSTICS for n in row}
    assert split == {'quant_missing_values', 'quant_imputed_values', 'quant_normalization', 'quant_pvalues'}
    assert split <= names
    from src.render import quant_cells
    assert names == set(quant_cells.BUILDERS)


def test_split_diagnostics_keeps_the_rest_of_its_row():
    rows = [['quant_table'], ['quant_pca', 'quant_diagnostics'], ['quant_go']]
    assert quant_view.split_diagnostics(rows) == [
        ['quant_table'], ['quant_pca'],
        ['quant_missing_values', 'quant_imputed_values'], ['quant_normalization', 'quant_pvalues'],
        ['quant_go'],
    ]


def test_saved_views_with_diagnostics_show_the_four_plots(tmp_path):
    path = quant_view.views_path(tmp_path)
    path.parent.mkdir(parents=True)
    own = [['quant_table'], ['quant_diagnostics']]
    path.write_text(json.dumps({
        'current': 'mine', 'layout': [own], 'side_by_side': False,
        'views': [{'name': 'mine', 'label': 'Mine', 'rows': own, 'in_bar': True}],
    }))
    fm = _Results()
    rows = quant_view.view_rows(fm, tmp_path, 'comparison_lfq')
    assert rows == [['quant_table'], *quant_view.QUANT_DIAGNOSTICS]
    assert 'quant_diagnostics' not in path.read_text()


def test_former_layout_with_diagnostics_is_split(tmp_path):
    fm = _Results(layout=[['quant_table'], ['quant_diagnostics']])
    assert quant_view.layout(fm) == [['quant_table'], *quant_view.QUANT_DIAGNOSTICS]
