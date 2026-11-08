"""The comparison views: their order, the QC of the comparison's kind only,
Identification locked for runs that were not identified, the default, and the
layout saved before the view bar."""

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


LABELS = ['Quantification', 'Identification', 'QC: label-free', 'QC: TMT, iTRAQ']


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
    select_view(path, quant_view.IDENTIFICATION_VIEW, QUANT_COMPONENTS, quant_view.presets(),
                quant_view.presets()[0]['layout'])
    assert quant_view.view_rows(fm, tmp_path, 'comparison_lfq') == by_name[quant_view.IDENTIFICATION_VIEW]
    assert quant_view.view_rows(fm, tmp_path, 'comparison_tmt') == by_name[quant_view.IDENTIFICATION_VIEW]


def test_each_comparison_offers_the_qc_of_its_kind(tmp_path):
    fm = _Results(kinds={'comparison_tmt': 'isobaric'})
    assert [p['label'] for p in quant_view.comparison_presets(fm, 'comparison_lfq')] == \
        ['Quantification', 'Identification', 'QC: label-free']
    assert [p['label'] for p in quant_view.comparison_presets(fm, 'comparison_tmt')] == \
        ['Quantification', 'Identification', 'QC: TMT, iTRAQ']
    # The QC picked for a TMT comparison is not shown for a label-free one
    path = quant_view.views_path(tmp_path)
    by_name = {p['name']: p['layout'] for p in quant_view.presets()}
    select_view(path, quant_view.ISOBARIC_QC_VIEW, QUANT_COMPONENTS, quant_view.presets(),
                quant_view.presets()[0]['layout'])
    assert quant_view.view_rows(fm, tmp_path, 'comparison_tmt') == by_name[quant_view.ISOBARIC_QC_VIEW]
    assert quant_view.view_rows(fm, tmp_path, 'comparison_lfq') == by_name[quant_view.QUANTIFICATION_VIEW]


def test_identification_needs_identified_runs():
    from src.layout_manager import missing_analyses
    from src.render.layout_editor.catalog import IDENTIFICATION, QUANT_LABEL_FREE_QC_LAYOUT

    by_name = {p['name']: p['layout'] for p in quant_view.presets()}
    assert missing_analyses(by_name[quant_view.IDENTIFICATION_VIEW], QUANT_COMPONENTS, set()) == [IDENTIFICATION]
    assert missing_analyses(by_name[quant_view.IDENTIFICATION_VIEW], QUANT_COMPONENTS, {IDENTIFICATION}) == []
    # The QC and the quantification views count features without identifications
    for name in (quant_view.QUANTIFICATION_VIEW, quant_view.LABEL_FREE_QC_VIEW, quant_view.ISOBARIC_QC_VIEW):
        assert missing_analyses(by_name[name], QUANT_COMPONENTS, set()) == [], name
    assert QUANT_LABEL_FREE_QC_LAYOUT == by_name[quant_view.LABEL_FREE_QC_VIEW]


class _Manifests:
    def __init__(self, manifests):
        self.manifests = manifests

    def get_results(self, dataset_id, names):
        return {'quant_manifest': self.manifests[dataset_id]}


def test_a_comparison_has_identification_when_its_runs_were_identified():
    from src.render.layout_editor.catalog import IDENTIFICATION

    fm = _Manifests({'a': {'tool': 'flashtnt'}, 'b': {'tool': 'flashdeconv'}, 'c': {}})
    assert quant_view.available(fm, 'a') == {IDENTIFICATION}
    assert quant_view.available(fm, 'b') == set()
    assert quant_view.available(fm, 'c') == set()


def test_the_runs_of_an_isobaric_comparison_are_its_fractions():
    fm = _Manifests({
        'lfq': {'runs': ['run1', 'run2']},
        'tmt': {'kind': 'isobaric', 'runs': ['sample1', 'sample2'], 'fractions': ['run1']},
    })
    assert quant_view.comparison_runs(fm, 'lfq') == ['run1', 'run2']
    assert quant_view.comparison_runs(fm, 'tmt') == ['run1']


def test_saved_views_get_quantification_first_once(tmp_path):
    path = quant_view.views_path(tmp_path)
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps({'current': 'qc-label-free', 'views': [
        {'name': 'identification', 'in_bar': True}, {'name': 'quantification', 'in_bar': True},
        {'name': 'custom-mine', 'label': 'Mine', 'rows': [['quant_table']], 'in_bar': True},
        {'name': 'qc-label-free', 'in_bar': False}, {'name': 'qc-isobaric', 'in_bar': True},
    ]}))
    quant_view.migrate_views(path)
    data = json.loads(path.read_text())
    assert [v['name'] for v in data['views']] == [
        'quantification', 'identification', 'custom-mine', 'qc-label-free', 'qc-isobaric']
    assert data['views'][3]['in_bar'] is False
    assert data['current'] == 'qc-label-free'
    # Once: an order the user gives it afterwards stays
    data['views'][0], data['views'][1] = data['views'][1], data['views'][0]
    path.write_text(json.dumps(data))
    quant_view.migrate_views(path)
    assert [v['name'] for v in json.loads(path.read_text())['views']][:2] == ['identification', 'quantification']


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


def test_a_qc_of_the_other_kind_comes_back_for_its_kind(tmp_path):
    fm = _Results(kinds={'comparison_tmt': 'isobaric'})
    path = quant_view.views_path(tmp_path)
    by_name = {p['name']: p['layout'] for p in quant_view.presets()}
    select_view(path, quant_view.ISOBARIC_QC_VIEW, QUANT_COMPONENTS, quant_view.presets(),
                quant_view.presets()[0]['layout'])
    assert quant_view.view_rows(fm, tmp_path, 'comparison_lfq') == by_name[quant_view.QUANTIFICATION_VIEW]
    assert quant_view.view_rows(fm, tmp_path, 'comparison_tmt') == by_name[quant_view.ISOBARIC_QC_VIEW]
    # No view of one's own was made of it
    names = [v['name'] for v in load_views(path, QUANT_COMPONENTS, quant_view.presets(),
                                           quant_view.presets()[0]['layout'])['views']]
    assert not any(n.startswith('custom-') for n in names)
    # A view picked meanwhile stays
    quant_view.view_rows(fm, tmp_path, 'comparison_lfq')
    select_view(path, quant_view.IDENTIFICATION_VIEW, QUANT_COMPONENTS, quant_view.presets(),
                quant_view.presets()[0]['layout'])
    assert quant_view.view_rows(fm, tmp_path, 'comparison_tmt') == by_name[quant_view.IDENTIFICATION_VIEW]
