"""The FLASHIda Method editor's layout and sense checks (src/render/ida_sections.py):
every key of FLASHIda's method schema has one place, optional parts switch on
and off the way FLASHIda reads a method, a key FLASHIda does not read with the
method as it is says why, and the defaults shown are FLASHIda's."""

from collections import Counter

import pytest

from src import flashida_method as fm
from src.render import ida_sections as sections


def _placed():
    return [p for spec in sections.PARTS.values() for paths in spec[3].values() for p in paths]


def test_every_schema_key_has_one_place():
    placed = _placed()
    # Quantify is the Quantification tab's own switch, not a field
    assert [k for k in fm.keys() if k not in placed] == ['quantification.enabled']
    assert [p for p in placed if p not in fm.keys()] == []
    assert [p for p, n in Counter(placed).items() if n > 1] == []


def test_every_part_is_in_a_tab():
    in_tabs = {part for _, parts in sections.TAB_PARTS.values() for part, _ in parts}
    assert in_tabs == set(sections.PARTS)
    assert list(sections.TAB_PARTS)[0] == 'Global'


def test_names():
    assert sections.name('deconvolution.min_charge') == 'Lowest charge'
    assert sections.name('ms_settings.ms2_quant.agc_target') == 'MS2 quant: AGC target'
    assert sections.name('ms_settings.ms2_quant.agc_target', short=True) == 'AGC target'
    assert sections.name('identification.characterization.exploration.ce_min') == 'MS3 sweep: Lowest energy'


@pytest.mark.parametrize('path, method, active', [
    ('quantification.labelling', {}, False),
    ('quantification.labelling', {'quantification': {'enabled': True}}, True),
    ('quantification.ms3_targets', {'quantification': {'enabled': True, 'scan': 'both'}}, True),
    ('quantification.ms3_targets', {'quantification': {'enabled': True}}, False),
    ('identification.refine.fdr', {}, False),
    ('identification.refine.fdr', {'identification': {'method': 'refine'}}, True),
    ('identification.characterization.max_targets',
     {'identification': {'characterization': {'mode': 'coverage'}}}, True),
    ('identification.characterization.exploration.ce_min',
     {'identification': {'characterization': {'mode': 'coverage'}}}, False),
    ('faims.max_cv_skip', {'faims': {'cv_values': [-50, -60]}}, True),
    ('faims.max_cv_skip', {'faims': {'cv_values': [-50]}}, False),
    ('files.inclusion_list', {'precursor_selection': {'targeting': 'inclusion'}}, True),
    ('files.inclusion_list', {}, False),
    # a key that depends on another one of its scan
    ('ms_settings.ms1.scan_rate', {'ms_settings': {'ms1': {'analyzer': 'Orbitrap'}}}, False),
    ('ms_settings.ms1.scan_rate', {'ms_settings': {'ms1': {'analyzer': 'IonTrap'}}}, True),
    ('ms_settings.ms1.resolution', {'ms_settings': {'ms1': {'analyzer': 'IonTrap'}}}, False),
    ('ms_settings.ms2.collision_energy', {'ms_settings': {'ms2': {'activation': 'ETD'}}}, False),
    ('ms_settings.ms2.reaction_time', {'ms_settings': {'ms2': {'activation': 'HCD'}}}, False),
    ('ms_settings.ms2.reaction_time', {'ms_settings': {'ms2': {'activation': 'EThcD'}}}, True),
    ('precursor_selection.max_charge_states', {'precursor_selection': {'precursor_charges': 'single'}}, False),
    ('deconvolution.min_charge', {}, True),
])
def test_inactive(path, method, active):
    assert (sections.inactive(method, path) is None) is active


def test_parts_switch_on_and_off_as_flashida_reads_them():
    method = {}
    assert 'faims' in sections.present(method)  # an omitted cv_values means FAIMS on at -50
    sections.switch(method, 'faims', False)
    assert fm.get(method, 'faims.cv_values') == [] and 'faims' not in sections.present(method)
    sections.switch(method, 'quant', True)
    assert {'quant', 'ms2q'} <= sections.present(method)
    sections.switch(method, 'ident', True)
    sections.switch(method, 'ms3plan', True)
    assert {'ident', 'ms3plan', 'ms3'} <= sections.present(method)
    sections.switch(method, 'ms3plan', False)  # MS3 is off when its block is gone
    assert 'characterization' not in method['identification']
    assert sections.why_not_add('sweep3', sections.present(method))  # needs MS3


def test_hints_name_what_does_not_fit():
    method = {'identification': {'method': 'flashtnt'},
              'ms_settings': {'ms2': {'activation': 'ETD', 'reaction_time': 0}}}
    found = sections.hints(method)
    assert any('protein database' in h for h in found['ident'])
    assert any('Reaction time' in h for h in found['ms2'])


def test_defaults_are_flashidas():
    assert sections.default_value('ms_settings.ms1.resolution') == 240000
    assert sections.default_of({}, 'ms_settings.ms1.resolution') == '240 000'
    # an MSn scan takes MS1's value where the schema says so
    assert sections.default_of({}, 'ms_settings.ms2.rf_lens') == '60 (from MS1)'
    assert sections.default_of({'ms_settings': {'ms1': {'rf_lens': 45}}}, 'ms_settings.ms2.rf_lens') == '45 (from MS1)'
    assert sections.default_of({}, 'precursor_selection.rt_window') == '180 s'


def test_fixed_values_are_offered():
    assert sections.choices('ms_settings.ms3.activation') == sections.ACTIVATIONS
    assert sections.choices('ms_settings.ms1.data_type') == ['Centroid', 'Profile']
    assert 240000 in sections.choices('ms_settings.ms2.resolution')
    assert sections.choices('precursor_selection.exploration.activations') == sections.ACTIVATIONS
    method = {'quantification': {'conditions': [{'name': 'ctrl'}, {'name': 'drug'}]}}
    assert sections.choices('quantification.enriched_in', method) == ['either', 'ctrl', 'drug']
    assert sections.choices('ms_settings.ms3.analyzer') is None  # the schema lists those itself


def test_summary_lines():
    method = fm.builtin_presets()['LC-MS real-time ID, MS2 only']
    assert sections.line(method, 'ms2').startswith('HCD 35 NCE · Orbitrap')
    assert 'no reference' in sections.line(method, 'ident')


@pytest.mark.parametrize('method, scans', [
    ({}, []),
    ({'ms_settings': {'ms2': {'activation': 'HCD'}}}, []),
    ({'ms_settings': {'ms2': {'activation': 'ETD'}}}, ['MS2']),
    ({'ms_settings': {'ms2': {'activation': 'EThcD'}}}, ['MS2']),
    # a scan FLASHIda does not send with the method as it is does not count
    ({'ms_settings': {'ms3': {'activation': 'ETD'}}}, []),
    ({'ms_settings': {'ms3': {'activation': 'ETD'}},
      'identification': {'characterization': {'mode': 'coverage'}}}, ['MS3']),
    ({'ms_settings': {'ms2_quant': {'activation': 'ETciD'}}}, []),
    ({'ms_settings': {'ms2_quant': {'activation': 'ETciD'}}, 'quantification': {'enabled': True}},
     ['Quant MS2']),
    # an energy sweep over an ETD-family activation, scored, is sent too
    ({'precursor_selection': {'exploration': {'metric': 'remaining_precursor',
                                              'activations': ['HCD', 'EThcD']}}}, ['MS2 sweep']),
    ({'precursor_selection': {'exploration': {'metric': 'none', 'activations': ['ETD']}}}, []),
])
def test_etd_scans(method, scans):
    assert sections.etd_scans(method) == scans
    assert sections.etd_scans(None) == []
