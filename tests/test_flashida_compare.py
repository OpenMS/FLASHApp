"""Two FLASHIda methods compared (src/flashida_compare.py): only settings that
change the acquisition count, a key left out is its default, and settings of A
are taken over into B."""

import copy

import pytest

from src import flashida_compare as fc
from src import flashida_method as fm
from src.render import ida_sections as sections


@pytest.fixture
def method():
    return copy.deepcopy(fm.builtin_presets()['LC-MS real-time ID, MS2 only'])


def paths(rows):
    return [r.path for r in rows]


def test_a_method_does_not_differ_from_itself(method):
    assert fc.compare(method, copy.deepcopy(method)) == []


def test_migrate_keeps_the_method_until_there_is_a_migrator(method):
    migrated, notes = fc.migrate(method)
    assert migrated == method and migrated is not method and notes == []


def test_a_key_left_out_is_its_default(method):
    b = copy.deepcopy(method)
    fm.put(b, 'precursor_selection.rank_by', 'qscore')  # its default
    a = copy.deepcopy(method)
    fm.drop(a, 'precursor_selection.rank_by')
    assert fc.compare(a, b) == []


def test_only_what_changes_the_acquisition(method):
    b = copy.deepcopy(method)
    fm.put(b, 'global.method_name', 'Other')  # FLASHIda does not act on it
    fm.put(b, 'ms_settings.ms1.scan_rate', 'Turbo')  # not read with an Orbitrap
    fm.put(b, 'precursor_selection.max_precursors', 3)
    assert paths(fc.compare(method, b)) == ['precursor_selection.max_precursors']


def test_rows_follow_the_tabs_and_name_their_step(method):
    b = copy.deepcopy(method)
    fm.put(b, 'ms_settings.ms2.collision_energy', 30)
    fm.put(b, 'deconvolution.min_charge', 3)
    rows = fc.compare(method, b)
    assert [(r.tab, r.path) for r in rows] == [('Deconvolution', 'deconvolution.min_charge'),
                                               ('MS2', 'ms_settings.ms2.collision_energy')]
    assert (rows[1].a, rows[1].b) == (35, 30)


def test_switches_are_rows(method):
    b = copy.deepcopy(method)
    sections.switch(b, 'quant', True)
    fm.drop(b, 'identification')
    rows = {r.path: r for r in fc.compare(method, b)}
    assert rows['quantification.enabled'].tab == 'Quantification'
    assert (rows['identification'].a, rows['identification'].b) == ('on', 'off')
    assert fc.switch_label('identification') == 'Identify in real time'


def test_taking_settings_over(method):
    b = copy.deepcopy(method)
    fm.put(b, 'precursor_selection.max_precursors', 3)
    fm.drop(b, 'precursor_selection.rt_window')
    new = fc.adopt(b, method, paths(fc.compare(method, b)))
    assert fc.compare(method, new) == []
    assert fm.get(b, 'precursor_selection.max_precursors') == 3  # B itself is left as it was


def test_a_key_a_leaves_out_goes_from_b(method):
    a = copy.deepcopy(method)
    fm.drop(a, 'identification.threshold')
    new = fc.adopt(method, a, ['identification.threshold'])
    assert fm.get(new, 'identification.threshold', missing=None) is None


def test_identification_and_ms3_come_or_go_as_a_whole(method):
    b = copy.deepcopy(method)
    fm.drop(b, 'identification')
    assert fc.compare(method, fc.adopt(b, method, ['identification'])) == []
    assert 'identification' not in fc.adopt(method, b, ['identification'])
    with_ms3 = copy.deepcopy(method)
    sections.switch(with_ms3, 'ms3plan', True)
    assert paths(fc.compare(method, with_ms3)) == ['identification.characterization.mode']
    off = fc.adopt(with_ms3, method, ['identification.characterization.mode'])
    assert 'characterization' not in off['identification']  # no half MS3 block left


def test_unknown_keys_are_compared_as_they_are(method):
    b = copy.deepcopy(method)
    b['future_section'] = {'key': 1}
    rows = fc.compare(method, b)
    assert [(r.tab, r.path, r.a, r.b) for r in rows] == [(fc.UNKNOWN, 'future_section.key', None, 1)]
    assert 'future_section' not in fc.adopt(b, method, ['future_section.key'])


def test_all_settings_marks_the_differences(method):
    b = copy.deepcopy(method)
    fm.put(b, 'precursor_selection.max_precursors', 3)
    every = fc.all_settings(method, b)
    assert len(every) == len(fm.keys())
    marked = [p for _, p, _, _, differs in every if differs]
    assert marked == ['precursor_selection.max_precursors']
    # left out keys show the value FLASHIda works with
    row = next(r for r in every if r[1] == 'ms_settings.ms2.rf_lens')
    assert row[2] == fm.get(method, 'ms_settings.ms1.rf_lens')
