"""FLASHIda methods: keys from FLASHIda's schema, defaults that load, presets
kept per computer and the runs table."""
import json

from src import flashida_method as fm


def test_keys_come_from_the_schema():
    keys = fm.keys()
    assert keys['precursor_selection.rank_by'].enum == ['qscore', 'intensity', 'none']
    assert keys['global.duration'].unit == 'min'
    # scan slots resolve the shared scan shape, and ms1 takes only its own keys
    assert 'ms_settings.ms2.collision_energy' in keys
    assert 'ms_settings.ms1.collision_energy' not in keys
    assert 'deconvolution' not in keys  # objects are walked, not keys


def test_defaults_and_presets_load():
    for name, method in fm.builtin_presets().items():
        assert fm.check(method) == [], name
    assert fm.get(fm.defaults(), 'deconvolution.tol') == [10, 10, 10]


def test_check_finds_what_flashida_refuses():
    method = fm.defaults()
    method['bogus'] = 1
    fm.put(method, 'precursor_selection.rank_by', 'loudest')
    fm.put(method, 'deconvolution.tol', [10, 10])
    problems = fm.check(method)
    assert not any('bogus' in p for p in problems)  # unknown, not refused
    assert any(p.startswith('precursor_selection.rank_by') for p in problems)
    assert any(p.startswith('deconvolution.tol') for p in problems)
    assert fm.check([]) != []


def test_unknown_keys_are_kept_and_do_not_block():
    method = fm.defaults()
    method['bogus'] = 1
    fm.put(method, 'deconvolution.new_key', {'a': 1})
    fm.put(method, 'ms_settings.ms1.new_ms1_key', 3)
    assert fm.check(method) == []
    assert fm.unknown_keys(method) == ['bogus', 'deconvolution.new_key', 'ms_settings.ms1.new_ms1_key']
    merged = fm.merge(fm.defaults(), method)
    assert merged['bogus'] == 1 and fm.get(merged, 'deconvolution.new_key') == {'a': 1}
    fm.put(method, 'precursor_selection.rank_by', 'loudest')
    assert [p.split(':')[0] for p in fm.check(method)] == ['precursor_selection.rank_by']
    assert fm.unknown_keys(fm.defaults()) == []


def test_get_put_drop_and_changes():
    m = {}
    fm.put(m, 'a.b.c', 1)
    assert fm.get(m, 'a.b.c') == 1 and fm.get(m, 'a.x') is None
    fm.drop(m, 'a.b.c')
    assert m == {}
    base = {'x': {'y': 1, 'z': [1]}}
    new = {'x': {'y': 2, 'z': [1]}, 'w': True}
    assert fm.changes(base, new) == [('w', None, True), ('x.y', 1, 2)]


def test_user_presets_are_kept_per_computer(tmp_path):
    assert fm.user_presets(tmp_path) == {}
    fm.save_user_preset(tmp_path, 'Histone MS3', {'global': {'duration': 60}})
    assert fm.user_presets(tmp_path) == {'Histone MS3': {'global': {'duration': 60}}}
    presets = fm.all_presets(tmp_path)
    assert list(presets)[0] == fm.DEFAULTS_PRESET and 'Histone MS3' in presets
    # a user preset never hides a built-in one
    fm.save_user_preset(tmp_path, fm.DEFAULTS_PRESET, {})
    assert fm.all_presets(tmp_path)[fm.DEFAULTS_PRESET] == fm.defaults()
    fm.delete_user_preset(tmp_path, 'Histone MS3')
    assert 'Histone MS3' not in fm.user_presets(tmp_path)


def test_runs_table(tmp_path):
    assert fm.load_runs(tmp_path) == []
    runs = [{'name': 'a', 'preset': 'p', 'method': {}}, {'name': 'b', 'method': 'broken'}]
    fm.save_runs(tmp_path, runs)
    assert fm.load_runs(tmp_path) == runs[:1]
    (tmp_path / fm.RUNS_FILE).write_text('not json')
    assert fm.load_runs(tmp_path) == []


def test_run_names():
    assert fm.name_problem('', []) != ''
    assert fm.name_problem('a/b', []) != ''
    assert fm.name_problem('a.raw', []) != ''
    assert fm.name_problem('Sample', ['sample']) != ''
    assert fm.name_problem('Sample_2', ['sample']) == ''
    assert fm.method_file_name('a:b') == 'a_b.json'
