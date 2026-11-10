"""FLASHIda method folder: methods as files, presets saved once, favorites,
runs that name their method file, moving and duplicating runs, durations."""
import json

import pytest

from src import flashida_folder as ff
from src import flashida_method as fm


def test_methods_are_the_json_files_of_the_folder(tmp_path):
    assert ff.list_methods(tmp_path / 'missing') == []
    ff.write(tmp_path, 'b_method', {'global': {}})
    ff.write(tmp_path, 'A_method', {'global': {}})
    (tmp_path / ff.FAVORITES_FILE).write_text('[]')
    (tmp_path / 'notes.txt').write_text('x')
    assert ff.list_methods(tmp_path) == ['A_method', 'b_method']  # sorted ignoring case, no dot files
    with pytest.raises(FileExistsError):
        ff.write(tmp_path, 'A_method', {})
    ff.write(tmp_path, 'A_method', {'x': 1}, overwrite=True)
    assert ff.read(tmp_path, 'A_method') == {'x': 1}
    (tmp_path / 'bad.json').write_text('[1]')
    with pytest.raises(ValueError, match='bad.json'):
        ff.read(tmp_path, 'bad')


def test_names(tmp_path):
    assert ff.stem('Histone MS2.json') == 'Histone MS2'
    assert ff.stem('a/b') == 'a_b'
    assert ff.name_problem('') and ff.name_problem('a:b') and ff.name_problem('.hidden')
    assert ff.name_problem('Histone MS2') == ''


def test_presets_are_saved_once_and_never_overwrite(tmp_path):
    ff.write(tmp_path, 'Wash', {'mine': True})
    written = ff.seed_presets(tmp_path, fm.builtin_presets())
    assert 'Wash' not in written and fm.DEFAULTS_PRESET in written
    assert ff.read(tmp_path, 'Wash') == {'mine': True}
    assert ff.seed_presets(tmp_path, fm.builtin_presets()) == []


def test_favorites_come_first(tmp_path):
    assert ff.toggle_favorite(tmp_path, 'B') is True
    assert ff.toggle_favorite(tmp_path, 'C') is True
    assert ff.favorites(tmp_path) == ['B', 'C']
    assert ff.ordered(['A', 'B', 'C', 'D'], ff.favorites(tmp_path)) == ['B', 'C', 'A', 'D']
    assert ff.toggle_favorite(tmp_path, 'B') is False
    assert json.loads((tmp_path / ff.FAVORITES_FILE).read_text()) == ['C']


def test_runs_get_their_method_file_as_it_is_now(tmp_path):
    ff.write(tmp_path, 'M', {'global': {'duration': 30}})
    run = {'name': 'r', 'method_file': 'M'}
    assert ff.run_method(run, tmp_path) == {'global': {'duration': 30}}
    ff.write(tmp_path, 'M', {'global': {'duration': 45}}, overwrite=True)
    assert ff.duration(ff.run_method(run, tmp_path)) == 45
    assert ff.run_method({'name': 'r', 'method_file': 'Gone'}, tmp_path) is None
    # a table saved before the method folder keeps the method in the run
    assert ff.run_method({'name': 'r', 'method': {'a': 1}}, tmp_path) == {'a': 1}
    problems = ff.method_problems([run, {'name': 'x', 'method_file': 'Gone'}, {'name': 'y'}], tmp_path)
    assert [i for i, _ in problems] == [1, 2]


def test_methods_kept_in_the_table_move_to_the_folder(tmp_path):
    ff.write(tmp_path, 'Same', {'a': 1})
    runs = [{'name': 'r1', 'preset': 'p', 'method': {'a': 1}},   # equal to a file: uses it
            {'name': 'r2', 'preset': 'Same', 'method': {'b': 2}},  # name taken: numbered
            {'name': 'r3', 'method': {'c': 3}},                     # no preset: the run's name
            {'name': 'r4', 'method_file': 'Same'}]
    out, written = ff.adopt_embedded(runs, tmp_path)
    assert out == [{'name': 'r1', 'method_file': 'Same'}, {'name': 'r2', 'method_file': 'Same (2)'},
                   {'name': 'r3', 'method_file': 'r3'}, {'name': 'r4', 'method_file': 'Same'}]
    assert written == ['Same (2)', 'r3'] and ff.read(tmp_path, 'Same (2)') == {'b': 2}


def test_duplicate_counts_the_replicate_up():
    assert ff.next_name('PK_BSA_R1', ['PK_BSA_R1', 'PK_BSA_R2']) == 'PK_BSA_R3'
    assert ff.next_name('A_rep01', ['A_rep01']) == 'A_rep02'
    assert ff.next_name('Blank', ['Blank', 'blank_2']) == 'Blank_3'
    runs = [{'name': 'A_R1', 'method_file': 'M'}, {'name': 'B_R1', 'method_file': 'N'}]
    out = ff.duplicate(runs, 0)
    assert [r['name'] for r in out] == ['A_R1', 'A_R2', 'B_R1'] and out[1]['method_file'] == 'M'
    assert runs[0]['name'] == 'A_R1' and len(runs) == 2  # the table given is left as it was


def test_move():
    runs = ['a', 'b', 'c']
    assert ff.move(runs, 0, 1) == ['b', 'a', 'c']
    assert ff.move(runs, 2, -1) == ['a', 'c', 'b']
    assert ff.move(runs, 0, -1) == runs and ff.move(runs, 2, 1) == runs


def test_total_duration_and_end(tmp_path):
    ff.write(tmp_path, 'Long', {'global': {'duration': 120}})
    ff.write(tmp_path, 'Wash', {'global': {}})
    runs = [{'name': 'a', 'method_file': 'Long'}, {'name': 'b', 'method_file': 'Long'},
            {'name': 'w', 'method_file': 'Wash'}]
    assert ff.total_minutes(runs, tmp_path) == (240.0, 1)
    assert ff.format_minutes(240) == '4 h 00 min' and ff.format_minutes(45) == '45 min'
    assert ff.format_end(90, now=0)  # a weekday and a time
