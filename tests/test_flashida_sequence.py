"""FLASHIda sequences: export and import of the runs table, raw file names
from an Xcalibur sequence, and run names made from a pattern."""
import json

import pytest

from src import flashida_folder as ff
from src import flashida_method as fm
from src import flashida_sequence as fs


def _runs(*names):
    return [{'name': n, 'method_file': 'M'} for n in names]


def test_export_and_import_round_trip(tmp_path):
    ff.write(tmp_path, 'M1', {'global': {'duration': 60}})
    runs = [{'name': 'A_R1', 'method_file': 'M1'}, {'name': 'A_R2', 'method_file': 'M1'}]
    text = fs.export_runs(runs, tmp_path)
    data = json.loads(text)
    assert data['format'] == fs.FORMAT and data['version'] == 2
    # the method travels with the file, for a PC whose folder lacks it
    assert data['runs'][0]['method'] == {'global': {'duration': 60}}
    parsed = fs.parse_runs(text.encode('utf-8-sig'))  # as uploaded, with a BOM
    assert [r['method_file'] for r in parsed] == ['M1', 'M1']
    # imported on another PC: the missing method file is written
    other = tmp_path / 'other'
    other.mkdir()
    new, written, missing = fs.into_folder(parsed, other)
    assert new == runs and written == ['M1'] and missing == []
    assert ff.read(other, 'M1') == {'global': {'duration': 60}}
    # a file without the method, naming one the folder lacks: listed as missing
    _, _, missing = fs.into_folder([{'name': 'x', 'method_file': 'Gone'}], other)
    assert missing == ['Gone']


def test_import_a_file_of_the_first_version(tmp_path):
    """Files exported before the method folder carry the method and its preset."""
    v1 = json.dumps({'format': fs.FORMAT, 'version': 1, 'runs': [
        {'name': 'A', 'preset': 'Wash', 'method': {'global': {'duration': 5}}}]})
    new, written, _ = fs.into_folder(fs.parse_runs(v1), tmp_path)
    assert new == [{'name': 'A', 'method_file': 'Wash'}] and written == ['Wash']


def test_csv_export():
    text = fs.export_csv([{'name': 'A_R1', 'method_file': 'M1'}, {'name': 'wash, 1', 'method_file': 'Wash'}])
    assert text.splitlines() == ['Row,Raw file name,Method', '1,A_R1,M1', '2,"wash, 1",Wash']


def test_xcalibur_csv_export_and_import():
    runs = [{'name': 'A_R1', 'method_file': 'M1'}, {'name': 'B_R1', 'method_file': 'M2'}]
    text = fs.export_xcalibur_csv(runs, r'C:\Data', r'C:\Xcalibur\methods\TD')
    lines = text.splitlines()
    assert lines[0] == 'Bracket Type=4,' and lines[1].split(',') == fs.XCALIBUR_COLUMNS
    assert lines[2].startswith(r'Unknown,A_R1,,C:\Data,C:\Xcalibur\methods\TD,')
    # read back: the names, without methods (Xcalibur knows none of ours)
    assert fs.csv_runs(text.encode('utf-8')) == ([{'name': 'A_R1', 'method_file': None},
                                                  {'name': 'B_R1', 'method_file': None}], False)


def test_the_table_csv_is_read_back_with_its_methods():
    runs = [{'name': 'A_R1', 'method_file': 'M1'}, {'name': 'wash, 1', 'method_file': 'Wash'}]
    assert fs.csv_runs(fs.export_csv(runs)) == (runs, True)
    with pytest.raises(ValueError):
        fs.csv_runs('a,b\n1,2\n')


@pytest.mark.parametrize('text', ['not json', '{"format": "other", "runs": []}', '{"runs": 3}',
                                  '[{"name": "a"}]', '[1]', '[{"name": "a", "method_file": ""}]'])
def test_import_refuses_what_is_not_a_sequence(text):
    with pytest.raises(ValueError):
        fs.parse_runs(text)


def test_append_or_replace_and_duplicates():
    old, new = _runs('A_R1', 'Wash_1'), _runs('a_r1', 'B_R1', 'B_R1')
    runs, dup = fs.combine(old, new, 'append')
    assert [r['name'] for r in runs] == ['A_R1', 'Wash_1', 'a_r1', 'B_R1', 'B_R1']
    assert dup == ['A_R1', 'B_R1']  # compared as Windows does, ignoring case
    runs, dup = fs.combine(old, new[:2], 'replace')
    assert [r['name'] for r in runs] == ['a_r1', 'B_R1'] and dup == []
    with pytest.raises(ValueError):
        fs.combine(old, new, 'merge')
    assert fs.duplicate_names(_runs('', '', 'x')) == []  # unnamed runs are not duplicates


def _sld(names):
    blobs = ''.join('\x00Finnigan{"ReInject":{"OriginalRawFileNameWithoutExtension":%s,'
                    '"DisplayRowNumber":%d}}' % (json.dumps(n), i) for i, n in enumerate(names, 1))
    return b'\x01\x02' + blobs.encode('utf-16-le')


def test_names_from_any_xcalibur_sequence():
    names = fs.sld_names(_sld(['PK_BSA_R1', 'PK_wash_1']))
    assert names == ['PK_BSA_R1', 'PK_wash_1']
    assert fs.runs_from_names(names, 'Wash') == [{'name': 'PK_BSA_R1', 'method_file': 'Wash'},
                                                 {'name': 'PK_wash_1', 'method_file': 'Wash'}]


def test_samples_with_their_own_values():
    samples = fs.parse_samples('BSA\n\n Ecoli ; comment=top3; replicates=3\nHeLa;comment=', 'c', 2)
    assert samples == [{'sample': 'BSA', 'comment': 'c', 'replicates': 2},
                       {'sample': 'Ecoli', 'comment': 'top3', 'replicates': 3},
                       {'sample': 'HeLa', 'comment': '', 'replicates': 2}]
    for bad in ['BSA; replicates=two', 'BSA; replicates=0', 'BSA; colour=red', '; comment=x']:
        with pytest.raises(ValueError, match='Line 1'):
            fs.parse_samples(bad)


def test_pattern_names_and_order():
    samples = fs.parse_samples('BSA\nEcoli; comment=top3; replicates=3', '', 2)
    names = [n['name'] for n in fs.pattern_names(samples, 'PK', '20261110', '', 'R{n}', 'sample')]
    assert names == ['PK_20261110_BSA_R1', 'PK_20261110_BSA_R2', 'PK_20261110_Ecoli_top3_R1',
                     'PK_20261110_Ecoli_top3_R2', 'PK_20261110_Ecoli_top3_R3']
    names = [n['name'] for n in fs.pattern_names(samples, '', '', 'x', 'rep{nn}', 'replicate')]
    assert names == ['BSA_rep01_x', 'Ecoli_top3_rep01_x', 'BSA_rep02_x', 'Ecoli_top3_rep02_x',
                     'Ecoli_top3_rep03_x']


def test_a_wash_every_few_runs_but_not_after_the_last():
    samples = fs.parse_samples('A\nB', '', 2)
    runs = fs.pattern_names(samples, 'PK', '20261110', wash_every=2)
    assert [(n['name'], n['wash']) for n in runs] == [
        ('PK_20261110_A_R1', False), ('PK_20261110_A_R2', False), ('PK_20261110_wash_1', True),
        ('PK_20261110_B_R1', False), ('PK_20261110_B_R2', False)]
    # an empty placeholder leaves no stray separator
    assert fs.wash_run_name('{prefix}_{date}_wash_{n}', '', '20261110', '', 3) == '20261110_wash_3'
    assert fs.unknown_placeholders('{prefix}_{sample}_{n}') == ['sample']


def test_wash_preset_is_built_in():
    presets = fm.builtin_presets()
    assert fs.WASH_PRESET in presets
    assert fm.get(presets[fs.WASH_PRESET], 'global.method_name') == 'Wash'


def test_runs_of_several_kinds_between_the_samples():
    samples = fs.parse_samples('A\nB\nC', '', 2)
    inserts = [{'kind': 'wash', 'every': 2, 'name': fs.DEFAULT_WASH_NAME},
               {'kind': 'qc', 'every': 3, 'name': fs.INSERT_KINDS['qc'][1]}]
    runs = fs.pattern_names(samples, 'PK', '20261113', inserts=inserts)
    assert [(n['name'], n['kind']) for n in runs] == [
        ('PK_20261113_A_R1', 'sample'), ('PK_20261113_A_R2', 'sample'), ('PK_20261113_wash_1', 'wash'),
        ('PK_20261113_B_R1', 'sample'), ('PK_20261113_QC_1', 'qc'), ('PK_20261113_B_R2', 'sample'),
        ('PK_20261113_wash_2', 'wash'), ('PK_20261113_C_R1', 'sample'), ('PK_20261113_C_R2', 'sample')]
    # two kinds on the same place come in the order given
    both = fs.with_inserts(['a', 'b', 'c'], [{'kind': 'qc', 'every': 2}, {'kind': 'wash', 'every': 2}])
    assert [n['kind'] for n in both] == ['sample', 'sample', 'qc', 'wash', 'sample']


def test_inserted_runs_take_the_first_free_number():
    samples = fs.parse_samples('A\nB\nC', '', 1)
    inserts = [{'kind': 'wash', 'every': 1, 'name': fs.DEFAULT_WASH_NAME}]
    runs = fs.pattern_names(samples, '', '20261113', inserts=inserts,
                            taken=['20261113_wash_1', '20261113_WASH_2'])
    assert [n['name'] for n in runs if n['wash']] == ['20261113_wash_3', '20261113_wash_4']
    # a name without {n} cannot be numbered: it stays, and the table marks the duplicate
    fixed = fs.with_inserts(['a', 'b', 'c'], [{'kind': 'blank', 'every': 1, 'name': 'blank'}])
    assert [n['name'] for n in fixed if n['kind'] == 'blank'] == ['blank', 'blank']


def test_random_order_shuffles_only_the_samples():
    samples = fs.parse_samples('A\nB\nC\nD', '', 3)
    inserts = [{'kind': 'wash', 'every': 3, 'name': '{date}_wash_{n}'}]
    one = fs.pattern_names(samples, '', 'd', order='random', seed=1, inserts=inserts)
    again = fs.pattern_names(samples, '', 'd', order='random', seed=1, inserts=inserts)
    other = fs.pattern_names(samples, '', 'd', order='random', seed=2, inserts=inserts)
    in_order = fs.pattern_names(samples, '', 'd', inserts=inserts)
    assert one == again and one != other
    for runs in (one, other):
        # the washes keep their places, the same sample runs are all there
        assert [i for i, n in enumerate(runs) if n['wash']] == [3, 7, 11]
        assert sorted(n['name'] for n in runs) == sorted(n['name'] for n in in_order)
