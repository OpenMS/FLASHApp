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
