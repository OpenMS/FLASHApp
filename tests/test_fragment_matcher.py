"""The mass lists the Match fragments page reads (src/fragment_matcher.py)."""

import json

import pytest

from src import fragment_matcher as fm


def test_plain_list_with_header_and_intensities():
    parsed = fm.parse_mass_list('mass\tintensity\n1043.5 2.5e4\n2511.37,18000\n\n# note\n980.1;10\n')
    assert parsed['masses'] == [1043.5, 2511.37, 980.1]
    assert parsed['intensities'] == [25000.0, 18000.0, 10.0]
    assert parsed['skipped'] == 1
    assert parsed['precursor_mass'] is None


def test_semicolon_separated_on_one_line():
    parsed = fm.parse_mass_list('1043.5; 2511.37;980.1;\n')
    assert parsed['masses'] == [1043.5, 2511.37, 980.1]
    assert parsed['intensities'] is None
    parsed = fm.parse_mass_list('1043.5 25000; 2511.37 18000')
    assert parsed['masses'] == [1043.5, 2511.37]
    assert parsed['intensities'] == [25000.0, 18000.0]


def test_masses_only_and_mixed_columns():
    parsed = fm.parse_mass_list('100.5\n200.25 7\n')
    assert parsed['masses'] == [100.5, 200.25]
    assert parsed['intensities'] is None


def test_msalign_first_spectrum():
    text = (
        'BEGIN IONS\nID=0\nSCANS=12\nPRECURSOR_MASS=9529.18\n'
        '1043.51\t25000\t2\n2511.37\t18000\t3\nEND IONS\n'
        'BEGIN IONS\nPRECURSOR_MASS=1.0\n5.0\t1\t1\nEND IONS\n'
    )
    parsed = fm.parse_mass_list(text.encode())
    assert parsed['masses'] == [1043.51, 2511.37]
    assert parsed['intensities'] == [25000.0, 18000.0]
    assert parsed['precursor_mass'] == 9529.18


def test_no_masses():
    with pytest.raises(fm.MassListError):
        fm.parse_mass_list('mass\nabc\n')
    with pytest.raises(fm.MassListError):
        fm.parse_mass_list('BEGIN IONS\nPRECURSOR_MASS=5\nEND IONS\n')


def test_protonated_masses_become_neutral():
    assert fm.to_neutral([1000.0], fm.PROTONATED) == [pytest.approx(1000.0 - 1.007276466812)]
    assert fm.to_neutral([1000.0], fm.NEUTRAL) == [1000.0]


def test_store_and_load(tmp_path):
    data = fm.store_input(tmp_path, [3.0, 1.0, 2.0], [30, 10, 20], 9.5,
                          mass_type=fm.AVERAGE, sequence='PEPTIDE', name='x')
    assert data['masses'] == [1.0, 2.0, 3.0]
    assert data['intensities'] == [10.0, 20.0, 30.0]
    assert fm.load_input(tmp_path) == json.loads(json.dumps(data))
    row = fm.spectrum_row(data)
    assert row == {'index': 0, 'MonoMass': [1.0, 2.0, 3.0], 'PrecursorMass': 9.5}
    assert fm.input_stamp(data) != fm.input_stamp({**data, 'mass_type': fm.MONOISOTOPIC})
    with pytest.raises(ValueError):
        fm.matcher_input([1.0], mass_type='other')


def test_clear_input(tmp_path):
    assert fm.clear_input(tmp_path) is False
    fm.store_input(tmp_path, [1000.0, 2000.0])
    assert fm.load_input(tmp_path) is not None
    assert fm.clear_input(tmp_path) is True
    assert fm.load_input(tmp_path) is None


def test_load_missing_or_broken(tmp_path):
    assert fm.load_input(tmp_path) is None
    fm.input_path(tmp_path).parent.mkdir(parents=True)
    fm.input_path(tmp_path).write_text('{')
    assert fm.load_input(tmp_path) is None


def test_average_fragment_masses():
    from src.render.sequence import getFragmentDataFromSeq
    mono = getFragmentDataFromSeq('PEPTIDE', fixed_mods={'fixed_modifications': []})
    avg = getFragmentDataFromSeq('PEPTIDE', fixed_mods={'fixed_modifications': []}, average=True)
    assert mono['theoretical_mass'] == pytest.approx(799.36, abs=0.01)
    assert avg['theoretical_mass'] == pytest.approx(799.82, abs=0.01)
    assert avg['fragment_masses_b'][1][0] > mono['fragment_masses_b'][1][0]


def _matches(proforma, mass, names, cleavages=50.0, tolerance=10):
    return {
        'spectrum': 0,
        'tolerance': {'value': tolerance, 'unit': 'ppm'},
        'fragments': [
            {'name': n, 'ion': n[0], 'kind': 'N' if n[0] in 'abc' else 'C',
             'mass': 100.0 + i, 'massIndex': i}
            for i, n in enumerate(names)
        ] + [{'name': 'b2-y3', 'ion': 'by', 'kind': 'I', 'mass': 5.0, 'massIndex': None}],
        'proteoform': {'proforma': proforma, 'sequence': 'PEPTIDE', 'length': 7,
                       'theoreticalMass': mass, 'cleavagePercentage': cleavages},
    }


def test_hypothesis_row_from_the_view_matches():
    row = fm.hypothesis(_matches('PEPT[+79.97]IDE', 1000.0, ['b2', 'b3', 'y2', 'b3']), 1000.01)
    # Internal fragments are not counted; an ion matched on two masses is
    # two rows, as in the view's table
    assert row['matched_fragments'] == 4
    assert (row['n_terminal'], row['c_terminal']) == (3, 1)
    assert row['proforma'] == 'PEPT[+79.97]IDE'
    assert row['delta_da'] == pytest.approx(0.01)
    assert row['delta_ppm'] == pytest.approx(10.0)
    assert row['tolerance'] == '10 ppm'
    # No precursor, no Δ; nothing matched yet or an older view, no row
    assert fm.hypothesis(_matches('PEPTIDE', 1000.0, ['b2']))['delta_da'] is None
    assert fm.hypothesis(None) is None
    assert fm.hypothesis({'fragments': []}) is None


def test_hypotheses_kept_per_mass_list(tmp_path):
    data = fm.store_input(tmp_path, [100.0, 200.0], precursor_mass=1000.0)
    first = fm.hypothesis(_matches('PEPTIDE', 1000.5, ['b2']), 1000.0)
    second = fm.hypothesis(_matches('PEPT[+79.97]IDE', 1000.0, ['b2', 'y3', 'y4'], 70), 1000.0)
    fm.add_hypothesis(tmp_path, data, first)
    fm.add_hypothesis(tmp_path, data, second)
    rows = fm.load_hypotheses(tmp_path, data)
    assert [r['proforma'] for r in rows] == ['PEPTIDE', 'PEPT[+79.97]IDE']
    assert fm.best_hypothesis(rows) == 1

    # The same proteoform tried again replaces its row in place
    again = fm.hypothesis(_matches('PEPTIDE', 1000.5, ['b2', 'b3', 'b4', 'y5'], 80), 1000.0)
    rows = fm.add_hypothesis(tmp_path, data, again)
    assert [r['matched_fragments'] for r in rows] == [4, 3]
    assert fm.best_hypothesis(rows) == 0

    assert [r['proforma'] for r in fm.remove_hypotheses(tmp_path, data, [0])] == ['PEPT[+79.97]IDE']

    # Another mass list starts an empty table
    other = fm.store_input(tmp_path, [300.0])
    assert fm.load_hypotheses(tmp_path, other) == []
    assert fm.best_hypothesis([]) is None
