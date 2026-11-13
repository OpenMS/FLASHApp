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


def _isotope_peaks(sequence, charges=(3, 4, 5)):
    """m/z and intensity of a peptide's isotope pattern at several charges."""
    import pyopenms as oms
    peptide = oms.AASequence.fromString(sequence)
    mono = peptide.getMonoWeight()
    pattern = peptide.getFormula().getIsotopeDistribution(oms.CoarseIsotopePatternGenerator(8))
    mz, intensities = [], []
    for z in charges:
        for k, peak in enumerate(pattern.getContainer()):
            mz.append((mono + k * 1.00235 + z * fm.PROTON) / z)
            intensities.append(peak.getIntensity() * 1e6)
    return mono, mz, intensities


def test_deconvolve_mz_finds_the_monoisotopic_mass():
    mono, mz, intensities = _isotope_peaks('MNKTQLIDVIAEKAELSKTQAKAAL')
    found = fm.deconvolve_mz(mz, intensities)
    assert len(found['masses']) == 1
    assert found['masses'][0] == pytest.approx(mono, abs=0.01)
    assert found['intensities'][0] > 0


def test_deconvolve_mz_without_peaks():
    assert fm.deconvolve_mz([]) == {'masses': [], 'intensities': []}


class _Response:
    def __init__(self, text):
        self.text = text


def test_uniprot_entry(monkeypatch):
    from src.workflow import uniprot
    asked = []

    def get(path, params):
        asked.append(path)
        return _Response('>sp|P0ACF0|DBHA_ECOLI DNA-binding protein HU-alpha '
                         'OS=Escherichia coli (strain K12) OX=83333 GN=hupA PE=1 SV=2\nMNKT\nQLID\n')

    monkeypatch.setattr(uniprot, '_get', get)
    entry = fm.uniprot_entry('sp|P0ACF0|DBHA_ECOLI')
    assert entry == {'accession': 'P0ACF0', 'entry': 'DBHA_ECOLI',
                     'name': 'DNA-binding protein HU-alpha',
                     'organism': 'Escherichia coli (strain K12)', 'sequence': 'MNKTQLID'}
    fm.uniprot_entry('p0acf0-2')
    assert asked == ['/uniprotkb/P0ACF0.fasta', '/uniprotkb/P0ACF0-2.fasta']


def test_uniprot_entry_errors(monkeypatch):
    from src.workflow import uniprot

    def missing(path, params):
        raise uniprot.UniProtError('UniProt answered 404')

    monkeypatch.setattr(uniprot, '_get', missing)
    with pytest.raises(fm.UniProtLookupError, match='not a UniProt accession'):
        fm.uniprot_entry('hello')
    with pytest.raises(fm.UniProtLookupError, match='has no entry Q9ZZZ9'):
        fm.uniprot_entry('Q9ZZZ9')


def test_range_fragments_with_and_without_the_shift():
    """A range of residues 3-5 (1-based): an ion holding part of it gets both
    masses (ambiguous), one holding all of it only the shifted mass."""
    from src.render.sequence import getFragmentDataFromSeq
    plain = getFragmentDataFromSeq('PEPTIDEK', fixed_mods={'fixed_modifications': []})
    ranged = getFragmentDataFromSeq('PEPTIDEK', fixed_mods={'fixed_modifications': []},
                                    modifications=[(3, 5, 16.0)])

    def kinds(ion):
        out = []
        for base, masses in zip(plain[f'fragment_masses_{ion}'], ranged[f'fragment_masses_{ion}']):
            if len(masses) == 2:
                out.append('both')
            else:
                out.append('shifted' if masses[0] == pytest.approx(base[0] + 16.0) else 'plain')
        return out

    # b1..b8: b1, b2 before the range; b3, b4 cut inside it; b5.. hold all of it
    assert kinds('b') == ['plain', 'plain', 'both', 'both', 'shifted', 'shifted', 'shifted', 'shifted']
    # y1..y8: y1-y3 after the range; y4, y5 cut inside it; y6.. hold all of it
    assert kinds('y') == ['plain', 'plain', 'plain', 'both', 'both', 'shifted', 'shifted', 'shifted']
