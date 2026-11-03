"""Modification library (src/render/modifications.py) and its use for the
fixed modifications of the FLASHDeconv sequence view."""

import pytest

pytest.importorskip("pyopenms")

from pyopenms import AASequence  # noqa: E402

from src.render import modifications  # noqa: E402
from src.render.sequence import getFragmentDataFromSeq  # noqa: E402


def test_library_has_unimod_entries_with_sites():
    lib = modifications.library()
    assert len(lib) > 1000
    cam = modifications.by_id('Carbamidomethyl (C)')
    assert cam['site'] == 'C'
    assert cam['unimod'] == 'UniMod:4'
    assert cam['mass'] == pytest.approx(57.021464, abs=1e-5)
    assert modifications.by_id('Acetyl (Protein N-term)')['site'] == modifications.N_TERM
    assert modifications.by_id('Amidated (C-term)')['site'] == modifications.C_TERM
    # ids are unique
    assert len({e['id'] for e in lib}) == len(lib)


def test_search_by_site_text_accession_and_mass():
    at_n = modifications.search(site='N')
    assert at_n and all(e['site'] == 'N' for e in at_n)
    # glycans are part of the library
    assert any(e['id'].startswith('HexNAc') for e in at_n)
    hits = [e['id'] for e in modifications.search('carbamidomethyl', site='C')]
    assert hits[0] == 'Carbamidomethyl (C)' and 'Phospho (C)' not in hits
    assert 'Phospho (S)' in [e['id'] for e in modifications.search('UniMod:21', site='S')]
    assert 'Oxidation (M)' in [e['id'] for e in modifications.search('15.995', site='M')]
    assert len(modifications.search(limit=3)) == 3


def test_legacy_names_map_to_library_ids():
    assert modifications.legacy_ids('Carbamidomethyl (+57)', 'L-methionine sulfoxide (+16)') == [
        'Carbamidomethyl (C)', 'Oxidation (M)',
    ]
    assert modifications.legacy_ids(None, None) == []
    assert modifications.legacy_ids('No modification', 'unknown') == []


def test_apply_fixed_residue_and_terminal_mods():
    protein = AASequence.fromString('MCKCSM')
    protein, residues, applied = modifications.apply_fixed(
        protein, ['Carbamidomethyl (C)', 'Acetyl (Protein N-term)', 'Amidated (C-term)', 'Phospho (Y)'],
    )
    assert protein.toString() == '.(Acetyl)MC(Carbamidomethyl)KC(Carbamidomethyl)SM.(Amidated)'
    assert residues == ['C']
    assert [(a['id'], a['count']) for a in applied] == [
        ('Carbamidomethyl (C)', 2), ('Acetyl (Protein N-term)', 1), ('Amidated (C-term)', 1), ('Phospho (Y)', 0),
    ]


def test_apply_fixed_first_mod_per_residue_wins_and_terminal_origin():
    protein, _, applied = modifications.apply_fixed(
        AASequence.fromString('ACK'), ['Carbamidomethyl (C)', 'Carboxymethyl (C)', 'Met->Hse (C-term M)'],
    )
    assert protein.toString() == 'AC(Carbamidomethyl)K'
    assert [a['count'] for a in applied] == [1, 0, 0]


def test_unknown_ids_are_skipped():
    protein, residues, applied = modifications.apply_fixed(AASequence.fromString('AC'), ['Nope (C)'])
    assert protein.toString() == 'AC' and residues == [] and applied == []


def test_modified_mass():
    plain = AASequence.fromString('ACDC').getMonoWeight()
    assert modifications.modified_mass('ACDC', ['Carbamidomethyl (C)']) == pytest.approx(plain + 2 * 57.021464, abs=1e-4)


def test_sequence_data_uses_library_ids_and_legacy_names():
    by_ids = getFragmentDataFromSeq('ACDKC', fixed_mods={'fixed_modifications': ['Carbamidomethyl (C)']})
    legacy = getFragmentDataFromSeq(
        'ACDKC', fixed_mods={'fixed_mod_cysteine': 'Carbamidomethyl (+57)', 'fixed_mod_methionine': None},
    )
    plain = getFragmentDataFromSeq('ACDKC', fixed_mods={'fixed_modifications': []})
    assert by_ids['fixed_modifications'] == ['C'] == legacy['fixed_modifications']
    assert by_ids['theoretical_mass'] == pytest.approx(legacy['theoretical_mass'])
    assert by_ids['theoretical_mass'] - plain['theoretical_mass'] == pytest.approx(2 * 57.021464, abs=1e-4)
    # b2 = A + C(Carbamidomethyl)
    assert by_ids['fragment_masses_b'][1][0] - plain['fragment_masses_b'][1][0] == pytest.approx(57.021464, abs=1e-4)


def test_terminal_mod_shifts_prefix_fragments():
    acetyl = getFragmentDataFromSeq('ACDK', fixed_mods={'fixed_modifications': ['Acetyl (Protein N-term)']})
    plain = getFragmentDataFromSeq('ACDK', fixed_mods={'fixed_modifications': []})
    assert acetyl['fragment_masses_b'][0][0] - plain['fragment_masses_b'][0][0] == pytest.approx(42.010565, abs=1e-4)
    assert acetyl['fragment_masses_y'][0][0] == pytest.approx(plain['fragment_masses_y'][0][0])
