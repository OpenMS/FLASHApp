"""Fragment masses of a FLASHTnT identification with a modification window
(``src/parse/tnt.py`` → ``src/render/sequence.py``).

FLASHTnT gives each modification a window of residues (``ModStart``,
``ModEnd``, 1-based in the protein). A window of one residue is a localized
modification: every ion holding that residue carries the shift, and none is
ambiguous. A longer window: the ions cut inside it may or may not carry it.
"""
import os
import sys

import pandas as pd
import pytest

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(PROJECT_ROOT)

pytest.importorskip("pyopenms")

from src.parse.tnt import _proteoform_fragments, _row_modifications  # noqa: E402

# Protein MKPEPTIDEK; the proteoform is residues 3-10 (PEPTIDEK)
PROTEIN = 'MKPEPTIDEK'
START, END = 2, 9


def _row(starts, ends, masses):
    return pd.Series({
        'ModCount': len(masses), 'ModMass': ';'.join(str(m) for m in masses),
        'ModStart': ';'.join(str(s) for s in starts), 'ModEnd': ';'.join(str(e) for e in ends),
        'ModID': None,
    })


def _kinds(row, ion):
    """Per ion of the proteoform: 'plain', 'shifted' or 'both' (ambiguous)."""
    plain = _proteoform_fragments(PROTEIN, START, END, [], [], [])
    modified = _proteoform_fragments(PROTEIN, START, END, *_row_modifications(row)[:3])
    out = []
    for base, masses in zip(plain[f'fragment_masses_{ion}'], modified[f'fragment_masses_{ion}']):
        if len(masses) == 2:
            out.append('both')
        else:
            out.append('shifted' if masses[0] == pytest.approx(base[0] + 79.966) else 'plain')
    return out


def test_a_one_residue_window_shifts_the_ions_holding_it():
    # Phospho on T6 of the protein, the 4th residue of the proteoform
    row = _row([6], [6], [79.966])
    assert _row_modifications(row)[:3] == ([6], [6], [79.966])
    # b1..b8: b1-b3 hold P,E,P; b4 on hold the T: shifted, none ambiguous
    assert _kinds(row, 'b') == ['plain'] * 3 + ['shifted'] * 5
    # y1..y8: y1-y4 hold I,D,E,K; y5 on hold the T
    assert _kinds(row, 'y') == ['plain'] * 4 + ['shifted'] * 4


def test_a_longer_window_makes_the_ions_cut_inside_it_ambiguous():
    # the shift somewhere in T6-D8 (proteoform residues 4-6)
    row = _row([6], [8], [79.966])
    # b4, b5 end inside the window; b6 on hold all of it
    assert _kinds(row, 'b') == ['plain'] * 3 + ['both'] * 2 + ['shifted'] * 3
    # y3, y4 start inside it (residues 6..8 and 5..8 of the proteoform); y5 on hold all of it
    assert _kinds(row, 'y') == ['plain'] * 2 + ['both'] * 2 + ['shifted'] * 4


def test_no_modifications():
    row = pd.Series({'ModCount': 0, 'ModMass': None, 'ModStart': None, 'ModEnd': None, 'ModID': None})
    assert _row_modifications(row) == ([], [], [], [])
    assert _kinds(row, 'b') == ['plain'] * 8
