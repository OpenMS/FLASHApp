"""Fragments around unlocalized modifications (src/render/fragment_ambiguity.py)."""

import polars as pl

from src.render.fragment_ambiguity import settle, settle_fragment_masses

SHIFT = -2.0


def _entry():
    # Proteoform of 8 residues starting at database residue 3 (0-based 2); a
    # shift somewhere on residues 3..5 of the proteoform (database 4..6, 0-based
    # 4..6 -> start 4, end 6). Masses as the parser stores them: both values
    # for a cleavage in the range (including its ends), the shifted one alone
    # after it.
    prefix = []
    for i in range(1, 9):
        base = 100.0 * i
        if 3 <= i <= 5:
            prefix.append([base, base + SHIFT])
        elif i > 5:
            prefix.append([base + SHIFT])
        else:
            prefix.append([base])
    suffix = []
    for k in range(8):
        first = 8 - k
        base = 1000.0 + k
        if 3 <= first <= 5:
            suffix.append([base, base + SHIFT])
        elif first < 3:
            suffix.append([base + SHIFT])
        else:
            suffix.append([base])
    entry = {'start_position': 2, 'modifications': [{'start': 4, 'end': 6, 'mass_diff': SHIFT, 'labels': ''}]}
    for ion in 'abc':
        entry[f'fragment_masses_{ion}'] = [list(m) for m in prefix]
    for ion in 'xyz':
        entry[f'fragment_masses_{ion}'] = [list(m) for m in suffix]
    return entry


def test_prefix_fragments_inside_the_range_have_no_mass():
    settled = settle_fragment_masses(_entry())['fragment_masses_b']
    assert settled[:2] == [[100.0], [200.0]]
    # b3, b4 end inside residues 3..5: the shift may or may not be in them
    assert settled[2] == [] and settled[3] == []
    # b5 covers the whole range, so it carries the shift
    assert settled[4] == [500.0 + SHIFT]
    assert settled[7] == [800.0 + SHIFT]


def test_suffix_fragments_inside_the_range_have_no_mass():
    settled = settle_fragment_masses(_entry())['fragment_masses_y']
    # y1..y3 start at residues 8..6, after the range
    assert settled[:3] == [[1000.0], [1001.0], [1002.0]]
    # y4, y5 start at residues 5 and 4, inside it
    assert settled[3] == [] and settled[4] == []
    # y6 starts at residue 3 and covers the whole range
    assert settled[5] == [1005.0 + SHIFT]


def test_localized_modifications_and_no_modifications_are_unchanged():
    entry = _entry()
    entry['modifications'] = [{'start': 4, 'end': 4, 'mass_diff': SHIFT, 'labels': ''}]
    assert settle_fragment_masses(entry)['fragment_masses_c'] == entry['fragment_masses_c']
    entry['modifications'] = []
    assert settle_fragment_masses(entry)['fragment_masses_z'] == entry['fragment_masses_z']


def test_settled_masses_are_settled_again_unchanged():
    once = {**_entry(), **settle_fragment_masses(_entry())}
    assert settle_fragment_masses(once) == settle_fragment_masses(_entry())


def test_settle_on_a_frame():
    entry = _entry()
    lf = pl.DataFrame([{'proteoform_index': 0, **entry}]).lazy()
    row = settle(lf).collect().row(0, named=True)
    assert row['fragment_masses_b'][2] == []
    assert row['fragment_masses_y'][5] == [1005.0 + SHIFT]
    assert row['proteoform_index'] == 0
