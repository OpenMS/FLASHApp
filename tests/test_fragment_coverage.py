"""Fragment ion coverage of a proteoform (src/render/fragment_coverage.py)."""
from pathlib import Path

import polars as pl
import pytest

from src.render import fragment_coverage as fc

EXAMPLE = Path('example-data/workspaces/default/cache/files/example_lysate_identified/sequence_data.pq')


def _fragment(ion, kind, first, last, charge, offset=10):
    return {'name': f'{ion}{last - first + 1}', 'ion_type': ion, 'kind': kind, 'first': first,
            'last': last, 'start': first + offset, 'end': last + offset, 'charge': charge}


def _entry():
    return {
        'sequence': list('X' * 10 + 'PEPTIDEKR' + 'X' * 5),
        'proteoform_start': 10, 'proteoform_end': 18,
        'modifications': [{'start': 12, 'end': 12, 'mass_diff': 79.97, 'labels': 'Phospho'}],
        'annotated_fragments': [
            _fragment('b', 'N', 0, 3, 2), _fragment('b', 'N', 0, 3, 3), _fragment('b', 'N', 0, 3, 3),
            _fragment('y', 'C', 5, 8, 1), _fragment('b', 'I', 2, 5, 1),
        ],
    }


def test_residues_are_the_proteoform():
    assert ''.join(fc.residues(_entry())) == 'PEPTIDEKR'


def test_one_bar_per_fragment_with_its_charge_states():
    rows = fc.rows(7, _entry())
    by_track = {(r['track'], r['start'], r['end']): r for r in rows}
    # Terminal ions cover their prefix or suffix, internal ones their stretch
    assert set(by_track) == {('b ions', 1, 4), ('y ions', 6, 9), (fc.INTERNAL, 3, 6)}
    assert by_track[('b ions', 1, 4)]['count'] == 2
    assert all(r['reference'] == 7 for r in rows)


def test_a_fragment_names_the_modifications_it_holds():
    rows = {(r['track'], r['start']): r for r in fc.rows(0, _entry())}
    assert rows[('b ions', 1)]['sequence'] == 'b4 (Phospho)'
    assert rows[('y ions', 6)]['sequence'] == 'y4'
    assert rows[(fc.INTERNAL, 3)]['mod_positions'] == [3]


def test_a_series_covers_growing_prefixes():
    entry = {'sequence': list('PEPTIDEKR'), 'annotated_fragments': [
        _fragment('b', 'N', 0, k, 1, offset=0) for k in range(1, 8)]}
    spans = sorted((r['start'], r['end']) for r in fc.rows(0, entry))
    assert spans == [(1, k + 1) for k in range(1, 8)]


def test_mass_shift_without_one_name():
    label = fc._modification_label({'labels': 'A; B', 'mass_diff': 15.995})
    assert label == '+15.99 Da'


def test_frame_schema_and_ids():
    frame = fc.coverage_rows([(0, _entry()), (1, _entry())])
    assert frame.schema == pl.Schema(fc.SCHEMA)
    assert frame['match_id'].to_list() == list(range(6))


def test_no_fragments_gives_an_empty_frame():
    frame = fc.coverage_rows([(0, {'sequence': list('PEPTIDE')})])
    assert frame.height == 0 and frame.schema == pl.Schema(fc.SCHEMA)


def test_track_order_puts_internal_last():
    assert fc.track_order()[-1] == fc.INTERNAL


@pytest.mark.skipif(not EXAMPLE.exists(), reason='example workspace not present')
def test_example_lysate_bars_lie_on_their_proteoform():
    sequences = pl.scan_parquet(EXAMPLE)
    frame = fc.from_sequence_data(sequences)
    assert frame.height > 0
    lengths = {row['proteoform_index']: len(fc.residues(row))
               for row in sequences.collect().iter_rows(named=True)}
    for row in frame.iter_rows(named=True):
        assert 1 <= row['start'] <= row['end'] <= lengths[row['reference']]
