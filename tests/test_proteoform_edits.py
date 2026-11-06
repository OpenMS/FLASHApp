"""Edits of identified proteoforms made in the sequence view."""

from src.render import proteoform_edits as pe


SEARCH = {
    'sequence': list('MKSTAK'),
    'start_position': 1,
    'end_position': 5,
    'observed_mass': 700.0,
    'coverage': [0, 1, 1, 0, 0, 0],
    'maxCoverage': 2,
    'modifications': [
        {'start': 2, 'end': 4, 'mass_diff': 79.9663, 'labels': 'Phospho'},
        {'start': 1, 'end': 1, 'mass_diff': 42.0106, 'labels': ''},
    ],
}


def test_edit_starts_from_the_search_result():
    edits = pe.apply({}, {'action': 'edit', 'nonce': 1}, 7, SEARCH)
    assert edits[pe.NONCE] == 1
    edit = edits[7]
    assert edit['sequence'] == 'MKSTAK'
    assert (edit['start_position'], edit['end_position']) == (1, 5)
    # Each modification on the first residue of its range, in sequence order
    assert [(m['start'], m['end'], m['mass_diff']) for m in edit['modifications']] == [
        (1, 1, 42.0106), (2, 2, 79.9663),
    ]
    assert edit['modifications'][1]['labels'] == 'Phospho'


def test_an_edit_is_applied_once():
    edits = pe.apply({}, {'action': 'edit', 'nonce': 1}, 7, SEARCH)
    assert pe.apply(edits, {'action': 'edit', 'nonce': 1}, 7, SEARCH) is None
    assert pe.apply(edits, None, 7, SEARCH) is None
    assert pe.apply(edits, {'action': 'edit'}, 7, SEARCH) is None


def test_update_and_reset():
    edits = pe.apply({}, {'action': 'edit', 'nonce': 1}, 7, SEARCH)
    update = {
        'action': 'update', 'nonce': 2, 'sequence': 'mkstak',
        'start_position': 1, 'end_position': 5,
        'modifications': [
            {'start': 3, 'end': 3, 'mass_diff': 79.9663, 'labels': 'Phospho (T)'},
            {'start': 4, 'end': 4, 'mass_diff': 0, 'labels': ''},
            {'start': 99, 'end': 99, 'mass_diff': 1.0, 'labels': ''},
            {'start': -1, 'end': -1, 'mass_diff': 42.0106, 'labels': 'Acetyl'},
        ],
    }
    edits = pe.apply(edits, update, 7, SEARCH)
    assert edits[7]['sequence'] == 'MKSTAK'
    # Zero shifts and positions off the sequence are dropped
    assert [m['start'] for m in edits[7]['modifications']] == [-1, 3]
    # Not residues: refused
    assert pe.apply(edits, {**update, 'nonce': 3, 'sequence': 'MK1'}, 7, SEARCH) is None
    # Positions beyond a shorter sequence fall back to the whole sequence
    short = pe.apply(edits, {**update, 'nonce': 4, 'sequence': 'MKS'}, 7, SEARCH)
    assert (short[7]['start_position'], short[7]['end_position']) == (0, 2)
    reset = pe.apply(edits, {'action': 'reset', 'nonce': 5}, 7, SEARCH)
    assert 7 not in reset and reset[pe.NONCE] == 5


def test_entry_of_an_edited_proteoform():
    edits = pe.apply({}, {'action': 'edit', 'nonce': 1}, 7, SEARCH)
    seen = []

    def fragments(sequence):
        seen.append(sequence)
        return {'theoretical_mass': 600.0, 'fragment_masses_b': [[1.0]] * len(sequence)}

    entry = pe.entry(edits[7], SEARCH, fragments)
    # Fragments of the proteoform's residues only
    assert seen == ['KSTAK']
    assert entry['edited'] is True
    assert entry['annotated_fragments'] is None
    assert entry['observed_mass'] == 700.0
    assert entry['sequence'] == list('MKSTAK')
    # Same sequence: the tag coverage stays
    assert entry['coverage'] == SEARCH['coverage'] and entry['maxCoverage'] == 2

    changed = dict(edits[7], sequence='MKSAAK')
    entry = pe.entry(changed, SEARCH, fragments)
    assert entry['coverage'] == [0.0] * 6 and entry['maxCoverage'] == 0


class _Files:
    def __init__(self, stored=None):
        self.stored = stored

    def result_exists(self, dataset_id, name):
        return self.stored is not None

    def get_results(self, dataset_id, names):
        return {names[0]: self.stored}


def test_load():
    assert pe.load(_Files(), 'ds') == {}
    assert pe.load(_Files({7: {'sequence': 'M'}}), 'ds') == {7: {'sequence': 'M'}}
