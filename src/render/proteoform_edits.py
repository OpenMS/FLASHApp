"""Edits of identified proteoforms made in the sequence view.

The identification sequence view lets the user change the sequence and the
modifications of the shown proteoform (openms-insight's ``edit_identifier``).
Every edit arrives as a dict with an ``action`` and a fresh ``nonce``:

- ``'edit'``: start from the search result, its modifications each placed on
  one residue (the first of its range; the view then marks where the
  fragments allow it);
- ``'update'``: the edited ``sequence``, ``start_position``, ``end_position``
  and ``modifications``;
- ``'reset'``: back to the search result.

The edits of a dataset are kept per proteoform in the workspace cache (under
``RESULT``), so they outlast a rerun and a reload. Streamlit-free.
"""

import re

RESULT = 'proteoform_edits'
NONCE = '_nonce'
_RESIDUES = re.compile(r'^[A-Z]+$')


def load(file_manager, dataset_id):
    """The edits of a dataset: {proteoform index: edit}, plus the last nonce."""
    if file_manager.result_exists(dataset_id, RESULT):
        edits = file_manager.get_results(dataset_id, [RESULT])[RESULT]
        if isinstance(edits, dict):
            return dict(edits)
    return {}


def _clean_modifications(modifications, length):
    """Modifications of an update, each on one residue (or terminus) of a
    sequence of `length` residues."""
    cleaned = []
    for m in modifications or []:
        try:
            start = int(m['start'])
            mass = float(m['mass_diff'])
        except (KeyError, TypeError, ValueError):
            continue
        if mass == 0 or start < -1 or start > length:
            continue
        cleaned.append({
            'start': start, 'end': start, 'mass_diff': mass,
            'labels': str(m.get('labels') or ''),
        })
    return sorted(cleaned, key=lambda m: m['start'])


def from_search(entry):
    """The edit that starts from a search result `entry` (the sequence view
    entry of the proteoform)."""
    sequence = ''.join(entry.get('sequence') or [])
    modifications = []
    for m in entry.get('modifications') or []:
        if m is None or m.get('mass_diff') is None:
            continue
        modifications.append({
            'start': int(m['start']), 'end': int(m['start']),
            'mass_diff': float(m['mass_diff']), 'labels': str(m.get('labels') or ''),
        })
    return {
        'sequence': sequence,
        'start_position': _position(entry.get('start_position'), -1),
        'end_position': _position(entry.get('end_position'), -1),
        'modifications': _clean_modifications(modifications, len(sequence)),
    }


def _position(value, default):
    try:
        if value is None or value != value:
            return default
        return int(value)
    except (TypeError, ValueError):
        return default


def apply(edits, edit, proteoform, entry):
    """Apply `edit` (the value the view wrote) to the edits of a dataset for
    the shown `proteoform`, whose search result is `entry`.

    Returns the new edits, or None when there is nothing to change (no edit,
    one applied before, or an invalid one)."""
    if not isinstance(edit, dict) or proteoform is None:
        return None
    nonce = edit.get('nonce')
    if nonce is None or nonce == edits.get(NONCE):
        return None
    action = edit.get('action')
    new = dict(edits)
    new[NONCE] = nonce
    key = int(proteoform)
    if action == 'reset':
        new.pop(key, None)
    elif action == 'edit':
        if entry is None:
            return None
        new[key] = from_search(entry)
    elif action == 'update':
        sequence = str(edit.get('sequence') or '').upper()
        if not _RESIDUES.match(sequence):
            return None
        start = _position(edit.get('start_position'), 0)
        end = _position(edit.get('end_position'), len(sequence) - 1)
        if start >= len(sequence) or end >= len(sequence):
            start, end = 0, len(sequence) - 1
        new[key] = {
            'sequence': sequence, 'start_position': start, 'end_position': end,
            'modifications': _clean_modifications(edit.get('modifications'), len(sequence)),
        }
    else:
        return None
    return new


def entry(edit, search_entry, fragment_data):
    """The sequence view entry of an edited proteoform.

    `fragment_data(sequence)` gives the theoretical fragment masses and mass
    of an unmodified sequence (see ``getFragmentDataFromSeq``); the view adds
    the modifications itself. The observed mass stays that of the search, and
    the tag coverage stays while the sequence is unchanged."""
    sequence = edit['sequence']
    start = edit['start_position']
    end = edit['end_position']
    first = 0 if start < 0 else start
    last = len(sequence) - 1 if end < 0 else end
    data = dict(fragment_data(sequence[first:last + 1]))
    data['sequence'] = list(sequence)
    data['start_position'] = start
    data['end_position'] = end
    data['observed_mass'] = search_entry.get('observed_mass') if search_entry else None
    data['modifications'] = [dict(m) for m in edit['modifications']]
    data['edited'] = True
    data['annotated_fragments'] = None
    data['annotation_summary'] = None
    unchanged = search_entry is not None and list(search_entry.get('sequence') or []) == list(sequence)
    if unchanged and search_entry.get('coverage') is not None:
        data['coverage'] = list(search_entry['coverage'])
        data['maxCoverage'] = search_entry.get('maxCoverage')
    else:
        data['coverage'] = [0.0] * len(sequence)
        data['maxCoverage'] = 0
    return data
