"""Fragment ion coverage of a proteoform, drawn with the template's sequence
coverage panel (src/viz/identification/sequence_coverage.py).

Every fragment ion FLASHAnnotator matched in the proteoform's spectrum is a
bar over the residues it holds, one track per ion type: a terminal ion
(a, b, c from the N-terminus, x, y, z from the C-terminus) covers its whole
prefix or suffix, an internal fragment the stretch between its two cleavages. A fragment seen at several charges is one mark whose count
is its number of charge states. A fragment's label names the modifications it
holds. Streamlit-free.
"""

import polars as pl

# Track order: N-terminal ions, then C-terminal ions, then internal fragments
ION_ORDER = ['a', 'b', 'c', 'x', 'y', 'z', 'z+1', 'z+2']
INTERNAL = 'Internal fragments'

SCHEMA = {
    'match_id': pl.Int64, 'reference': pl.Int64, 'track': pl.String, 'start': pl.Int64,
    'end': pl.Int64, 'sequence': pl.String, 'count': pl.Int64,
    'mod_positions': pl.List(pl.Int64), 'mod_labels': pl.List(pl.String),
}


def track_of(fragment):
    """The track of a fragment: 'b ions', 'y ions', ... or the internal ones."""
    if fragment.get('kind') == 'I':
        return INTERNAL
    return f"{fragment.get('ion_type')} ions"


def track_order():
    return [f'{ion} ions' for ion in ION_ORDER] + [INTERNAL]


def _offset(entry):
    """Database position (0-based) of the proteoform's first residue."""
    for fragment in entry.get('annotated_fragments') or []:
        if fragment.get('start') is not None and fragment.get('first') is not None:
            return int(fragment['start']) - int(fragment['first'])
    start = entry.get('proteoform_start')
    return int(start) if start is not None and start >= 0 else 0


def _end(entry, offset):
    """Database position (0-based, inclusive) of the proteoform's last residue."""
    sequence = entry.get('sequence') or []
    end = entry.get('proteoform_end')
    if end is not None and end >= offset:
        return min(int(end), len(sequence) - 1)
    return len(sequence) - 1


def residues(entry):
    """The residues of the proteoform (a list of tokens), or []."""
    sequence = list(entry.get('sequence') or [])
    if not sequence:
        return []
    offset = _offset(entry)
    return sequence[offset:_end(entry, offset) + 1]


def _modification_label(modification):
    """The name of a modification when FLASHTnT names one, else its mass shift."""
    names = [n.strip() for n in str(modification.get('labels') or '').split(';') if n.strip()]
    if len(names) == 1:
        return names[0]
    mass = modification.get('mass_diff')
    return f'{float(mass):+.2f} Da' if mass is not None else 'modified'


def modifications(entry):
    """(position, label) of each modification, 1-based on the proteoform."""
    offset = _offset(entry)
    marks = []
    for modification in entry.get('modifications') or []:
        start = modification.get('start')
        if start is None:
            continue
        position = int(start) - offset + 1
        if position < 1:
            continue
        label = _modification_label(modification)
        end = modification.get('end')
        if end is not None and int(end) > int(start):
            label = f'{label} ({position}–{int(end) - offset + 1})'
        marks.append((position, label))
    return marks


def rows(proteoform_index, entry):
    """The bars of one proteoform (see the module docstring)."""
    fragments = entry.get('annotated_fragments') or []
    marks = modifications(entry)
    bars = {}
    for fragment in fragments:
        if fragment.get('first') is None or fragment.get('last') is None:
            continue
        first, last = int(fragment['first']) + 1, int(fragment['last']) + 1
        start, end = first, last
        key = (track_of(fragment), start, end)
        bar = bars.setdefault(key, {'name': fragment.get('name') or '', 'charges': set(),
                                    'held': [l for p, l in marks if first <= p <= last]})
        bar['charges'].add(fragment.get('charge'))
    records = []
    for (track, start, end), bar in bars.items():
        inside = [(p, l) for p, l in marks if start <= p <= end]
        records.append({
            'reference': int(proteoform_index), 'track': track, 'start': start, 'end': end,
            'sequence': ' '.join([bar['name']] + [f'({l})' for l in bar['held']]), 'count': len(bar['charges']),
            'mod_positions': [p for p, _ in inside], 'mod_labels': [l for _, l in inside],
        })
    return records


def coverage_rows(entries):
    """Bars of every proteoform, ``entries`` being (proteoform_index, entry)
    pairs; ``reference`` is the proteoform the sequence view shows."""
    records = []
    for index, entry in entries:
        records.extend(rows(index, entry))
    for match_id, record in enumerate(records):
        record['match_id'] = match_id
    return pl.DataFrame(records, schema=SCHEMA)


def from_sequence_data(sequences):
    """``coverage_rows`` of a sequence data frame (see sequence_data_store)."""
    names = sequences.collect_schema().names()
    # The Viewer's sequence frame names the proteoform's ends as
    # ProteinSequenceView reads them
    sequences = sequences.rename({
        k: v for k, v in (('start_position', 'proteoform_start'), ('end_position', 'proteoform_end'))
        if k in names and v not in names
    })
    columns = ['proteoform_index', 'sequence', 'proteoform_start', 'proteoform_end',
               'modifications', 'annotated_fragments']
    frame = sequences.select([c for c in columns if c in sequences.collect_schema().names()])
    if isinstance(frame, pl.LazyFrame):
        frame = frame.collect()
    return coverage_rows((row['proteoform_index'], row) for row in frame.iter_rows(named=True))


# Charge states per fragment: few, so four shades of the matched-fragment blue
COUNT_BINS = [(1, '1', '#c7d3fb'), (2, '2', '#8ea6f5'), (3, '3', '#5c7ff0'),
              (float('inf'), '4+', '#3157e9')]
NOTES = ['A terminal ion covers its whole prefix (a, b, c) or suffix (x, y, z); '
         'an internal fragment the residues between its two cleavages.']
