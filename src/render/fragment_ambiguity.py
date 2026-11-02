"""Theoretical fragment masses around modifications FLASHTnT could not localize.

FLASHTnT reports a modification it could not place on one residue as a range
(ModStart..ModEnd): the mass shift sits on one of those residues. The parser
(src/render/sequence.py, getFragmentDataFromSeq) stores, for every fragment
whose cleavage falls in that range, both the mass with and without the shift.
Matching both lets the Sequence View show fragments inside the range, often
contradicting each other (one with the shift before a cleavage, another
without it after) and the range itself, which FLASHTnT would have narrowed
had such a fragment been evidence.

A fragment says nothing about where the shift is only when its cleavage
falls strictly inside the range, so those fragments get no theoretical mass.
A fragment that covers the whole range carries the shift. Kept free of
Streamlit and pyOpenMS.
"""

import polars as pl

PREFIX_IONS = ('a', 'b', 'c')
SUFFIX_IONS = ('x', 'y', 'z')
FRAGMENT_COLUMNS = [f'fragment_masses_{ion}' for ion in PREFIX_IONS + SUFFIX_IONS]


def _ranges(modifications, start_position):
    """Ranges (1-based residues of the proteoform) of the modifications that
    span more than one residue. `start` / `end` of a modification are 0-based
    in the database sequence, `start_position` the proteoform's first residue
    there (negative when FLASHTnT did not determine it)."""
    offset = max(start_position if start_position is not None else 0, 0)
    ranges = []
    for modification in modifications or []:
        start, end = modification.get('start'), modification.get('end')
        if start is None or end is None or end <= start:
            continue
        ranges.append((start - offset + 1, end - offset + 1))
    return ranges


def _settled(masses, ambiguous):
    """The masses of one fragment: none when it lies in a range, else the
    last one stored (the parser appends the shifted mass after the unshifted
    one, and a fragment covering a whole range carries the shift)."""
    if ambiguous:
        return []
    if masses is None or len(masses) <= 1:
        return masses
    return [masses[-1]]


def settle_fragment_masses(entry):
    """`fragment_masses_<ion>` of one proteoform (a dict with those lists,
    `modifications` and `start_position`) without the fragments whose
    cleavage lies inside an unlocalized modification's range."""
    ranges = _ranges(entry.get('modifications'), entry.get('start_position'))
    out = {}
    for ion in PREFIX_IONS + SUFFIX_IONS:
        column = f'fragment_masses_{ion}'
        fragments = entry.get(column)
        if not ranges or fragments is None:
            out[column] = fragments
            continue
        settled = []
        length = len(fragments)
        for index, masses in enumerate(fragments):
            if ion in PREFIX_IONS:
                # Residues 1..last; the shift is on one of start..end
                last = index + 1
                ambiguous = any(start <= last < end for start, end in ranges)
            else:
                # Residues first..length
                first = length - index
                ambiguous = any(start < first <= end for start, end in ranges)
            settled.append(_settled(masses, ambiguous))
        out[column] = settled
    return out


def settle(lf):
    """Apply settle_fragment_masses to every row of a frame of sequence data
    (the columns of insight_data.proteoform_sequences)."""
    schema = lf.collect_schema()
    columns = [c for c in FRAGMENT_COLUMNS if c in schema]
    if not columns or 'modifications' not in schema or 'start_position' not in schema:
        return lf
    dtype = pl.Struct({c: schema[c] for c in columns})
    return lf.with_columns(
        pl.struct(columns + ['modifications', 'start_position'])
        .map_elements(lambda entry: {c: settle_fragment_masses(entry)[c] for c in columns},
                      return_dtype=dtype)
        .alias('_settled')
    ).drop(columns).unnest('_settled')
