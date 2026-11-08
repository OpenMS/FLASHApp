"""Fragment masses matched against a sequence without a FLASHDeconv run.

The Match fragments page (content/Tools/FragmentMatcher.py) takes a list of
observed fragment masses, pasted or read from a file (plain text, tsv, csv or
a TopFD/FLASHDeconv .msalign), and shows them in the sequence view against a
sequence the user types or pastes there. Streamlit-free: this module parses
and stores the input; the page draws it.

Other pages open the matcher with a mass list of their own through
``store_input`` (and then switch to the page), e.g. the masses of a selected
spectrum and the sequence of an identification.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

PROTON = 1.007276466812

MONOISOTOPIC = 'monoisotopic'
AVERAGE = 'average'
NEUTRAL = 'M'
PROTONATED = 'MH+'

# Starting fragment tolerance per mass type: average masses of large
# fragments are seldom closer than a few tenths of a Dalton
DEFAULT_SETTINGS = {
    MONOISOTOPIC: {'tolerance': 10, 'tolerance_unit': 'ppm'},
    AVERAGE: {'tolerance': 0.5, 'tolerance_unit': 'Da'},
}

_NUMBER = re.compile(r'^[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?$')


class MassListError(ValueError):
    """The text holds no mass the matcher can use."""


def _numbers(line):
    fields = [f for f in re.split(r'[\s,;]+', line.strip()) if f]
    values = []
    for field in fields:
        if not _NUMBER.match(field):
            return None
        values.append(float(field))
    return values


def parse_mass_list(text):
    """Masses (and intensities, when given) from pasted text or a file.

    Plain lists: one mass per line, optionally followed by its intensity
    (and more columns, which are ignored), separated by spaces, tabs, commas
    or semicolons; lines that are not numbers (a header, comments) are
    skipped. An .msalign file: the peaks of its first spectrum, and its
    PRECURSOR_MASS.

    Returns {'masses': [...], 'intensities': [...] or None,
    'precursor_mass': float or None, 'skipped': number of lines left out}.
    """
    if isinstance(text, bytes):
        text = text.decode('utf-8', errors='replace')
    lines = text.splitlines()
    if any(line.strip().upper() == 'BEGIN IONS' for line in lines):
        return _parse_msalign(lines)
    masses, intensities, skipped = [], [], 0
    for line in lines:
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        values = _numbers(line)
        if not values or values[0] <= 0:
            skipped += 1
            continue
        masses.append(values[0])
        intensities.append(values[1] if len(values) > 1 else None)
    if not masses:
        raise MassListError('No masses found: give one mass per line, optionally followed by its intensity.')
    return {
        'masses': masses,
        'intensities': intensities if all(i is not None for i in intensities) else None,
        'precursor_mass': None,
        'skipped': skipped,
    }


def _parse_msalign(lines):
    masses, intensities, precursor, inside, skipped = [], [], None, False, 0
    for line in lines:
        stripped = line.strip()
        if stripped.upper() == 'BEGIN IONS':
            inside = True
            continue
        if stripped.upper() == 'END IONS':
            if masses:
                break
            inside = False
            continue
        if not inside or not stripped:
            continue
        if '=' in stripped:
            key, _, value = stripped.partition('=')
            if key.strip().upper() == 'PRECURSOR_MASS' and precursor is None:
                try:
                    precursor = float(value.split(':')[0])
                except ValueError:
                    pass
            continue
        values = _numbers(stripped)
        if not values or values[0] <= 0:
            skipped += 1
            continue
        masses.append(values[0])
        intensities.append(values[1] if len(values) > 1 else None)
    if not masses:
        raise MassListError('The .msalign file holds no peaks.')
    return {
        'masses': masses,
        'intensities': intensities if all(i is not None for i in intensities) else None,
        'precursor_mass': precursor if precursor and precursor > 0 else None,
        'skipped': skipped,
    }


def to_neutral(masses, given_as=NEUTRAL):
    """Neutral masses from masses given as M or as MH+ (one proton added)."""
    if given_as == PROTONATED:
        return [m - PROTON for m in masses]
    return list(masses)


def matcher_input(masses, intensities=None, precursor_mass=None, mass_type=MONOISOTOPIC,
                  sequence=None, name=None):
    """The matcher's input as stored: neutral masses, sorted, with intensities."""
    if mass_type not in (MONOISOTOPIC, AVERAGE):
        raise ValueError(f'mass_type must be {MONOISOTOPIC!r} or {AVERAGE!r}')
    pairs = sorted(
        zip(masses, intensities if intensities is not None else [None] * len(masses)),
        key=lambda p: p[0],
    )
    return {
        'masses': [float(m) for m, _ in pairs],
        'intensities': None if intensities is None else [float(i) for _, i in pairs],
        'precursor_mass': float(precursor_mass) if precursor_mass else None,
        'mass_type': mass_type,
        'sequence': sequence or None,
        'name': name or None,
    }


def input_path(workspace):
    return Path(workspace, 'fragment-matcher', 'input.json')


def store_input(workspace, masses, intensities=None, precursor_mass=None,
                mass_type=MONOISOTOPIC, sequence=None, name=None):
    """Save the masses (and a sequence) the Match fragments page shows.

    For another page that opens the matcher: call this, then
    ``st.switch_page`` to the page. ``sequence`` (plain one-letter codes)
    replaces the sequence shown when given."""
    data = matcher_input(masses, intensities, precursor_mass, mass_type, sequence, name)
    path = input_path(workspace)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))
    return data


def load_input(workspace):
    """The stored input, or None."""
    path = input_path(workspace)
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return None


def spectrum_row(data):
    """The one row the sequence view matches: its masses and precursor mass."""
    return {
        'index': 0,
        'MonoMass': list(data['masses']),
        'PrecursorMass': data.get('precursor_mass') or 0.0,
    }


def input_stamp(data):
    """A short hash of the input, naming the view's cache."""
    import hashlib
    payload = json.dumps({k: data.get(k) for k in ('masses', 'precursor_mass', 'mass_type')})
    return hashlib.sha1(payload.encode()).hexdigest()[:12]
