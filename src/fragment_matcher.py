"""Fragment masses matched against a sequence without a FLASHDeconv run.

The Match fragments page (content/Tools/FragmentMatcher.py) takes a list of
observed fragment masses, pasted or read from a file (plain text, tsv, csv or
a TopFD/FLASHDeconv .msalign), and shows them in the sequence view against a
sequence the user types or pastes there. Streamlit-free: this module parses
and stores the input; the page draws it.

Other pages open the matcher with a mass list of their own through
``store_input`` (and then switch to the page), e.g. the masses of a selected
spectrum and the sequence of an identification.

A list of m/z values (a fragment spectrum's peaks) is turned into masses by
``deconvolve_mz``, FLASHDeconv's own spectral deconvolution in pyOpenMS.
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
    skipped. Text holding a single line with semicolons is a semicolon-
    separated list of masses instead. An .msalign file: the peaks of its first spectrum, and its
    PRECURSOR_MASS.

    Returns {'masses': [...], 'intensities': [...] or None,
    'precursor_mass': float or None, 'skipped': number of lines left out}.
    """
    if isinstance(text, bytes):
        text = text.decode('utf-8', errors='replace')
    lines = text.splitlines()
    if any(line.strip().upper() == 'BEGIN IONS' for line in lines):
        return _parse_msalign(lines)
    data = [line for line in lines if line.strip() and not line.lstrip().startswith('#')]
    if len(data) == 1 and ';' in data[0]:
        # All on one line: semicolons separate the masses (each optionally
        # followed by its intensity), as in '1043.5; 2511.37 18000; 980.1'
        lines = data[0].split(';')
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


def clear_input(workspace):
    """Forget the stored input; True if there was one."""
    path = input_path(workspace)
    if not path.exists():
        return False
    path.unlink()
    return True


def deconvolve_mz(mz, intensities=None, min_charge=1, max_charge=None):
    """Neutral monoisotopic masses of a fragment spectrum given as m/z peaks,
    by FLASHDeconv (pyOpenMS ``FLASHDeconvAlgorithm``, its default settings).

    ``intensities`` default to 1 per peak (FLASHDeconv scores isotope
    patterns by intensity, so real intensities give better masses).
    ``max_charge`` None keeps FLASHDeconv's default.

    Returns {'masses': [...], 'intensities': [...]} ordered by mass; empty
    lists when no isotope pattern is found.
    """
    import pyopenms as oms

    peaks = sorted(zip(mz, intensities or [1.0] * len(mz)))
    if not peaks:
        return {'masses': [], 'intensities': []}
    spectrum = oms.MSSpectrum()
    spectrum.setMSLevel(2)
    spectrum.setRT(1.0)
    spectrum.setNativeID('scan=2')
    spectrum.set_peaks(([p[0] for p in peaks], [float(p[1] or 1.0) for p in peaks]))
    precursor = oms.Precursor()
    precursor.setMZ(peaks[-1][0])
    spectrum.setPrecursors([precursor])
    # FLASHDeconv reads a run: one MS1 scan before the fragment spectrum
    survey = oms.MSSpectrum()
    survey.setMSLevel(1)
    survey.setRT(0.9)
    survey.setNativeID('scan=1')
    survey.set_peaks(([peaks[-1][0]], [1.0]))
    run = oms.MSExperiment()
    run.addSpectrum(survey)
    run.addSpectrum(spectrum)

    algorithm = oms.FLASHDeconvAlgorithm()
    params = algorithm.getParameters()
    params.setValue('SD:min_charge', int(min_charge))
    if max_charge is not None:
        params.setValue('SD:max_charge', int(max_charge))
    algorithm.setParameters(params)
    deconvolved, _ = algorithm.run(run)
    found = []
    for spec in deconvolved:
        if spec.getOriginalSpectrum().getMSLevel() != 2:
            continue
        found.extend((pg.getMonoMass(), pg.getIntensity()) for pg in spec)
    found.sort()
    return {'masses': [m for m, _ in found], 'intensities': [i for _, i in found]}


class UniProtLookupError(ValueError):
    """An accession that is not one, or a UniProt entry that cannot be had."""


def uniprot_entry(text):
    """The UniProt entry of an accession (``P0ACF0``, an isoform ``P0ACF0-2`` or a
    FASTA identifier ``sp|P0ACF0|DBHA_ECOLI``): {'accession', 'entry', 'name',
    'organism', 'sequence'}, fetched with the app's UniProt client."""
    from src.render.structure import uniprot_accession
    from src.workflow import uniprot

    text = (text or '').strip()
    # Accessions are upper case; a typed p0acf0 means P0ACF0
    if '|' not in text:
        text = text.upper()
    accession = uniprot_accession(text)
    if accession is None:
        raise UniProtLookupError(f"'{text}' is not a UniProt accession (e.g. P0ACF0 or P0ACF0-2)")
    isoform = re.search(rf'{accession}(-\d+)', text.upper())
    accession += isoform.group(1) if isoform else ''
    try:
        fasta = uniprot._get(f'/uniprotkb/{accession}.fasta', {}).text
    except uniprot.UniProtError as e:
        message = str(e)
        if '400' in message or '404' in message:
            message = f'UniProt has no entry {accession}'
        raise UniProtLookupError(message) from e
    lines = fasta.strip().splitlines()
    if not lines or not lines[0].startswith('>') or len(lines) < 2:
        raise UniProtLookupError(f'UniProt has no entry {accession}')
    header = lines[0]
    found = re.match(r'>\w+\|[^|]+\|(\S+)\s+(.*?)(?:\s+OS=(.*?))?(?:\s+OX=.*)?$', header)
    return {
        'accession': accession,
        'entry': found.group(1) if found else accession,
        'name': found.group(2) if found else header[1:],
        'organism': (found.group(3) or '') if found else '',
        'sequence': ''.join(lines[1:]).replace(' ', ''),
    }


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


# Proteoforms tried against one mass list ------------------------------------
#
# The sequence view writes the fragments it matched, with the proteoform it
# matched them against, to the fragment_matches selection. The page keeps the
# proteoforms the user adds as hypotheses, one row each, for the mass list
# they were tried on (a new list starts an empty table).

def hypotheses_path(workspace):
    return Path(workspace, 'fragment-matcher', 'hypotheses.json')


def hypothesis(matches, precursor_mass=None):
    """One row for the proteoform of a fragment_matches selection, or None.

    ``matches`` is what the sequence view wrote (``proteoform`` and
    ``fragments``, see openms-insight's ``fragment_matches_proteoform``).
    """
    if not isinstance(matches, dict) or not isinstance(matches.get('proteoform'), dict):
        return None
    proteoform = matches['proteoform']
    if not proteoform.get('sequence'):
        return None
    fragments = [f for f in matches.get('fragments') or [] if f.get('kind') in ('N', 'C')]
    theoretical = float(proteoform.get('theoreticalMass') or 0)
    tolerance = matches.get('tolerance') or {}
    row = {
        'proforma': proteoform.get('proforma') or proteoform['sequence'],
        'length': int(proteoform.get('length') or len(proteoform['sequence'])),
        'theoretical_mass': theoretical,
        # Counted as the view's fragment table does (an ion matched on two
        # masses, e.g. one isotope off, is two rows there)
        'matched_fragments': len(fragments),
        'n_terminal': sum(f.get('kind') == 'N' for f in fragments),
        'c_terminal': sum(f.get('kind') == 'C' for f in fragments),
        'cleavages': float(proteoform.get('cleavagePercentage') or 0),
        'tolerance': f"{tolerance.get('value')} {tolerance.get('unit')}" if tolerance else '',
        'delta_da': None,
        'delta_ppm': None,
    }
    if precursor_mass and theoretical:
        delta = float(precursor_mass) - theoretical
        row['delta_da'] = delta
        row['delta_ppm'] = delta / theoretical * 1e6
    return row


def load_hypotheses(workspace, data):
    """The hypotheses kept for this input (``load_input``), oldest first."""
    path = hypotheses_path(workspace)
    if data is None or not path.exists():
        return []
    try:
        stored = json.loads(path.read_text())
    except (OSError, ValueError):
        return []
    if stored.get('input') != input_stamp(data):
        return []
    return list(stored.get('rows') or [])


def _store_hypotheses(workspace, data, rows):
    path = hypotheses_path(workspace)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({'input': input_stamp(data), 'rows': rows}))
    return rows


def add_hypothesis(workspace, data, row):
    """Keep ``row``; one tried again (same ProForma and tolerance) replaces
    its earlier result in place."""
    rows = load_hypotheses(workspace, data)
    key = (row['proforma'], row.get('tolerance'))
    for i, kept in enumerate(rows):
        if (kept.get('proforma'), kept.get('tolerance')) == key:
            rows[i] = row
            break
    else:
        rows.append(row)
    return _store_hypotheses(workspace, data, rows)


def remove_hypotheses(workspace, data, indices):
    rows = [r for i, r in enumerate(load_hypotheses(workspace, data)) if i not in set(indices)]
    return _store_hypotheses(workspace, data, rows)


def best_hypothesis(rows):
    """Index of the row with the most matched fragments (ties: more
    cleavages, then the smaller |Δ mass|), or None."""
    if not rows:
        return None
    def rank(item):
        _, row = item
        delta = row.get('delta_da')
        return (row['matched_fragments'], row['cleavages'],
                -abs(delta) if delta is not None else 0)
    return max(enumerate(rows), key=rank)[0]
