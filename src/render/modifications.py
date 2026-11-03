"""Modification library for the proteoform sequence, from pyOpenMS ModificationsDB.

ModificationsDB holds the Unimod modifications (glycans such as
``Hex(5)HexNAc(4) (N)`` among them). This module lists and searches them by
site and applies a set of fixed modifications to an ``AASequence``: a residue
modification to every residue of its kind, a terminal one to the proteoform's
terminus. It is Streamlit-free.
"""

from functools import lru_cache

from pyopenms import AASequence, ModificationsDB, ResidueModification

N_TERM = 'N-term'
C_TERM = 'C-term'

_TERM = ResidueModification.TermSpecificity
_N_TERMINAL = (_TERM.N_TERM, _TERM.PROTEIN_N_TERM)
_C_TERMINAL = (_TERM.C_TERM, _TERM.PROTEIN_C_TERM)

# Fixed modifications the former Sequence Input page offered before the library,
# by stored name: (residue, mass shift)
LEGACY_FIXED_MODS = {
    'Carbamidomethyl (+57)': ('C', 57.021464),
    'Carboxymethyl (+58)': ('C', 58.005479),
    'Xlink:Disulfide (-1 per C)': ('C', -1.007825),
    'L-methionine sulfoxide (+16)': ('M', 15.994915),
    'L-methionine sulfone (+32)': ('M', 31.989829),
}


def _text(value):
    """pyOpenMS returns some names as bytes."""
    if not isinstance(value, bytes):
        return str(value)
    try:
        return value.decode('utf-8')
    except UnicodeDecodeError:
        return value.decode('latin-1')  # e.g. the "®" in Unimod's TMT names


def _site(mod):
    """'N-term' / 'C-term' for terminal modifications, else the residue."""
    term = mod.getTermSpecificity()
    if term in _N_TERMINAL:
        return N_TERM
    if term in _C_TERMINAL:
        return C_TERM
    return _text(mod.getOrigin())


@lru_cache(maxsize=1)
def library():
    """Every modification as a dict, sorted by name then site.

    Keys: ``id`` (the full id, e.g. ``'Carbamidomethyl (C)'``), ``name``,
    ``description``, ``site`` (a residue letter, ``'N-term'`` or ``'C-term'``),
    ``origin`` (the residue a terminal modification needs, ``'X'`` for any),
    ``protein_terminal``, ``mass`` (monoisotopic mass shift) and ``unimod``.
    """
    db = ModificationsDB()
    entries = {}
    for i in range(db.getNumberOfModifications()):
        mod = db.getModification(i)
        full_id = _text(mod.getFullId())
        if full_id in entries:
            continue
        entries[full_id] = {
            'id': full_id,
            'name': _text(mod.getId()),
            'description': _text(mod.getFullName()),
            'site': _site(mod),
            'origin': _text(mod.getOrigin()),
            'protein_terminal': mod.getTermSpecificity() in (_TERM.PROTEIN_N_TERM, _TERM.PROTEIN_C_TERM),
            'mass': mod.getDiffMonoMass(),
            'unimod': _text(mod.getUniModAccession()),
        }
    return tuple(sorted(entries.values(), key=lambda e: (e['name'].lower(), e['site'], e['id'])))


def by_id(mod_id):
    """The library entry with this full id, or None."""
    for entry in library():
        if entry['id'] == mod_id:
            return entry
    return None


def search(query='', site=None, limit=None):
    """Library entries at ``site`` (a residue, 'N-term', 'C-term' or None for
    all) whose id, description or Unimod accession contains ``query`` (case
    insensitive). A numeric query matches mass shifts within 0.01 Da."""
    query = (query or '').strip()
    mass = None
    try:
        mass = float(query)
    except ValueError:
        pass
    needle = query.lower()
    hits = []
    for entry in library():
        if site is not None and entry['site'] != site:
            continue
        if mass is not None:
            if abs(entry['mass'] - mass) > 0.01:
                continue
        elif needle and not any(
            needle in text.lower() for text in (entry['id'], entry['description'], entry['unimod'])
        ):
            continue
        hits.append(entry)
        if limit is not None and len(hits) >= limit:
            break
    return hits


def legacy_ids(fixed_mod_cysteine=None, fixed_mod_methionine=None):
    """Library ids for the Cys/Met fixed modifications saved by earlier
    versions of the former Sequence Input page."""
    db = ModificationsDB()
    ids = []
    for name in (fixed_mod_cysteine, fixed_mod_methionine):
        if name not in LEGACY_FIXED_MODS:
            continue
        residue, mass = LEGACY_FIXED_MODS[name]
        mod = db.getBestModificationByDiffMonoMass(mass, 0.001, residue, _TERM.ANYWHERE)
        if mod is not None:
            ids.append(_text(mod.getFullId()))
    return ids


def fixed_ids(stored, picked=None):
    """Fixed modification ids for a saved sequence (the former Sequence Input page's
    record): ``picked`` (ids chosen in the sequence view's dialog, unknown ones
    dropped) when not None, else the saved ids, else those of the Cys/Met
    names saved before the library."""
    if picked is not None:
        return [m for m in picked if by_id(m) is not None]
    saved = stored.get('fixed_modifications')
    if saved is not None:
        return list(saved)
    return legacy_ids(stored.get('fixed_mod_cysteine'), stored.get('fixed_mod_methionine'))


def sites_in(sequence, entry):
    """0-based residue indices of ``sequence`` (a plain one-letter string) that
    ``entry`` modifies; for a terminal modification the terminal residue."""
    if not sequence:
        return []
    if entry['site'] == N_TERM:
        return [0] if entry['origin'] in ('X', sequence[0]) else []
    if entry['site'] == C_TERM:
        last = len(sequence) - 1
        return [last] if entry['origin'] in ('X', sequence[last]) else []
    return [i for i, aa in enumerate(sequence) if aa == entry['site']]


def apply_fixed(protein, mod_ids):
    """Apply the fixed modifications ``mod_ids`` to ``protein`` (an AASequence,
    changed in place).

    A residue modification goes on every residue of its kind that carries no
    modification yet, so the first of two modifications for one residue wins.
    Returns ``(protein, residues, applied)``: the residue letters that carry a
    fixed modification, and per id ``{'id', 'site', 'mass', 'count'}`` with the
    number of sites it was applied to.
    """
    db = ModificationsDB()
    plain = protein.toUnmodifiedString()
    residues = []
    applied = []
    for mod_id in mod_ids:
        entry = by_id(mod_id)
        if entry is None:
            continue
        mod = db.getModification(mod_id)
        count = 0
        if entry['site'] == N_TERM:
            if sites_in(plain, entry) and not protein.hasNTerminalModification():
                protein.setNTerminalModification(mod)
                count = 1
        elif entry['site'] == C_TERM:
            if sites_in(plain, entry) and not protein.hasCTerminalModification():
                protein.setCTerminalModification(mod)
                count = 1
        else:
            for index in sites_in(plain, entry):
                if protein.getResidue(index).isModified():
                    continue
                protein.setModification(index, mod)
                count += 1
            if count and entry['site'] not in residues:
                residues.append(entry['site'])
        applied.append({'id': mod_id, 'site': entry['site'], 'mass': entry['mass'], 'count': count})
    return protein, residues, applied


def modified_mass(sequence, mod_ids):
    """Monoisotopic mass of ``sequence`` with the fixed modifications applied."""
    protein, _, _ = apply_fixed(AASequence.fromString(sequence), mod_ids)
    return protein.getMonoWeight()
