"""Mass shifts of the identified proteoforms, and the Unimod modifications they match.

FLASHTnT reports, for every PrSM, the mass shifts it placed on the sequence
(``ModMass``, several separated by ``;``): variable modifications and masses
the search could not explain. Over a whole run the same shift recurs on many
PrSMs (an N-terminal acetylation, an oxidation, a sodium adduct); this module
groups the shifts, matches each group to the Unimod modifications of that mass
(pyOpenMS ModificationsDB) and keeps the PrSMs of each group, so the Viewer can
show which modifications a sample carries and which proteoforms carry them.

Intact masses are less precise than peptide masses, and a proteoform mass can
be off by one or two 13C peaks; shifts are therefore grouped and matched with
wider tolerances than for peptides, and a match after an isotope offset is
named as such.

Streamlit-free, so tests call it.
"""
from __future__ import annotations

import numpy as np
import polars as pl
from pyopenms import ModificationsDB, ResidueModification

C13_DELTA = 1.0033548378

# Sorted shifts further apart than this start a new group
CLUSTER_GAP = 0.05
# Unimod matches must be this close to the group's median shift
MATCH_TOLERANCE = 0.05
# Histogram bin width
BIN_WIDTH = 0.05
# Isotope offsets tried besides the shift itself
ISOTOPE_OFFSETS = (-1, 1)
# Groups labelled in the histogram, most PrSMs first
LABELLED_PEAKS = 10

# Shown first when several modifications share a mass
COMMON_MODIFICATIONS = [
    'Acetyl', 'Oxidation', 'Phospho', 'Methyl', 'Dimethyl', 'Trimethyl', 'Deamidated',
    'Amidated', 'Formyl', 'Carbamidomethyl', 'Sodium', 'Potassium', 'Cation:Na', 'Cation:K',
    'Dehydrated', 'Ammonia-loss', 'Gln->pyro-Glu', 'Glu->pyro-Glu', 'Dioxidation', 'Trioxidation',
    'Carbamyl', 'Nitro', 'Hydroxylation', 'Succinyl', 'Malonyl', 'Crotonyl', 'GG', 'Myristoyl',
    'Palmitoyl', 'Farnesyl', 'Lipoyl', 'Biotin', 'HexNAc', 'Hex', 'Sulfo', 'Disulfide',
]

SHIFT_SCHEMA = {
    'shift_id': pl.Int64, 'shift': pl.Float64, 'prsms': pl.Int64, 'proteoforms': pl.Int64,
    'interpretation': pl.String, 'unimod': pl.Float64, 'deviation': pl.Float64,
    'residues': pl.String, 'other': pl.String,
}
PRSM_SCHEMA = {
    'row_id': pl.Int64, 'shift_id': pl.Int64, 'shift': pl.Float64, 'index': pl.Int64,
    'scan_idx': pl.Int64, 'Scan': pl.Int64, 'proteoform_id': pl.Int64,
    'accession': pl.String, 'RetentionTime': pl.Float64, 'Score': pl.Float64,
    'site': pl.String, 'Proforma': pl.String,
}
HISTOGRAM_SCHEMA = {
    'delta': pl.Float64, 'prsms': pl.Int64, 'shift_id': pl.Int64,
    'peak': pl.Boolean, 'label': pl.String,
}


def _masses(value) -> list[float]:
    """The shifts of one PrSM: ``ModMass`` as a number or as ``;``-separated text."""
    if value is None:
        return []
    if isinstance(value, (int, float)):
        return [] if np.isnan(value) else [float(value)]
    masses = []
    for part in str(value).split(';'):
        try:
            mass = float(part)
        except ValueError:
            continue
        if not np.isnan(mass):
            masses.append(mass)
    return masses


def _positions(value) -> list:
    if value is None:
        return []
    if isinstance(value, (int, float)):
        return [] if np.isnan(value) else [int(value)]
    out = []
    for part in str(value).split(';'):
        try:
            out.append(int(float(part)))
        except ValueError:
            out.append(None)
    return out


def _site(start, end) -> str:
    """Residue (or range, when FLASHTnT could not localize the shift) in the
    protein's numbering, as FLASHTnT reports it (1-based)."""
    if start is None:
        return ''
    if end is None or end == start:
        return str(start)
    return f'{start}–{end}'


def prsm_shifts(hits: pl.DataFrame) -> pl.DataFrame:
    """One row per mass shift of a PrSM (PRSM_SCHEMA without ids), from the
    hierarchy PrSMs (``hierarchy_data.prsms``)."""
    rows = []
    columns = set(hits.columns)
    for hit in hits.iter_rows(named=True):
        masses = _masses(hit.get('ModMass'))
        if not masses:
            continue
        starts = _positions(hit.get('ModStart')) if 'ModStart' in columns else []
        ends = _positions(hit.get('ModEnd')) if 'ModEnd' in columns else []
        for k, mass in enumerate(masses):
            start = starts[k] if k < len(starts) else None
            end = ends[k] if k < len(ends) else None
            rows.append({
                'shift': mass,
                'index': hit.get('index'),
                'scan_idx': hit.get('scan_idx'),
                'Scan': hit.get('Scan'),
                'proteoform_id': hit.get('proteoform_id'),
                'accession': hit.get('accession') or '',
                'RetentionTime': hit.get('RetentionTime'),
                'Score': hit.get('Score'),
                'site': _site(start, end),
                'Proforma': hit.get('Proforma') or '',
            })
    schema = {k: v for k, v in PRSM_SCHEMA.items() if k not in ('row_id', 'shift_id')}
    return pl.DataFrame(rows, schema=schema)


def unimod_matches(shift: float, tolerance: float = MATCH_TOLERANCE) -> list[dict]:
    """Unimod modifications within ``tolerance`` of ``shift``, one per name,
    common ones first, then the closest."""
    anywhere = ResidueModification.TermSpecificity.ANYWHERE
    db = ModificationsDB()
    try:
        # pyOpenMS 3.6 returns the ids
        full_ids = db.searchModificationsByDiffMonoMass(float(shift), float(tolerance), '', anywhere)
    except TypeError:
        # pyOpenMS 3.5 and earlier fill a list passed in
        full_ids = []
        db.searchModificationsByDiffMonoMass(full_ids, float(shift), float(tolerance), '', anywhere)
    by_name = {}
    for full_id in full_ids:
        full_id = full_id.decode() if isinstance(full_id, bytes) else str(full_id)
        mod = db.getModification(full_id)
        name = mod.getId()
        # Run-time placeholders and PSI-MOD duplicates of Unimod entries
        if not name or name.startswith(('CUSTOM', '[', 'MOD:')):
            continue
        entry = by_name.setdefault(name, {'name': name, 'delta': mod.getDiffMonoMass(), 'residues': set()})
        term = mod.getTermSpecificity()
        origin = mod.getOrigin()
        if term == ResidueModification.TermSpecificity.N_TERM:
            entry['residues'].add('N-term')
        elif term == ResidueModification.TermSpecificity.C_TERM:
            entry['residues'].add('C-term')
        elif origin and origin != 'X':
            entry['residues'].add(origin)
    matches = []
    for entry in by_name.values():
        residues = entry['residues']
        entry['residues'] = ''.join(sorted(r for r in residues if len(r) == 1)) + ''.join(
            f' {r}' for r in sorted(r for r in residues if len(r) > 1)
        )
        entry['deviation'] = shift - entry['delta']
        matches.append(entry)
    rank = {name: i for i, name in enumerate(COMMON_MODIFICATIONS)}
    matches.sort(key=lambda m: (rank.get(m['name'], len(rank)), abs(m['deviation']), m['name']))
    return matches


def _isotope_label(n: int) -> str:
    peaks = 'peak' if abs(n) == 1 else 'peaks'
    return f'{n:+d} isotope {peaks}'


def interpret(shift: float) -> dict:
    """Best interpretation of a shift. Common modifications come first, matched
    to the shift itself, then after one missed or extra isotope peak; then any
    other Unimod entry of that mass, likewise. Unimod holds an entry for
    almost every mass (amino acid substitutions), so a common modification
    one isotope peak off is the likelier reading of an intact mass."""
    common = set(COMMON_MODIFICATIONS)
    candidates = []
    for rank_offset, n in enumerate((0, *ISOTOPE_OFFSETS)):
        for match in unimod_matches(shift - n * C13_DELTA):
            candidates.append((0 if match['name'] in common else 1, rank_offset, n, match))
    if not candidates:
        return {'interpretation': '', 'unimod': None, 'deviation': None, 'residues': '', 'other': ''}
    candidates.sort(key=lambda c: (c[0], c[1], abs(c[3]['deviation'])))
    _, _, offset, best = candidates[0]
    name = best['name'] if offset == 0 else f"{best['name']} ({_isotope_label(offset)})"
    unimod = best['delta'] + offset * C13_DELTA
    others = []
    for _, _, n, match in candidates[1:]:
        label = match['name'] if n == 0 else f"{match['name']} ({_isotope_label(n)})"
        if label not in others and len(others) < 5:
            others.append(label)
    return {
        'interpretation': name,
        'unimod': round(unimod, 4),
        'deviation': round((shift - unimod) * 1000, 1),
        'residues': best['residues'].strip(),
        'other': ', '.join(others),
    }


def _groups(values: np.ndarray) -> list[np.ndarray]:
    """Indices of ``values`` grouped into shifts: sorted values further apart
    than CLUSTER_GAP start a new group."""
    if len(values) == 0:
        return []
    order = np.argsort(values, kind='stable')
    breaks = np.where(np.diff(values[order]) > CLUSTER_GAP)[0] + 1
    return [group for group in np.split(order, breaks) if len(group)]


def summary(shifts: pl.DataFrame) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """Groups of shifts, their histogram and their PrSMs, linked by ``shift_id``.

    Returns:
        groups: SHIFT_SCHEMA, one row per group, most PrSMs first.
        histogram: one row per non-empty bin, the tallest bin of each of the
            LABELLED_PEAKS largest groups marked (``peak``) and labelled.
        prsms: PRSM_SCHEMA, the input rows with their group and a row id.
    """
    if not shifts.height:
        return (
            pl.DataFrame(schema=SHIFT_SCHEMA),
            pl.DataFrame(schema=HISTOGRAM_SCHEMA),
            pl.DataFrame(schema=PRSM_SCHEMA),
        )
    values = shifts['shift'].to_numpy().astype(float)
    proteoforms = shifts['proteoform_id'].to_list()
    rows = []
    member_of = np.full(len(values), -1, dtype=np.int64)
    for members in _groups(values):
        median = float(np.median(values[members]))
        rows.append({
            'members': members,
            'shift': round(median, 4),
            'prsms': len(members),
            'proteoforms': len({proteoforms[i] for i in members}),
            **interpret(median),
        })
    rows.sort(key=lambda r: (-r['prsms'], r['shift']))
    for shift_id, row in enumerate(rows):
        row['shift_id'] = shift_id
        member_of[row.pop('members')] = shift_id
    groups = pl.DataFrame(rows, schema=SHIFT_SCHEMA)

    prsms = shifts.with_columns(pl.Series('shift_id', member_of)).with_row_index('row_id').with_columns(
        pl.col('row_id').cast(pl.Int64)
    ).select(list(PRSM_SCHEMA)).cast(PRSM_SCHEMA)

    bins = np.round(values / BIN_WIDTH).astype(np.int64)
    histogram = (
        pl.DataFrame({'bin': bins, 'shift_id': member_of})
        .group_by('bin')
        .agg(pl.len().cast(pl.Int64).alias('prsms'), pl.col('shift_id').max())
        .sort('bin')
        .with_columns((pl.col('bin') * BIN_WIDTH).alias('delta'))
    )
    labels = {}
    for row in rows[:LABELLED_PEAKS]:
        name = row['interpretation'] or '?'
        labels[row['shift_id']] = f"{row['shift']:+.2f} {name}"
    peak_bins = (
        histogram.filter(pl.col('shift_id').is_in(list(labels)))
        .sort(['prsms', 'bin'], descending=[True, False])
        .unique('shift_id', keep='first')
        .select('bin')
    )['bin'].to_list()
    histogram = histogram.with_columns(
        pl.col('bin').is_in(peak_bins).alias('peak'),
        pl.when(pl.col('bin').is_in(peak_bins))
        .then(pl.col('shift_id').replace_strict(labels, default='', return_dtype=pl.String))
        .otherwise(pl.lit(''))
        .alias('label'),
    ).select(list(HISTOGRAM_SCHEMA)).cast(HISTOGRAM_SCHEMA)
    return groups, histogram, prsms
