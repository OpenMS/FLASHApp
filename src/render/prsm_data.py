"""Input frame of the PrSM mass table.

Every deconvolved mass of the identified spectrum, with the fragment ion it
explains (or none), the theoretical mass and the mass error, so the list can
be filtered to all, matched or not matched masses. This module builds that
list for FLASHTnT results.

Masses come from the spectrum of each proteoform (the `sequence_view` and
`mass_table` results); fragments are matched with the rule of the Sequence
View and the fragment error plot (validation_data.match_fragments). Where
FLASHAnnotator envelopes are stored, a deconvolved mass counts as matched
by an annotated fragment whose observed mass lies within the tolerance.
Like insight_data, the module is free of Streamlit and openms-insight.
"""

import numpy as np
import polars as pl

from src.render import insight_data, validation_data

PRSM_MASS_SCHEMA = {
    'proteoform_index': pl.Int64,
    'mass_idx': pl.Int64,
    'mass_id': pl.Int64,
    'MonoMass': pl.Float64,
    'SumIntensity': pl.Float64,
    'charges': pl.String,
    'QScore': pl.Float64,
    'matched': pl.Boolean,
    'status': pl.String,
    'ions': pl.String,
    'position': pl.Int64,
    'theoretical_mass': pl.Float64,
    'mass_error': pl.Float64,
    'ppm_error': pl.Float64,
}

MATCHED = 'Matched'
NOT_MATCHED = 'Not matched'


def _settings(file_manager, dataset_id):
    settings = {}
    if insight_data.result_path(file_manager, dataset_id, 'settings') is not None:
        settings = file_manager.get_results(dataset_id, ['settings'])['settings'] or {}
    tolerance = float(settings.get('tolerance') or 10.0)
    ions = validation_data.PREFIX_IONS + validation_data.SUFFIX_IONS
    ion_types = [i for i in (settings.get('ion_types') or ['b', 'c', 'y', 'z']) if i in ions]
    return tolerance, ion_types


def _nearest(masses, value, tolerance_ppm):
    """Index of the mass in `masses` nearest to `value` within the tolerance,
    or None."""
    if value is None or masses.size == 0:
        return None
    nearest = int(np.argmin(np.abs(masses - value)))
    if abs(masses[nearest] - value) <= abs(value) * tolerance_ppm * 1e-6:
        return nearest
    return None


def mass_matches(entry, masses, tolerance_ppm, ion_types):
    """Matches of one proteoform, keyed by the index of the deconvolved mass.

    Returns {mass index: [match, ...]} with the match dicts of
    validation_data (fragment, residue, theoretical_mass, observed_mass,
    ppm_error)."""
    values = [float(m) for m in (masses if masses is not None else []) if m is not None]
    masses = np.asarray(values, dtype=float)
    annotated = entry.get('annotated_fragments')
    by_mass = {}
    if annotated:
        for match in validation_data.annotated_fragment_errors(annotated):
            index = _nearest(masses, match.get('observed_mass'), tolerance_ppm)
            if index is not None and match.get('ppm_error') is not None:
                by_mass.setdefault(index, []).append(match)
        return by_mass
    # match_fragments returns the observed values themselves, so the index is
    # an exact lookup
    positions = {}
    for index, mass in enumerate(masses):
        positions.setdefault(float(mass), index)
    for match in validation_data.match_fragments(entry, values, tolerance_ppm, ion_types):
        index = positions.get(float(match['observed_mass']))
        if index is not None:
            by_mass.setdefault(index, []).append(match)
    return by_mass


def _row(proteoform_index, mass_idx, mass, details, matches):
    row = {
        'proteoform_index': proteoform_index,
        'mass_idx': mass_idx,
        'mass_id': mass_idx + 1,
        'MonoMass': float(mass),
        'SumIntensity': details.get('SumIntensity'),
        'charges': details.get('charges'),
        'QScore': details.get('QScore'),
        'matched': bool(matches),
        'status': MATCHED if matches else NOT_MATCHED,
        'ions': None, 'position': None, 'theoretical_mass': None,
        'mass_error': None, 'ppm_error': None,
    }
    if matches:
        # Every ion the mass explains is named; the numbers are the closest one's
        best = min(matches, key=lambda m: abs(m['ppm_error']))
        names = []
        for match in sorted(matches, key=lambda m: abs(m['ppm_error'])):
            if match['fragment'] not in names:
                names.append(match['fragment'])
        residue = best.get('residue')
        row.update({
            'ions': ', '.join(str(n) for n in names),
            # 1-based residue next to the bond, as in the Sequence View
            'position': None if residue is None else int(residue) + 1,
            'theoretical_mass': best['theoretical_mass'],
            'mass_error': float(mass) - best['theoretical_mass'],
            'ppm_error': best['ppm_error'],
        })
    return row


def _mass_details(file_manager, dataset_id):
    """{(scan_idx, mass_idx): intensity, charge range and QScore} from the
    mass table, when the results have it."""
    if insight_data.result_path(file_manager, dataset_id, 'mass_table') is None:
        return {}
    table = insight_data.mass_table(file_manager, dataset_id).select(
        'scan_idx', 'mass_idx', 'SumIntensity', 'MinCharges', 'MaxCharges', 'QScore'
    ).collect()
    details = {}
    for row in table.iter_rows(named=True):
        low, high = row['MinCharges'], row['MaxCharges']
        if low is None or high is None:
            charges = None
        elif low == high:
            charges = str(low)
        else:
            charges = f'{low}–{high}'
        details[(row['scan_idx'], row['mass_idx'])] = {
            'SumIntensity': None if row['SumIntensity'] is None else float(row['SumIntensity']),
            'charges': charges,
            'QScore': row['QScore'],
        }
    return details


def prsm_masses(file_manager, dataset_id):
    """One row per deconvolved mass of the spectrum of every proteoform, with
    the fragment ions it matches (PRSM_MASS_SCHEMA).

    `mass_idx` is the mass's index in the spectrum, the value of the viewer's
    `mass` selection; `mass_id` is the same, 1-based, as shown in the table."""
    tolerance, ion_types = _settings(file_manager, dataset_id)
    details = _mass_details(file_manager, dataset_id)
    sequences = insight_data.proteoform_sequences(file_manager, dataset_id).collect()
    rows = []
    for entry in sequences.iter_rows(named=True):
        entry.setdefault('proteoform_start', entry.get('start_position'))
        entry.setdefault('proteoform_end', entry.get('end_position'))
        masses = entry.get('MonoMass') or []
        matches = mass_matches(entry, masses, tolerance, ion_types)
        scan_idx = entry.get('scan_idx')
        for mass_idx, mass in enumerate(masses):
            if mass is None:
                continue
            rows.append(_row(
                entry['proteoform_index'], mass_idx, mass,
                details.get((scan_idx, mass_idx), {}), matches.get(mass_idx),
            ))
    return (
        pl.DataFrame(rows, schema=PRSM_MASS_SCHEMA)
        .sort(['proteoform_index', 'mass_idx'])
        .lazy()
    )
