"""Input frames for the PrSM validation views.

Three views check an identification or a deconvolved mass by eye, as TopMSV
and LcMsSpectator do:

- fragment_errors: the mass error (ppm) of every matched fragment of a
  proteoform. A drift or a band at +-1 Da / mass shows a wrong charge or
  monoisotopic mass that the coverage map hides.
- isotope_envelopes: the observed isotope envelope of a deconvolved mass
  (signal peaks of all charges, summed per isotope) against the averagine
  envelope FLASHDeconv scored it with.
- precursor_xics: the precursor mass of a fragmentation spectrum traced over
  the MS1 spectra, one trace per charge.

All three read only results FLASHDeconv and FLASHTnT already wrote (see
src/parse/); pyOpenMS supplies the averagine envelopes. Like insight_data,
the module is free of Streamlit and openms-insight.
"""

import numpy as np
import polars as pl

from src.render import insight_data

PROTON_MASS = 1.007276466621
# Mass difference of 13C and 12C, the isotope spacing FLASHDeconv uses
ISOTOPE_SPACING = 1.0033548378

# Precursor masses of the MS1 spectra within this tolerance count as the same
# species for the XIC
XIC_TOLERANCE_PPM = 10.0
# Isotopes of the theoretical envelope; the tail beyond is below 0.1 %
MAX_ISOTOPES = 200
MIN_RELATIVE_ABUNDANCE = 1e-3

PREFIX_IONS = ('a', 'b', 'c')
SUFFIX_IONS = ('x', 'y', 'z')


# --- Fragment mass errors ---------------------------------------------------

FRAGMENT_ERROR_SCHEMA = {
    'proteoform_index': pl.Int64,
    'ion_type': pl.String,
    'fragment': pl.String,
    'residue': pl.Int64,
    'theoretical_mass': pl.Float64,
    'observed_mass': pl.Float64,
    'ppm_error': pl.Float64,
}


def _sorted(values):
    masses = np.asarray([v for v in (values or []) if v is not None], dtype=float)
    masses.sort()
    return masses


def match_fragments(entry, observed_masses, tolerance_ppm, ion_types):
    """Matched fragments of one proteoform, as the Sequence View matches them:
    every theoretical fragment mass of a selected ion type against the
    deconvolved masses of the spectrum within `tolerance_ppm`.

    `entry` holds the sequence data of the proteoform (src/parse/tnt.py):
    `fragment_masses_<ion>` per fragment length, `proteoform_start`/`_end`
    (0-based in the database sequence, negative when not determined).
    Returns a list of dicts with the FRAGMENT_ERROR_SCHEMA fields (but
    proteoform_index); `residue` is the database residue next to the bond.
    """
    observed = _sorted(observed_masses)
    if observed.size == 0:
        return []
    sequence = entry.get('sequence') or []
    start_reported = entry.get('proteoform_start', entry.get('start_position', -1))
    end_reported = entry.get('proteoform_end', entry.get('end_position', -1))
    start_reported = -1 if start_reported is None else int(start_reported)
    end_reported = -1 if end_reported is None else int(end_reported)
    start = max(start_reported, 0)
    end = end_reported if end_reported >= 0 else len(sequence) - 1

    rows = []
    for ion in ion_types:
        # Fragments of a terminus FLASHTnT did not determine are not matched
        if ion in PREFIX_IONS and start_reported < 0:
            continue
        if ion in SUFFIX_IONS and end_reported < 0:
            continue
        for index, masses in enumerate(entry.get(f'fragment_masses_{ion}') or []):
            for theoretical in masses or []:
                if theoretical is None or theoretical <= 0:
                    continue
                window = theoretical * tolerance_ppm * 1e-6
                lo = np.searchsorted(observed, theoretical - window, side='left')
                hi = np.searchsorted(observed, theoretical + window, side='right')
                for mass in observed[lo:hi]:
                    rows.append({
                        'ion_type': ion,
                        'fragment': f'{ion}{index + 1}',
                        'residue': start + index if ion in PREFIX_IONS else end - index,
                        'theoretical_mass': float(theoretical),
                        'observed_mass': float(mass),
                        'ppm_error': float((mass - theoretical) / theoretical * 1e6),
                    })
    return rows


def annotated_fragment_errors(fragments):
    """Rows of FLASHAnnotator's accepted envelopes (see src/parse/annotator.py),
    one per envelope; internal fragments are typed 'internal'."""
    rows = []
    for fragment in fragments or []:
        kind = fragment.get('kind')
        if kind == 'N':
            residue = fragment.get('end')
        else:
            residue = fragment.get('start')
        rows.append({
            'ion_type': 'internal' if kind == 'I' else fragment.get('ion_type'),
            'fragment': (
                f"{fragment.get('name')} ({fragment.get('charge')}+)"
                if fragment.get('charge') else fragment.get('name')
            ),
            'residue': residue,
            'theoretical_mass': fragment.get('mono_mass'),
            'observed_mass': fragment.get('observed_mass'),
            'ppm_error': fragment.get('ppm_error'),
        })
    return rows


def fragment_errors(file_manager, dataset_id):
    """Matched fragments of every proteoform with their mass error (ppm).

    FLASHAnnotator envelopes are used where the results have them; otherwise
    fragments are matched against the deconvolved masses with the FLASHTnT
    fragment tolerance and ion types (the `settings` result)."""
    settings = {}
    if insight_data.result_path(file_manager, dataset_id, 'settings') is not None:
        settings = file_manager.get_results(dataset_id, ['settings'])['settings'] or {}
    tolerance = float(settings.get('tolerance') or 10.0)
    ion_types = [i for i in (settings.get('ion_types') or ['b', 'c', 'y', 'z'])
                 if i in PREFIX_IONS + SUFFIX_IONS]

    sequences = insight_data.proteoform_sequences(file_manager, dataset_id).collect()
    rows = []
    for entry in sequences.iter_rows(named=True):
        # proteoform_sequences renames the positions for Insight
        entry.setdefault('proteoform_start', entry.get('start_position'))
        entry.setdefault('proteoform_end', entry.get('end_position'))
        annotated = entry.get('annotated_fragments')
        if annotated:
            matched = annotated_fragment_errors(annotated)
        else:
            matched = match_fragments(entry, entry.get('MonoMass'), tolerance, ion_types)
        for row in matched:
            row['proteoform_index'] = entry['proteoform_index']
        rows.extend(matched)
    frame = pl.DataFrame(rows, schema=FRAGMENT_ERROR_SCHEMA)
    return (
        frame.filter(pl.col('ppm_error').is_not_null())
        .select(list(FRAGMENT_ERROR_SCHEMA))
        .sort(['proteoform_index', 'ion_type', 'theoretical_mass'])
        .lazy()
    )


# --- Isotope envelopes --------------------------------------------------------

ENVELOPE_SCHEMA = {
    'scan_idx': pl.Int64,
    'mass_idx': pl.Int64,
    'side': pl.String,
    'isotope': pl.Int64,
    'mass': pl.Float64,
    'abundance': pl.Float64,
}


class Averagine:
    """Averagine isotope envelopes from pyOpenMS, cached per 10 Da."""

    BIN = 10.0

    def __init__(self):
        import pyopenms
        self._generator = pyopenms.CoarseIsotopePatternGenerator(MAX_ISOTOPES)
        self._cache = {}

    def envelope(self, mono_mass):
        """Relative abundances (maximum 1) of isotopes 0, 1, ... of a
        protein of `mono_mass`; the tail below 0.1 % is cut."""
        key = int(mono_mass // self.BIN)
        if key not in self._cache:
            # The generator takes the average mass; for averagine it is about
            # 0.06 % above the monoisotopic mass
            distribution = self._generator.estimateFromPeptideWeight(
                (key + 0.5) * self.BIN * 1.0006
            )
            abundances = np.array([p.getIntensity() for p in distribution.getContainer()])
            if abundances.size == 0 or abundances.max() <= 0:
                self._cache[key] = np.array([1.0])
            else:
                abundances = abundances / abundances.max()
                last = np.nonzero(abundances >= MIN_RELATIVE_ABUNDANCE)[0].max()
                self._cache[key] = abundances[: last + 1]
        return self._cache[key]


def _observed_envelopes(masses_with_peaks):
    """Observed envelope rows (vectorised observed_envelope) of a frame with
    scan_idx, mass_idx, MonoMass and SignalPeaks."""
    return (
        masses_with_peaks.explode('SignalPeaks')
        .drop_nulls('SignalPeaks')
        .with_columns(
            pl.col('SignalPeaks').list.get(1, null_on_oob=True).alias('mz'),
            pl.col('SignalPeaks').list.get(2, null_on_oob=True).alias('intensity'),
            pl.col('SignalPeaks').list.get(3, null_on_oob=True).alias('charge'),
        )
        .filter(pl.col('charge') > 0)
        .with_columns(
            (((pl.col('mz') - PROTON_MASS) * pl.col('charge') - pl.col('MonoMass'))
             / ISOTOPE_SPACING).round(0).cast(pl.Int64).alias('isotope')
        )
        .group_by(['scan_idx', 'mass_idx', 'MonoMass', 'isotope'])
        .agg(pl.col('intensity').sum())
        .with_columns(
            (pl.col('intensity') / pl.col('intensity').max().over(['scan_idx', 'mass_idx']))
            .alias('abundance'),
            pl.lit('observed').alias('side'),
        )
    )


def _theoretical_envelopes(masses, averagine):
    """Averagine envelope rows for every (scan_idx, mass_idx, MonoMass)."""
    masses = masses.with_columns(
        (pl.col('MonoMass') // Averagine.BIN).cast(pl.Int64).alias('bin')
    )
    bins = masses.select('bin').unique().collect()['bin'].to_list()
    table = [
        (b, isotope, float(value))
        for b in bins
        for isotope, value in enumerate(averagine.envelope((b + 0.5) * Averagine.BIN))
    ]
    envelopes = pl.DataFrame(
        table, schema={'bin': pl.Int64, 'isotope': pl.Int64, 'abundance': pl.Float64},
        orient='row',
    ).lazy()
    return masses.join(envelopes, on='bin').with_columns(pl.lit('theoretical').alias('side'))


def _envelope_frame(masses_with_peaks, averagine):
    observed = _observed_envelopes(masses_with_peaks)
    # Only masses with signal peaks get a theoretical envelope
    with_signal = observed.select('scan_idx', 'mass_idx', 'MonoMass').unique()
    theoretical = _theoretical_envelopes(with_signal, averagine)
    columns = ['scan_idx', 'mass_idx', 'side', 'isotope', 'MonoMass', 'abundance']
    return (
        pl.concat([observed.select(columns), theoretical.select(columns)])
        .with_columns(
            (pl.col('MonoMass') + pl.col('isotope') * ISOTOPE_SPACING).alias('mass')
        )
        .select(list(ENVELOPE_SCHEMA))
    )


def isotope_envelopes(file_manager, dataset_id):
    """Observed and averagine envelope of every deconvolved mass, normalised
    to their highest isotope, keyed by (scan_idx, mass_idx).

    A fragmentation spectrum also gets rows with mass_idx -1: the envelope of
    its precursor mass in the precursor spectrum, as in the 3D S/N plot."""
    averagine = Averagine()
    masses = (
        insight_data.mass_table(file_manager, dataset_id)
        .select('scan_idx', 'mass_idx', 'MonoMass')
    )
    peaks = (
        insight_data.scan(file_manager, dataset_id, 'threedim_SN_plot')
        .select(pl.col('index').cast(pl.Int64).alias('scan_idx'), pl.col('SignalPeaks'))
        .with_columns(pl.int_ranges(0, pl.col('SignalPeaks').list.len()).alias('mass_idx'))
        .explode(['SignalPeaks', 'mass_idx'])
        .drop_nulls('mass_idx')
        .with_columns(pl.col('mass_idx').cast(pl.Int64))
    )
    with_peaks = masses.join(peaks, on=['scan_idx', 'mass_idx'], how='inner')
    envelopes = _envelope_frame(with_peaks, averagine).collect()

    precursors = _precursor_masses(file_manager, dataset_id)
    precursor_envelopes = (
        precursors.lazy()
        .join(
            envelopes.lazy().rename({'scan_idx': 'precursor_idx', 'mass_idx': 'precursor_mass_idx'}),
            on=['precursor_idx', 'precursor_mass_idx'],
        )
        .with_columns(pl.lit(-1, dtype=pl.Int64).alias('mass_idx'))
        .select(list(ENVELOPE_SCHEMA))
    )
    return (
        pl.concat([envelopes.lazy(), precursor_envelopes])
        .sort(['scan_idx', 'mass_idx', 'side', 'isotope'])
    )


def _precursor_masses(file_manager, dataset_id):
    """(scan_idx, precursor_idx, precursor_mass_idx) of every fragmentation
    spectrum whose precursor mass FLASHDeconv found in its precursor
    spectrum (within 0.01 Da, as the 3D S/N plot pairs them)."""
    scans = insight_data.scan_table(file_manager, dataset_id)
    threedim = insight_data.scan(file_manager, dataset_id, 'threedim_SN_plot')
    masses = (
        insight_data.mass_table(file_manager, dataset_id)
        .select(
            pl.col('scan_idx').alias('precursor_idx'),
            pl.col('mass_idx').alias('precursor_mass_idx'),
            'MonoMass',
        )
    )
    return (
        threedim.select(
            pl.col('index').cast(pl.Int64).alias('scan_idx'),
            pl.col('PrecursorScan').cast(pl.Int64),
        )
        .filter(pl.col('PrecursorScan') > 0)
        .join(scans.select(pl.col('index').alias('scan_idx'), 'PrecursorMass'), on='scan_idx')
        .join(
            scans.select(pl.col('Scan').alias('PrecursorScan'), pl.col('index').alias('precursor_idx'))
            .unique('PrecursorScan', keep='first'),
            on='PrecursorScan',
        )
        .join(masses, on='precursor_idx')
        .filter((pl.col('MonoMass') - pl.col('PrecursorMass')).abs() < 1e-2)
        .group_by('scan_idx')
        .agg(pl.col('precursor_idx').first(), pl.col('precursor_mass_idx').min())
        .sort('scan_idx')
        .collect()
    )


# --- Precursor XIC per charge -----------------------------------------------

XIC_SCHEMA = {
    'scan_idx': pl.Int64,
    'charge': pl.Int64,
    'series': pl.String,
    'rt': pl.Float64,
    'intensity': pl.Float64,
}


def precursor_xics(file_manager, dataset_id, tolerance_ppm=XIC_TOLERANCE_PPM):
    """For every fragmentation spectrum with a precursor mass: the intensity
    of that mass in every MS1 spectrum, per charge (summed signal peaks of the
    deconvolved masses within `tolerance_ppm`) and over all charges.

    One row per (scan_idx, charge, MS1 spectrum); charge 0 is the sum, named
    'All charges'. Rows come sorted by charge, then retention time, so the
    legend lists the charges in order."""
    scans = insight_data.scan_table(file_manager, dataset_id).collect()
    precursors = scans.filter(
        (pl.col('MSLevel') > 1) & (pl.col('PrecursorMass') > 0)
    ).select('index', 'PrecursorMass')
    ms1 = scans.filter(pl.col('MSLevel') == 1).select(
        pl.col('index').alias('ms1_idx'), pl.col('RT').cast(pl.Float64).alias('rt')
    )
    if precursors.is_empty() or ms1.is_empty():
        return pl.LazyFrame(schema=XIC_SCHEMA)

    ms1_masses = (
        insight_data.mass_table(file_manager, dataset_id)
        .select(pl.col('scan_idx').alias('ms1_idx'), 'mass_idx', 'MonoMass')
        .join(ms1.lazy(), on='ms1_idx')
        .collect()
        .sort('MonoMass')
    )
    if ms1_masses.is_empty():
        return pl.LazyFrame(schema=XIC_SCHEMA)

    # Masses of the MS1 spectra within the tolerance of each precursor mass
    sorted_masses = ms1_masses['MonoMass'].to_numpy()
    matches = []
    for scan_idx, precursor_mass in precursors.iter_rows():
        window = precursor_mass * tolerance_ppm * 1e-6
        lo = np.searchsorted(sorted_masses, precursor_mass - window, side='left')
        hi = np.searchsorted(sorted_masses, precursor_mass + window, side='right')
        for row in range(lo, hi):
            matches.append((scan_idx, row))
    if not matches:
        return pl.LazyFrame(schema=XIC_SCHEMA)
    matched = pl.DataFrame(matches, schema={'scan_idx': pl.Int64, 'row': pl.UInt32}, orient='row')
    matched = matched.join(
        ms1_masses.with_row_index('row'), on='row'
    ).select('scan_idx', 'ms1_idx', 'mass_idx', 'rt')

    peaks = (
        insight_data.scan(file_manager, dataset_id, 'threedim_SN_plot')
        .select(pl.col('index').cast(pl.Int64).alias('ms1_idx'), pl.col('SignalPeaks'))
        .join(matched.lazy().select('ms1_idx').unique(), on='ms1_idx', how='semi')
        .with_columns(pl.int_ranges(0, pl.col('SignalPeaks').list.len()).alias('mass_idx'))
        .explode(['SignalPeaks', 'mass_idx'])
        .drop_nulls('mass_idx')
        .with_columns(pl.col('mass_idx').cast(pl.Int64))
        .join(matched.lazy().select('ms1_idx', 'mass_idx').unique(), on=['ms1_idx', 'mass_idx'], how='semi')
        .explode('SignalPeaks')
        .drop_nulls('SignalPeaks')
        .select(
            'ms1_idx', 'mass_idx',
            pl.col('SignalPeaks').list.get(3, null_on_oob=True).cast(pl.Int64).alias('charge'),
            pl.col('SignalPeaks').list.get(2, null_on_oob=True).alias('intensity'),
        )
        .group_by(['ms1_idx', 'mass_idx', 'charge'])
        .agg(pl.col('intensity').sum())
    )
    per_charge = (
        matched.lazy()
        .join(peaks, on=['ms1_idx', 'mass_idx'])
        .group_by(['scan_idx', 'charge', 'ms1_idx', 'rt'])
        .agg(pl.col('intensity').sum())
    )
    total = (
        per_charge.group_by(['scan_idx', 'ms1_idx', 'rt'])
        .agg(pl.col('intensity').sum())
        .with_columns(pl.lit(0, dtype=pl.Int64).alias('charge'))
    )
    return (
        pl.concat([total.select(per_charge.collect_schema().names()), per_charge])
        .with_columns(
            pl.when(pl.col('charge') == 0)
            .then(pl.lit('All charges'))
            .otherwise(pl.format('{}+', pl.col('charge')))
            .alias('series'),
            pl.col('intensity').cast(pl.Float64),
        )
        .sort(['scan_idx', 'charge', 'rt'])
        .select(list(XIC_SCHEMA))
    )


# --- Characterization ----------------------------------------------------------

# Levels of the characterization badge in the protein table, best first
CHARACTERIZED = 'Characterized'
PARTIAL = 'Partial'
AMBIGUOUS = 'Ambiguous'
# Bonds explained by fragments, in % of the proteoform's bonds, at least for
# CHARACTERIZED / PARTIAL
CHARACTERIZED_MIN_CLEAVAGES = 25.0
PARTIAL_MIN_CLEAVAGES = 10.0

CHARACTERIZATION_SCHEMA = {
    'proteoform_index': pl.Int64,
    'Cleavages': pl.Float64,
    'LocalizedMods': pl.String,
    'Characterization': pl.String,
}


def explained_bonds(fragments):
    """Bonds explained by terminal fragments; bond b lies between database
    residues b and b + 1. A prefix fragment (a/b/c) ends at its residue, a
    suffix fragment (x/y/z) starts at it."""
    bonds = set()
    for fragment in fragments:
        ion = str(fragment['ion_type'])[:1]
        residue = fragment['residue']
        if residue is None:
            continue
        if ion in PREFIX_IONS:
            bonds.add(int(residue))
        elif ion in SUFFIX_IONS:
            bonds.add(int(residue) - 1)
    return bonds


def characterize(entry, fragments):
    """Characterization of one proteoform, ProSight's C-score traffic light in
    spirit but from FLASHTnT's own evidence:

    - Cleavages: % of the proteoform's bonds explained by terminal fragments.
    - LocalizedMods: modifications (mass shifts) FLASHTnT localized to one
      residue, of all of them ('k/n').
    - Characterization: CHARACTERIZED when both termini are determined, every
      modification is localized and at least CHARACTERIZED_MIN_CLEAVAGES % of
      the bonds are explained; AMBIGUOUS below PARTIAL_MIN_CLEAVAGES %;
      PARTIAL otherwise.
    """
    sequence = entry.get('sequence') or []
    start_reported = entry.get('proteoform_start')
    end_reported = entry.get('proteoform_end')
    start_reported = -1 if start_reported is None else int(start_reported)
    end_reported = -1 if end_reported is None else int(end_reported)
    start = max(start_reported, 0)
    end = end_reported if end_reported >= 0 else len(sequence) - 1
    total_bonds = max(end - start, 0)
    bonds = {b for b in explained_bonds(fragments) if start <= b < end}
    cleavages = 100.0 * len(bonds) / total_bonds if total_bonds else 0.0

    modifications = [m for m in (entry.get('modifications') or []) if m]
    localized = sum(1 for m in modifications if m.get('start') == m.get('end'))

    if cleavages < PARTIAL_MIN_CLEAVAGES:
        level = AMBIGUOUS
    elif (
        start_reported >= 0 and end_reported >= 0
        and localized == len(modifications)
        and cleavages >= CHARACTERIZED_MIN_CLEAVAGES
    ):
        level = CHARACTERIZED
    else:
        level = PARTIAL
    return {
        'Cleavages': round(cleavages, 1),
        'LocalizedMods': f'{localized}/{len(modifications)}' if modifications else '-',
        'Characterization': level,
    }


def characterization(file_manager, dataset_id):
    """Characterization (see characterize) of every proteoform."""
    errors = fragment_errors(file_manager, dataset_id).collect()
    by_proteoform = {
        key[0]: group.to_dicts()
        for key, group in errors.group_by(['proteoform_index'])
    }
    sequences = insight_data.proteoform_sequences(file_manager, dataset_id).select(
        'proteoform_index', 'sequence', 'start_position', 'end_position', 'modifications'
    ).collect()
    rows = []
    for entry in sequences.iter_rows(named=True):
        entry['proteoform_start'] = entry.pop('start_position')
        entry['proteoform_end'] = entry.pop('end_position')
        row = characterize(entry, by_proteoform.get(entry['proteoform_index'], []))
        row['proteoform_index'] = entry['proteoform_index']
        rows.append(row)
    return pl.DataFrame(rows, schema=CHARACTERIZATION_SCHEMA).lazy()
