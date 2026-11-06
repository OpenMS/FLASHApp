"""Theoretical isotope envelope of a deconvolved mass, drawn on its raw peaks.

For the selected mass of a spectrum, every charge state that has signal peaks
gets the averagine isotope envelope of the mass at that charge: one peak per
isotope at m/z (M + i * 1.00335) / z + proton mass, scaled to the observed
peaks of that charge by least squares. The annotated spectrum draws them as
circles over the raw peaks, so a wrong charge, a monoisotopic mass off by one
isotope or a missing envelope shows at a glance.

Signal peaks are FLASHDeconv's (3D S/N plot data): [peak index, m/z,
intensity, charge] per peak. pyOpenMS supplies the averagine envelopes
(validation_data.Averagine). Like insight_data, the module is free of
Streamlit and openms-insight.
"""

from collections import defaultdict

import numpy as np
import polars as pl

from src.render import insight_data
from src.render.validation_data import ISOTOPE_SPACING, PROTON_MASS, Averagine

# Theoretical isotopes below this fraction of the highest one are not drawn
MIN_DRAWN_ABUNDANCE = 0.01

_AVERAGINE = None


def _averagine():
    global _AVERAGINE
    if _AVERAGINE is None:
        _AVERAGINE = Averagine()
    return _AVERAGINE


def theoretical_peaks(mono_mass, signal_peaks, abundances, min_abundance=MIN_DRAWN_ABUNDANCE):
    """Theoretical peaks of `mono_mass` for each charge in `signal_peaks`:
    dicts with mz, intensity, charge and isotope, by charge then isotope.

    abundances: relative isotope abundances of the mass (isotope 0, 1, ...).
    """
    abundances = np.asarray(abundances, dtype=float)
    if mono_mass is None or not mono_mass > 0 or abundances.size == 0 or abundances.max() <= 0:
        return []
    abundances = abundances / abundances.max()
    # Observed intensity per (charge, isotope)
    observed = defaultdict(lambda: defaultdict(float))
    for peak in signal_peaks or []:
        if peak is None or len(peak) < 4:
            continue
        mz, intensity, charge = peak[1], peak[2], int(peak[3])
        if charge <= 0 or mz is None or intensity is None:
            continue
        isotope = int(round(((mz - PROTON_MASS) * charge - mono_mass) / ISOTOPE_SPACING))
        if 0 <= isotope < abundances.size:
            observed[charge][isotope] += float(intensity)
    peaks = []
    for charge in sorted(observed):
        isotopes = np.array(sorted(observed[charge]))
        measured = np.array([observed[charge][i] for i in isotopes])
        expected = abundances[isotopes]
        denominator = float(np.dot(expected, expected))
        if denominator <= 0:
            continue
        scale = float(np.dot(measured, expected)) / denominator
        if scale <= 0:
            continue
        for isotope in np.nonzero(abundances >= min_abundance)[0]:
            peaks.append({
                'mz': (mono_mass + isotope * ISOTOPE_SPACING) / charge + PROTON_MASS,
                'intensity': scale * abundances[isotope],
                'charge': charge,
                'isotope': int(isotope),
            })
    return peaks


def selected_mass_envelope(file_manager, dataset_id, scan_idx, mass_idx):
    """theoretical_peaks() of mass `mass_idx` of spectrum `scan_idx`, or []."""
    if scan_idx is None or mass_idx is None or mass_idx < 0:
        return []
    for name in ('threedim_SN_plot', 'mass_table'):
        if insight_data.result_path(file_manager, dataset_id, name) is None:
            return []
    masses = (
        insight_data.scan(file_manager, dataset_id, 'mass_table')
        .filter(pl.col('index').cast(pl.Int64) == int(scan_idx))
        .select('MonoMass')
        .collect()
    )
    signal = (
        insight_data.scan(file_manager, dataset_id, 'threedim_SN_plot')
        .filter(pl.col('index').cast(pl.Int64) == int(scan_idx))
        .select('SignalPeaks')
        .collect()
    )
    if masses.is_empty() or signal.is_empty():
        return []
    mono_masses = masses['MonoMass'].to_list()[0]
    peaks_per_mass = signal['SignalPeaks'].to_list()[0]
    if mono_masses is None or peaks_per_mass is None:
        return []
    if mass_idx >= len(mono_masses) or mass_idx >= len(peaks_per_mass):
        return []
    mono_mass = mono_masses[mass_idx]
    return theoretical_peaks(mono_mass, peaks_per_mass[mass_idx], _averagine().envelope(mono_mass))
