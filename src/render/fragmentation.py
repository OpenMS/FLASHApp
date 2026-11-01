"""Fragmentation method of a run, read from its mzML, for the Sequence View.

The Sequence View (openms-insight ProteinSequenceView) picks the ion series to
match from a fragmentation method, as ProSight Lite does: CID/HCD -> b, y;
ETD/ECD -> c, z•; EThcD -> b, c, y, z•; UVPD -> nine series. This module finds
the method most MS2 spectra of a run were acquired with. Streamlit-free.
"""
from collections import Counter
from functools import lru_cache
from pathlib import Path

import pyopenms as oms

_AM = oms.Precursor.ActivationMethod

# Collisional activation: b/y ions
COLLISIONAL = {
    _AM.CID, _AM.HCID, _AM.LCID, _AM.HCD, _AM.TRAP, _AM.PQD,
    _AM.BIRD, _AM.SORI, _AM.IMD, _AM.SID, _AM.PSD, _AM.INSOURCE, _AM.LIFT,
}
# Electron-based activation: c/z• ions
ELECTRON = {_AM.ETD, _AM.ECD}
# Combined electron transfer and collisional activation
COMBINED = {_AM.EThcD, _AM.ETciD}


def method_of(activation_methods):
    """Sequence View fragmentation method ('CID/HCD', 'ETD/ECD', 'EThcD') for
    the activation methods of one precursor, or None when they do not say
    (none given, or photodissociation, which can be IRMPD or UVPD)."""
    methods = set(activation_methods)
    if methods & COMBINED or (methods & ELECTRON and methods & COLLISIONAL):
        return 'EThcD'
    if methods & ELECTRON:
        return 'ETD/ECD'
    if methods & COLLISIONAL:
        return 'CID/HCD'
    return None


def run_method(experiment):
    """Most common fragmentation method of the MS2+ spectra of an MSExperiment,
    or None when no spectrum states a known one."""
    counts = Counter()
    for spectrum in experiment:
        if spectrum.getMSLevel() < 2:
            continue
        for precursor in spectrum.getPrecursors():
            method = method_of(precursor.getActivationMethods())
            if method is not None:
                counts[method] += 1
    if not counts:
        return None
    # Ties: the order above (collisional first) as the stable fallback
    order = ['CID/HCD', 'ETD/ECD', 'EThcD']
    return max(counts, key=lambda m: (counts[m], -order.index(m)))


@lru_cache(maxsize=64)
def _mzml_method(path, size, mtime_ns):
    options = oms.PeakFileOptions()
    options.setFillData(False)
    options.setMSLevels([2, 3, 4])
    mzml = oms.MzMLFile()
    mzml.setOptions(options)
    experiment = oms.MSExperiment()
    mzml.load(path, experiment)
    return run_method(experiment)


def mzml_method(path):
    """Fragmentation method of an mzML file (spectra metadata only), cached
    per file version. None when the file is missing or unreadable."""
    path = Path(path)
    try:
        stat = path.stat()
        return _mzml_method(str(path), stat.st_size, stat.st_mtime_ns)
    except Exception:  # pylint: disable=broad-except
        return None
