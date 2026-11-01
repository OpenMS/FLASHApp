"""Fragmentation method of each scan, read from its mzML, for the Sequence View.

The Sequence View (openms-insight ProteinSequenceView) picks the ion series to
match from a fragmentation method, as ProSight Lite does: CID/HCD -> b, y;
ETD/ECD -> c, z•; EThcD -> b, c, y, z•; UVPD -> nine series. This module reads
the method each spectrum was acquired with from the precursor activation
methods FLASHDeconv carries over from the uploaded mzML, so a run that mixes
methods (e.g. alternating HCD and ETD scans) gets the right series per scan.
Streamlit-free.
"""
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


def spectrum_method(spectrum):
    """Fragmentation method of one spectrum (None for MS1 or when unknown)."""
    if spectrum.getMSLevel() < 2:
        return None
    for precursor in spectrum.getPrecursors():
        method = method_of(precursor.getActivationMethods())
        if method is not None:
            return method
    return None


@lru_cache(maxsize=64)
def _mzml_scan_methods(path, size, mtime_ns):
    options = oms.PeakFileOptions()
    options.setFillData(False)
    mzml = oms.MzMLFile()
    mzml.setOptions(options)
    experiment = oms.MSExperiment()
    mzml.load(path, experiment)
    return tuple(spectrum_method(spectrum) for spectrum in experiment)


def mzml_scan_methods(path):
    """Fragmentation method of every spectrum of an mzML file, by spectrum
    index (None for MS1 and unknown methods). Reads spectra metadata only,
    cached per file version. An empty tuple when the file is missing or
    unreadable."""
    path = Path(path)
    try:
        stat = path.stat()
        return _mzml_scan_methods(str(path), stat.st_size, stat.st_mtime_ns)
    except Exception:  # pylint: disable=broad-except
        return ()

