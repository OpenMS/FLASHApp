"""Compare two mzML files spectrum by spectrum (temporary check)."""
import sys
import numpy as np
import pyopenms as oms

def load(path):
    e = oms.MSExperiment()
    oms.MzMLFile().load(path, e)
    return e

a, b = load(sys.argv[1]), load(sys.argv[2])
na, nb = a.getNrSpectra(), b.getNrSpectra()
print(f"spectra: {na} vs {nb}")
diff = {"level": 0, "rt": 0, "npeaks": 0, "mz": 0, "int": 0, "prec_mz": 0, "prec_z": 0, "polarity": 0, "native_id": 0}
max_mz = max_int_rel = max_rt = max_pmz = 0.0
examples = []
for i in range(min(na, nb)):
    s, t = a[i], b[i]
    if s.getMSLevel() != t.getMSLevel(): diff["level"] += 1
    if s.getNativeID() != t.getNativeID(): diff["native_id"] += 1
    if s.getInstrumentSettings().getPolarity() != t.getInstrumentSettings().getPolarity(): diff["polarity"] += 1
    d = abs(s.getRT() - t.getRT()); max_rt = max(max_rt, d)
    if d > 1e-4: diff["rt"] += 1
    ms, is_ = s.get_peaks(); mt, it = t.get_peaks()
    if len(ms) != len(mt):
        diff["npeaks"] += 1
        if len(examples) < 5: examples.append(f"#{i} {s.getNativeID()}: {len(ms)} vs {len(mt)} peaks")
    elif len(ms):
        dm = float(np.max(np.abs(ms - mt))); max_mz = max(max_mz, dm)
        if dm > 1e-6: diff["mz"] += 1
        rel = float(np.max(np.abs(is_ - it) / np.maximum(np.abs(is_), 1e-12))); max_int_rel = max(max_int_rel, rel)
        if rel > 1e-5: diff["int"] += 1
    ps, pt = s.getPrecursors(), t.getPrecursors()
    if len(ps) != len(pt): diff["prec_mz"] += 1
    for p, q in zip(ps, pt):
        dp = abs(p.getMZ() - q.getMZ()); max_pmz = max(max_pmz, dp)
        if dp > 1e-6: diff["prec_mz"] += 1
        if p.getCharge() != q.getCharge(): diff["prec_z"] += 1
lv = lambda e: [sum(1 for s in e if s.getMSLevel() == l) for l in (1, 2)]
print("MS1/MS2:", lv(a), "vs", lv(b))
print("spectra differing:", diff)
print(f"max |dRT| {max_rt:.2e} s, max |dm/z| {max_mz:.2e}, max rel. dIntensity {max_int_rel:.2e}, max |dprecursor m/z| {max_pmz:.2e}")
print("first native IDs:", a[0].getNativeID(), "|", b[0].getNativeID())
for x in examples: print(" ", x)
print("chromatograms:", a.getNrChromatograms(), "vs", b.getNrChromatograms())
s, t = a[0], b[0]
print("spectrum 0 meta A:", s.getKeys() if hasattr(s, "getKeys") else "", "type", s.getType())
print("spectrum 0 meta B:", t.getKeys() if hasattr(t, "getKeys") else "", "type", t.getType())
