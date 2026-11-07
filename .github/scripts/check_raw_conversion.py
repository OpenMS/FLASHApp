"""Converts a Thermo .raw file the way FLASHApp does (CommandExecutor.convert_to_mzml).

    python check_raw_conversion.py <app dir> <file.raw> <workflow dir>

Prints the outcome, the number of spectra in the mzML and the workflow log;
exits 1 when the conversion failed. Temporary check for the OpenMS 3.6 switch.
"""
import os
import re
import sys
import time
from pathlib import Path

app, raw, wdir = sys.argv[1:4]
sys.path.insert(0, app)
os.chdir(app)
from src.workflow.Logger import Logger  # noqa: E402
from src.workflow.ParameterManager import ParameterManager  # noqa: E402
from src.workflow.CommandExecutor import CommandExecutor  # noqa: E402

w = Path(wdir)
for sub in ("pids", "logs"):
    Path(w, sub).mkdir(parents=True, exist_ok=True)
ex = CommandExecutor(w, Logger(w), ParameterManager(w))
start = time.time()
out = ex.convert_to_mzml([raw])
took = time.time() - start
log = Path(w, "logs", "all.log")
print(log.read_text(errors="replace")[-6000:] if log.exists() else "(no log)")
if not out:
    print(f"RESULT: failed after {took:.0f} s")
    sys.exit(1)
text = Path(out[0]).read_text(errors="replace")
n = len(re.findall(r"<spectrum\b", text))
levels = re.findall(r'name="ms level" value="(\d)"', text)
print(f"RESULT: {n} spectra (MS1 {levels.count('1')}, MS2 {levels.count('2')}) in {took:.0f} s -> {out[0]}")
