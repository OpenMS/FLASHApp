"""The Viewer's modules must not load scipy.stats when imported.

Importing scipy.stats takes half a second or more (seconds on a cold disk),
and the first page load after every restart paid it before the view bar
appeared. The modules import it where they use it instead.
"""

import subprocess
import sys

import pytest


@pytest.mark.parametrize("module", [
    "src.render.compression",
    "src.parse.deconv",
    "src.parse.tnt",
])
def test_import_leaves_scipy_stats_unloaded(module):
    code = f"import sys, {module}; print('scipy.stats' in sys.modules)"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert out.stdout.strip().splitlines()[-1] == "False"
