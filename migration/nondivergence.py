#!/usr/bin/env python3
"""
Non-divergence gate (Phase 3): assert FLASHApp's grid/layout code is the SAME code
as the frozen streamlit-template module — i.e. FLASHApp reuses the template verbatim
and never forks it.

Reads file pairs from migration/units.yaml -> meta.nondivergence_pairs:
    [[flashapp_path, template_path], ...]

For each pair both files are normalized (strip trailing whitespace, drop blank lines
and full-line comments) and compared by SHA-256:

  * both missing   -> PENDING (not yet built; passes, prints a note)
  * one missing    -> FAIL (a side exists but its counterpart does not)
  * present + equal -> OK
  * present + diff  -> FAIL (FLASHApp has diverged from the template)

Exit 0 iff no pair FAILs.
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

try:
    import yaml
except ImportError:  # pragma: no cover
    sys.exit("nondivergence.py requires pyyaml  (pip install pyyaml)")

ROOT = Path(__file__).resolve().parent
CONFIG = ROOT / "units.yaml"


def _normalized_hash(path: Path) -> str:
    lines = []
    for raw in path.read_text().splitlines():
        line = raw.rstrip()
        stripped = line.lstrip()
        if not stripped or stripped.startswith("#"):
            continue
        lines.append(line)
    return hashlib.sha256("\n".join(lines).encode()).hexdigest()


def main() -> int:
    cfg = yaml.safe_load(CONFIG.read_text()) or {}
    pairs = (cfg.get("meta") or {}).get("nondivergence_pairs") or []

    if not pairs:
        print("[nondivergence] no pairs configured yet (template grid not frozen) -> PENDING")
        return 0

    failed = False
    for pair in pairs:
        a, b = Path(pair[0]), Path(pair[1])
        ea, eb = a.exists(), b.exists()
        if not ea and not eb:
            print(f"[nondivergence] PENDING (both missing): {a.name}")
            continue
        if ea != eb:
            print(f"[nondivergence] FAIL (one side missing): {a} exists={ea} | {b} exists={eb}")
            failed = True
            continue
        ha, hb = _normalized_hash(a), _normalized_hash(b)
        if ha == hb:
            print(f"[nondivergence] OK: {a.name} == template")
        else:
            print(f"[nondivergence] FAIL (diverged): {a}\n                          != {b}")
            failed = True

    print(f"\n[nondivergence] {'RED' if failed else 'GREEN'}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
