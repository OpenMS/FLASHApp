#!/usr/bin/env python3
"""
Convergence-harness driver for the OpenMS-Insight -> FLASHApp parity migration.

This script does NOT spawn review agents (that orchestration happens in the Claude
Code session via the Agent tool). Its jobs are:

  record  - append a structured review result to the phase ledger
  gate    - run the machine gate (pytest / npm build / parity-diff) for a phase
  report  - read the ledger, compute per-round cleanliness + the consecutive
            clean-round counter, print the evidence, and exit 0 iff converged

A phase is CONVERGED when >= `meta.convergence` (default 3) consecutive rounds are
clean, where a round is clean iff every unit has a `clean` review record in that
round AND every gate step recorded for that round passed.

Ledger: one JSON object per line in migration/review-log/phase-<N>.jsonl
Fields: ts, phase, round, kind(review|gate|note), unit, status, findings[], msg
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

try:
    import yaml
except ImportError:  # pragma: no cover
    sys.exit("run_review.py requires pyyaml  (pip install pyyaml)")

ROOT = Path(__file__).resolve().parent            # migration/
REPO = ROOT.parent                                # repo root
CONFIG = ROOT / "units.yaml"
LOGDIR = ROOT / "review-log"


# --------------------------------------------------------------------------- io
def load_config() -> dict:
    with open(CONFIG) as fh:
        return yaml.safe_load(fh)


def phase_cfg(cfg: dict, phase) -> dict:
    phases = cfg.get("phases", {})
    pc = phases.get(str(phase)) or phases.get(int(phase))
    if pc is None:
        sys.exit(f"phase {phase} not defined in {CONFIG}")
    return pc


def unit_ids(pc: dict) -> list:
    return [u["id"] for u in pc.get("units", [])]


def ledger_file(phase) -> Path:
    LOGDIR.mkdir(parents=True, exist_ok=True)
    return LOGDIR / f"phase-{phase}.jsonl"


def append(phase, entry: dict) -> dict:
    entry = {"ts": time.strftime("%Y-%m-%dT%H:%M:%S"), "phase": int(phase), **entry}
    with open(ledger_file(phase), "a") as fh:
        fh.write(json.dumps(entry) + "\n")
    return entry


def read_ledger(phase) -> list:
    fp = ledger_file(phase)
    if not fp.exists():
        return []
    rows = []
    for line in fp.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            rows.append(json.loads(line))
        except json.JSONDecodeError:
            pass
    return rows


# --------------------------------------------------------------------- commands
def cmd_record(args) -> int:
    findings = []
    for f in args.finding or []:
        parts = f.split("|", 2)
        findings.append(
            {
                "id": parts[0],
                "severity": parts[1] if len(parts) > 1 else "info",
                "desc": parts[2] if len(parts) > 2 else "",
                "status": "open",
            }
        )
    entry = append(
        args.phase,
        {
            "round": args.round,
            "kind": "review",
            "unit": args.unit,
            "status": args.status,
            "findings": findings,
            "msg": args.msg or "",
        },
    )
    extra = f" ({len(findings)} finding(s))" if findings else ""
    print(f"recorded: round {entry['round']} unit {entry['unit']} -> {entry['status']}{extra}")
    return 0


def cmd_gate(args) -> int:
    cfg = load_config()
    pc = phase_cfg(cfg, args.phase)
    steps = pc.get("gate", [])
    if not steps:
        print(f"[gate] no gate steps configured for phase {args.phase}")
        return 0
    all_ok = True
    print(f"=== machine gate: phase {args.phase} round {args.round} ===")
    for step in steps:
        name, cmd = step["name"], step["cmd"]
        cwd = step.get("cwd", str(REPO))
        print(f"\n--- gate step: {name} ---\n$ {cmd}   (cwd={cwd})")
        proc = subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True)
        tail = "\n".join((proc.stdout + proc.stderr).splitlines()[-15:])
        ok = proc.returncode == 0
        all_ok = all_ok and ok
        if tail:
            print(tail)
        print(f"--> {name}: {'PASS' if ok else 'FAIL'} (rc={proc.returncode})")
        append(
            args.phase,
            {"round": args.round, "kind": "gate", "unit": name,
             "status": "pass" if ok else "fail", "msg": tail[-2000:]},
        )
    print(f"\n=== machine gate: {'GREEN' if all_ok else 'RED'} ===")
    return 0 if all_ok else 1


def cmd_report(args) -> int:
    cfg = load_config()
    pc = phase_cfg(cfg, args.phase)
    units = unit_ids(pc)
    conv = int(cfg.get("meta", {}).get("convergence", 3))
    rows = read_ledger(args.phase)

    rounds = sorted({r["round"] for r in rows if r.get("round") is not None})
    review_status, gate_records, fstate = {}, {}, {}
    for r in rows:
        rd = r.get("round")
        if r.get("kind") == "review":
            review_status[(rd, r.get("unit"))] = r.get("status")
        elif r.get("kind") == "gate":
            gate_records.setdefault(rd, []).append(r.get("status"))
        for f in r.get("findings") or []:
            fstate[f["id"]] = f.get("status", "open")

    def gate_ok(rd) -> bool:
        recs = gate_records.get(rd, [])
        return bool(recs) and all(s == "pass" for s in recs)

    def round_clean(rd) -> bool:
        units_clean = all(review_status.get((rd, u)) == "clean" for u in units) if units else True
        return units_clean and gate_ok(rd)

    print(f"\n================ REVIEW REPORT — phase {args.phase} ================")
    print(f"units: {len(units)} | rounds: {rounds or '—'} | convergence target: {conv}\n")
    if rounds:
        header = "unit".ljust(30) + "".join(f"R{rd}".rjust(5) for rd in rounds)
        print(header)
        print("-" * len(header))
        for u in units:
            line = u.ljust(30)
            for rd in rounds:
                line += {"clean": "✓", "finding": "✗"}.get(review_status.get((rd, u)), "·").rjust(5)
            print(line)
        print("GATE".ljust(30) + "".join(
            ("✓" if gate_ok(rd) else ("✗" if rd in gate_records else "·")).rjust(5) for rd in rounds))
        print("ROUND CLEAN".ljust(30) + "".join(
            ("✓" if round_clean(rd) else "✗").rjust(5) for rd in rounds))

    streak = 0
    for rd in rounds:
        streak = streak + 1 if round_clean(rd) else 0
    converged = streak >= conv
    open_ids = sorted(fid for fid, st in fstate.items() if st == "open")

    print(f"\nOPEN FINDINGS: {len(open_ids)}" + (": " + ", ".join(open_ids) if open_ids else ""))
    print(f"CONSECUTIVE CLEAN ROUNDS: {streak} / {conv}")
    print("STATUS: " + ("CONVERGED" if converged else "NOT CONVERGED"))

    if args.tail:
        tail = rows[-args.tail:]
        if tail:
            print(f"\n---- ledger tail (last {len(tail)}) ----")
            for r in tail:
                print(json.dumps(r))
    print("=" * 64)
    return 0 if converged else 1


# ------------------------------------------------------------------------- main
def main() -> None:
    p = argparse.ArgumentParser(description="migration convergence harness")
    sub = p.add_subparsers(dest="cmd", required=True)

    pr = sub.add_parser("record", help="append a review result")
    pr.add_argument("--phase", required=True)
    pr.add_argument("--round", type=int, required=True)
    pr.add_argument("--unit", required=True)
    pr.add_argument("--status", required=True, choices=["clean", "finding"])
    pr.add_argument("--finding", action="append", help="ID|severity|desc (repeatable)")
    pr.add_argument("--msg")
    pr.set_defaults(fn=cmd_record)

    pg = sub.add_parser("gate", help="run the machine gate for a phase")
    pg.add_argument("--phase", required=True)
    pg.add_argument("--round", type=int, required=True)
    pg.set_defaults(fn=cmd_gate)

    rp = sub.add_parser("report", help="print rollup + convergence status")
    rp.add_argument("--phase", required=True)
    rp.add_argument("--tail", type=int, default=12)
    rp.set_defaults(fn=cmd_report)

    args = p.parse_args()
    sys.exit(args.fn(args) or 0)


if __name__ == "__main__":
    main()
