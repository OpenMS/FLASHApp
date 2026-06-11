# Migration review rollup — FLASHApp page rebuild (Phase 3)

> The live rollup matrix + the `CONSECUTIVE CLEAN ROUNDS: k / 3` counter are printed
> by `python migration/run_review.py report --phase 3`.

## Status

| Phase | Description | Converged? |
|------:|-------------|:----------:|
| 3 | Rebuild FLASHApp viewer pages from the frozen template (no divergence) | ⏳ not started |

Convergence target: **≥3 consecutive clean rounds** (every unit clean + machine gate
green, including the non-divergence check). Prereq: Phases 1 & 2 converged in
`OpenMS-Insight/migration/`.

## Units (see `units.yaml`)

- **Template (built & frozen first):** `template:grid`, `template:page`,
  `template:common`, `template:filemanager`.
- **FLASHApp rebuild (from frozen template):** `flashapp:schema`, `flashapp:builders`,
  `flashapp:deconv-viewer`, `flashapp:tnt-viewer`, `flashapp:quant-viewer`,
  `flashapp:nondivergence`.

Critics per unit: **template / original-FLASHApp parity / final**. Fixes land at the
**template level first**, then re-propagate to FLASHApp.
