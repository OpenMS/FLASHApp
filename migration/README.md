# Migration harness — FLASHApp page rebuild (Phase 3)

This directory tracks **Phase 3** of the migration: rebuild FLASHApp's viewer pages
on top of OpenMS-Insight, *via a reusable visualization template* in
`OpenMS/streamlit-template` — so the grid/layout layer is written once, frozen, and
imported unchanged by FLASHApp (no FLASHApp fork).

Phases 1 & 2 (OpenMS-Insight parity + simplification) are tracked in
`OpenMS-Insight/migration/`.

## Order of operations (freeze-then-propagate — no divergence)

1. **Schema prep** — `src/render/schema.py` post-processes FileManager caches into
   Insight-ready tidy parquet (stable IDs, exploded arrays, long-format density).
2. **Build & FREEZE the template** in `streamlit-template`:
   `src/view/grid.py` (`render_linked_grid` + `LayoutManager`),
   `content/visualization_template.py`, `src/common/common.py::show_linked_grid`,
   and the `src/workflow/FileManager.py` data-layer usage examples
   (store → `data_path` → Insight). Drive its review to ≥3 clean, then freeze.
3. **Rebuild FLASHApp** viewer pages from the frozen template: a builders factory
   (`comp_name -> BaseComponent`, `data_path=` parquet) + one `StateManager` per
   (tool, experiment); delete `src/render/{components,initialize,update,StateTracker}.py`.
4. **Iterate** three critics — template / original-FLASHApp parity / final — fixing at
   the **template level first**, then re-propagating, until ≥3 clean AND the
   **non-divergence gate** passes (FLASHApp grid == frozen template, by hash).

## Oracle (read-only)

`/home/user/FLASHApp/src/render/update.py` is the authoritative index→value selection
oracle; the old viewer pages `content/FLASH*/FLASH*Viewer.py` define the panels that
must all still render and cross-link.

## Files

- `units.yaml` — Phase 3 unit registry + gate definition + non-divergence file pairs.
- `run_review.py` — same convergence driver as OpenMS-Insight (`record`/`gate`/`report`).
- `nondivergence.py` — asserts FLASHApp's grid code is byte-identical (normalized) to
  the frozen template module.
- `review-log/phase-3.jsonl`, `REVIEW.md`, `specs/`.
