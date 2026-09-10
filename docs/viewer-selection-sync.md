# Viewer selection sync

How the grid cells of the FLASHDeconv / FLASHTnT Viewer share a selection, and the invariants that keep the tables populated. Background: the September 2026 report "table contents only load when I select the first element" (root cause and reproduction in OpenMS/FLASHApp#100 and the linked session report).

## Moving parts

- Every grid cell is its own Streamlit component instance (one iframe of the Vue bundle). Cells share selection state only by round-tripping it through Python: `src/render/render.py::render_component` sends `selection_store` (the `StateTracker` state plus `counter` and `id`) with every render, and a cell that changes a selection calls `Streamlit.setComponentValue` with its whole store, which Python merges in `src/render/StateTracker.py`.
- One `StateTracker` exists per (tool, experiment selector). `render_grid` replaces it with a fresh one (new random `id`, empty state) whenever the selected experiment changes. The Vue store clears its local selection when it sees a new `id`.
- `render_component` hashes each cell's payload; the Vue store skips re-parsing a render whose hash equals the previous one (the heatmap payload is large).
- Tabulator tables select their first row when they are (re)built and emit that selection. This is what seeds `scanIndex` (Scan Table) or `proteinIndex` (Protein Table) after a load; the Mass Table, spectra and 3D plot are filtered by that selection on the Python side and are empty until it arrives.

## Invariants

1. **The payload hash includes the tracker id** (`src/render/util.py::payload_hash`). Two runs of the same input produce byte-identical tables; without the id the Scan Table skipped the render after an experiment switch, was never rebuilt, never seeded the default row, and every dependent cell stayed empty. The Vue store additionally refuses the hash skip when the tracker id changed, so either side alone is enough.
2. **`None` never claims a key.** The frontend sends `null` for every unset field. `StateTracker.updateState` adopts unknown keys only when their value is not `None`; otherwise the first cell to report would own every key and a later first-time value carrying the same counter (the Scan Table's default row) would be dropped as a stale conflict.
3. **A click never deselects.** Tabulator runs with `selectableRows: 'highlight'` and a `rowClick` handler that leaves exactly the clicked row selected, so clicking the highlighted row re-emits the selection instead of clearing it. Old Tabulator instances are destroyed before a table is rebuilt.
4. **Empty data clears the spectra.** The spectrum component replaces its plot with the "No Data Available" placeholder when its scan data goes away, instead of keeping the previous experiment's spectrum on screen. It no longer writes `massIndex` to the shared store on data changes; the cells that change the scan (Scan Table, Protein Table, heatmap click) reset the mass selection to the first mass.

## Known caveat

`StateTracker` ignores a conflicting update whose `counter` is older than the server's, which drops a genuine click sent while another cell's update is in flight (window of tens of milliseconds locally, longer with network latency). Such drops are logged at debug level (`StateTracker: ignored update ...`). Tracked in OpenMS/FLASHApp#100.

## Tests and reproduction

- `tests/test_selection_clear.py` pins invariants 1 and 2 (pure Python, no Streamlit needed).
- Vue: `src/stores/__tests__/streamlit-data.spec.ts` in openms-streamlit-vue-component pins the hash/tracker rule (`npx vitest run`).
- End to end: switch the Viewer dropdown between two runs of the same file; the Mass Table must fill without a click.
