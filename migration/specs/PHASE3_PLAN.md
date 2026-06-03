# Phase 3 Plan — Rebuild FLASHApp viewers on OpenMS-Insight via a frozen `streamlit-template` grid

**Goal.** Re-implement the three FLASHApp visualization pages (FLASHDeconv, FLASHTnT,
FLASHQuant) on top of the parity-complete `openms-insight` package, through a
**single reusable grid module that lives in `OpenMS/streamlit-template` and is imported by
FLASHApp byte-for-byte unchanged**. The template is built and *frozen* first; FLASHApp then
rebuilds against the frozen module so `migration/nondivergence.py` is GREEN.

**Scope discipline.** This is a planning doc. The build order, exact signatures, per-component
tidy-parquet schemas, and the non-divergence mechanism below are the contract. Oracle behavior
to preserve = the current FLASHApp render layer (`src/render/*`) + the two `FLASH*LayoutManager`
pages + the three viewer pages. Everything the new design *deletes* is listed in §5.5.

---

## 0. Background: what the oracle does today (so we preserve it)

The current grid is a **bespoke Vue mega-component** (`js-component/`, declared in
`src/render/components.py::get_component_function`) that receives the *entire* per-panel
dataset plus a `selection_store`, and does selection/filtering Python-side every rerun:

- `render.py::render_grid(selected_data, layout_info_per_exp, file_manager, tool, identifier, grid_key)`
  iterates `layout_info_per_exp` (a list of rows, each row a list of `comp_name` strings, ≤3 cols),
  `st.columns(len(row))` per row, and for each cell:
  1. `initialize.py::initialize_data(comp_name, selected_data, file_manager, tool)` loads the
     cache(s) for that panel into a `(data_to_send, components, additional_data)` triple, keyed
     in `st.session_state['plot_data'][tool][identifier][comp_name]`.
  2. `render.py::render_component(...)` runs `update.py::update_data` then `filter_data`, hashes,
     and calls the single Vue component, then reconciles state via `StateTracker`.
- **Selection is index-based** (`update.py`): `selection_store['scanIndex']` slices `per_scan_data`
  by `.iloc[scanIndex:scanIndex+1]`; `massIndex` indexes into `SignalPeaks[massIndex]`;
  `proteinIndex` keys a `proteoform_scan_map`; heatmap zoom is `xRange/yRange`. This is the exact
  oracle we must reproduce with Insight's **value-based** `filters`/`interactivity`.
- `StateTracker.py` is a per-(tool,identifier) counter+id reconciler — **the local twin of
  `openms_insight.StateManager`** (compare `StateManager.update_from_vue`). It is replaced 1:1.
- Layout managers (`FLASHDeconvLayoutManager.py`, `FLASHTnTLayoutManager.py`) are ~330 lines each,
  **near-identical** apart from `COMPONENT_OPTIONS`/`COMPONENT_NAMES` and session-state key names.
  They edit a 3-level nested list (`[exp][row][col] = option-label`), enforce ≤3 columns, validate
  `"(... needed)"` dependencies, persist `{'layout': trimmed, 'side_by_side': bool}` to the
  FileManager under `('layout','layout')` (deconv) / `('flashtnt_layout','layout')` (tnt), and
  support JSON download / upload. **This duplication is the distillation target for `LayoutManager`.**

Data layer: `src/workflow/FileManager.py` is a SQLite-indexed results store keyed by
`(dataset_id, name_tag)`, writing DataFrames as `.pq` and everything else as `.pkl.gz`, with
`get_results(dataset_id, name_tags, use_pyarrow=, use_polars=)`. For `.pq` columns it returns a
**pandas DF** (default), a **polars LazyFrame** (`use_polars`), or a **pyarrow Dataset**
(`use_pyarrow`). It does *not* expose a "give me the parquet path" mode — Insight wants
`data_path=`. We add exactly that (§2).

OpenMS-Insight public surface we build on (from `openms_insight/__init__.py`, `core/base.py`,
`core/state.py`, README): 7 components subclassing `BaseComponent`, each
`Comp(cache_id, data=/data_path=, filters=, filter_defaults=, interactivity=, cache_path=, **cfg)`
and render-time `comp(key=, state_manager=, height=, **render_switches)`. `StateManager(session_key=)`
routes selections by identifier. Crucially: `data_path=` triggers **subprocess preprocessing** and
disk cache keyed by `cache_id`+config-hash; presentation (titles/labels/colors/thresholds) is
render-time and never rebuilds the cache.

---

## 1. `streamlit-template/src/view/grid.py` (NEW) — the single source of truth

A new package `streamlit-template/src/view/` (add `src/view/__init__.py`). `grid.py` is
**tool-agnostic**: it knows nothing about FLASHDeconv/TnT/Quant, scans, masses, or proteins.
It distills `render.py::render_grid` + both `FLASH*LayoutManager` classes into two public objects:
`render_linked_grid(...)` and `LayoutManager`.

### 1.1 `render_linked_grid` — exact signature

```python
# streamlit-template/src/view/grid.py
from typing import Callable, Dict, List, Optional, Sequence
import streamlit as st
from openms_insight import StateManager, BaseComponent

# A layout is the trimmed nested list the LayoutManager persists:
#   List[row], row = List[comp_name:str], <=3 entries per row.   (one experiment)
Layout = List[List[str]]
# `builders` maps a comp_name -> a zero-arg factory returning a *constructed* BaseComponent.
# Zero-arg so the grid can lazily build only the panels a given layout references, and so the
# factory closes over (dataset_id, file_manager, cache_path) on the FLASHApp side (see §5.2).
BuilderMap = Dict[str, Callable[[], BaseComponent]]


def render_linked_grid(
    layout: Layout,
    builders: BuilderMap,
    state_key: str,
    *,
    grid_key: str = "linked_grid",
    height: Optional[int] = None,
    column_heights: Optional[Dict[str, int]] = None,
    on_missing: str = "warn",          # "warn" | "error" | "skip"
) -> StateManager:
    """Render one experiment's linked grid.

    For each row in `layout`, open `st.columns(len(row))` (clamped to <=3, mirroring the
    oracle's hard cap) and, in each column, call `builders[comp_name]()` to construct the
    Insight component, then render it with a SHARED `StateManager(session_key=state_key)` and a
    per-cell Streamlit key `f"{grid_key}_{r}_{c}"`. The shared StateManager is what cross-links
    every panel in the grid: clicks (`interactivity`) write selections, other panels read them
    (`filters`). Returns the StateManager so callers can introspect/seed selections.

    Args
      layout       : trimmed nested list (rows of comp_names) for ONE experiment.
      builders     : comp_name -> () -> BaseComponent  (factory; see BuilderMap).
      state_key    : StateManager session_key. MUST be unique per (tool, experiment) so two
                     experiments shown together do not share selections (see §5.3). This is the
                     direct replacement for the oracle's (tool, identifier) StateTracker scoping.
      grid_key     : prefix for per-cell component keys (replaces oracle `grid_key`).
      height       : default px height passed to every comp's __call__ (None -> Insight default).
      column_heights: optional comp_name -> height override (e.g. heatmaps taller).
      on_missing   : behavior when a comp_name has no builder ("warn" st.warning + skip).
    """
```

**Render loop (the distilled `render_grid` inner body), reference implementation:**

```python
    sm = StateManager(session_key=state_key)
    n_rows = len(layout)
    for r, row in enumerate(layout):
        cols = st.columns(min(len(row), 3))           # <=3 columns, oracle invariant
        for c, comp_name in enumerate(row[:3]):
            factory = builders.get(comp_name)
            if factory is None:
                if on_missing == "error":
                    raise KeyError(f"No builder registered for component '{comp_name}'")
                if on_missing == "warn":
                    cols[c].warning(f"Unknown component: {comp_name}")
                continue
            h = (column_heights or {}).get(comp_name, height)
            with cols[c]:
                factory()(key=f"{grid_key}_{r}_{c}", state_manager=sm, height=h)
    return sm
```

Design notes that preserve oracle behavior:
- **State scoping.** The oracle nests `st.session_state['state_tracker'][tool][identifier]`. We
  achieve the same isolation purely through `StateManager(session_key=state_key)` —
  StateManager stores under `st.session_state[session_key]`, so distinct `state_key`s are fully
  independent (matches `render_grid`'s per-identifier tracker and its "dataset changed -> reset"
  behavior, which now falls out of cache_id+state_key changing per dataset, see §5.3).
- **Dataset-change reset.** The oracle wipes `plot_data`/`state_tracker` when `selected_data`
  changes. Equivalent here: the FLASHApp builders bake `dataset_id` into both `cache_id` and
  `state_key` (§5.2/§5.3), so selecting another experiment yields a fresh StateManager + fresh
  component caches automatically — no manual reset code in the template.
- **No data plumbing in the template.** Unlike `render_component`, the grid never touches data,
  hashing, or `update/filter`. All of that moved *into* each Insight component's `_preprocess`
  + `filters`/`interactivity`. The grid is pure layout + a shared StateManager. This is what
  makes it tool-agnostic and safe to freeze.
- **`@st.fragment`.** Do **not** decorate `render_linked_grid` itself (it opens `st.columns` for
  the caller's container). Individual Insight components already manage their own rerun via
  StateManager. (Side-by-side wrapping uses fragments at the page level — see §3/§5.3.)

### 1.2 `LayoutManager` — exact API (distillation of both `FLASH*LayoutManager`)

A class that owns the layout-editor UI + persistence, parameterized by the things that differ
between the two FLASH managers (component vocabulary, storage keys, session namespace).

```python
class LayoutManager:
    def __init__(
        self,
        component_options: List[str],      # human labels, e.g. "Scan table"
        component_names:   List[str],      # parallel internal names, e.g. "scan_table"
        *,
        store,                             # object with get/set/exists/remove (see Store proto)
        layout_id: str = "layout",         # FileManager dataset_id for the saved layout
        layout_tag: str = "layout",        # FileManager name_tag for the saved layout
        max_columns: int = 3,
        max_experiments: int = 5,
        session_prefix: str = "lm",        # namespaces all st.session_state keys
        download_name: str = "layout_settings.json",
        title: str = "Layout Manager",
    ): ...

    # --- persistence (replaces set_layout/get_layout in both managers) ---
    def get_layout(self) -> Optional[tuple[list, bool]]:
        """Return (layout_per_experiment, side_by_side) or None if unset.
        layout_per_experiment: List[experiment], experiment = List[row], row = List[comp_name]."""
    def set_layout(self, layout: list, side_by_side: bool = False) -> None: ...

    # --- label<->name transforms (oracle getTrimmed/getExpanded) ---
    def trim(self, expanded: list) -> list:     # labels -> internal names, drop empties
    def expand(self, trimmed: list) -> list:    # internal names -> labels

    # --- validation (oracle validateSubmittedLayout: non-empty + "(... needed)" deps) ---
    def validate(self, layout: Optional[list] = None) -> str:   # '' if OK else message

    # --- the whole editor page (renders edit/saved modes, buttons, upload/download, tips) ---
    def render(self) -> None:
        """Draw the full Layout Manager page exactly like the oracle: experiment count
        selector, per-experiment expanders with add-column(+)/add-row(+)/delete(x) controls,
        the <=3-column cap, side-by-side checkbox (shown when #experiments==2), Save/Edit/
        Reset, JSON download (disabled while invalid) + JSON upload, success/error toasts, tips."""

    # --- extension hook for FLASHDeconv's dynamic "Sequence view" option ---
    def add_options(self, options: List[str], names: List[str]) -> None:
        """Append (label, name) pairs at runtime (oracle setSequenceView: adds Sequence/Internal
        options once an input sequence exists)."""
```

`Store` protocol (so the template does not import FLASHApp's FileManager — it accepts any object
implementing the 4 calls; FLASHApp passes its FileManager, the template demo passes the template
FileManager from §2):

```python
class Store(Protocol):
    def get_results(self, dataset_id: str, name_tags: list) -> dict: ...
    def store_data(self, dataset_id: str, name_tag: str, data) -> None: ...
    def result_exists(self, dataset_id: str, name_tag: str) -> bool: ...
    def remove_results(self, dataset_id: str) -> None: ...
```

**Why a class, not free functions:** the two oracle managers are 95% duplicated free-function
modules whose only real differences are the vocab lists and the `*_tagger` session-key suffix.
Folding them into one class parameterized by `component_options/names`, `layout_id/tag`, and
`session_prefix` removes the duplication while keeping the exact UI/JSON-format/validation
behavior. The deconv manager becomes `LayoutManager(DECONV_OPTIONS, DECONV_NAMES, store=fm,
layout_id="layout", session_prefix="deconv"); lm.add_options(...); lm.render()`; the tnt manager
becomes the same with TNT vocab, `layout_id="flashtnt_layout"`, `session_prefix="tnt"`.

**Behavioral invariants to preserve (verbatim from the oracle):**
- ≤3 columns per row (column "+"/delete "x"/row "+" controls).
- `"<Component> (X needed)"` dependency validation (`X` must also be present in the same exp).
- Saved JSON is the **trimmed internal-name** nested list (so old saved layouts keep loading).
- `side_by_side` only offered when exactly 2 experiments; persisted alongside the layout.
- "If nothing is set, default layout is used in the Viewer" (Viewer supplies `DEFAULT_LAYOUT`).

> **Freeze point.** Once §1 + §2 + §3 land and tests pass, this file is *frozen*: FLASHApp must
> consume it unchanged (§6). Register the pair in `units.yaml`.

---

## 2. `streamlit-template/src/workflow/FileManager.py` — results-store data layer

The template's current FileManager (180 lines) only does path munging — it has **no caching/store
API at all**. Port the richer FLASHApp FileManager (SQLite-indexed `(dataset_id, name_tag)` store
with parquet/pickle) into the template, and add the one method Insight needs: **return the parquet
path** so it can be handed to `data_path=`.

### 2.1 What to port (verbatim from FLASHApp `src/workflow/FileManager.py`)

Bring over, unchanged in behavior: `__init__(workflow_dir, cache_path=None)`, the SQLite
`_connect_to_sql`/`__getstate__`/`__setstate__`, `_add_column/_add_entry`, `store_data` (+ the
`_store_data` parquet/pickle split and `row_group_size`), `parquet_sink` contextmanager,
`store_file`, `get_results_list`, `get_results(..., use_pyarrow=, use_polars=)`, `result_exists`,
`remove_results`, `clear_cache`, `get_display_name`, `rename_dataset`. Keep the existing
`get_files`/`_set_type`/`_set_dir` path helpers (the template's current contract) so existing
template pages still work.

### 2.2 NEW method — parquet path for Insight `data_path=`

```python
def get_results(self, dataset_id, name_tags, partial=False,
                use_pyarrow=False, use_polars=False, as_path=False):
    """... existing behavior ... PLUS:
    as_path=True  -> for parquet (.pq) columns, return the str path to the parquet file
                     (NOT a loaded frame), so it can be passed straight to an Insight
                     component's data_path=. Pickle (.pkl.gz) columns still load + return
                     the object (there is no path contract for non-tabular data)."""
```

Implementation: in the data-column branch, when `as_path=True` and `file_path.suffix == '.pq'`,
set `results[c] = str(file_path)` instead of reading it. (Mutually exclusive with
`use_pyarrow`/`use_polars`; if more than one is set, precedence `as_path > use_pyarrow >
use_polars > pandas`, documented in the docstring.) This is the `get_results(..., use_pyarrow=True)`-
style API the prompt calls for, generalized to "give me the path".

Convenience wrapper (sugar used pervasively by the FLASHApp builders in §5.2):

```python
def result_path(self, dataset_id: str, name_tag: str) -> str:
    """Return the on-disk parquet path for a single (dataset_id, name_tag), or raise KeyError.
    Equivalent to get_results(dataset_id, [name_tag], as_path=True)[name_tag]."""
```

### 2.3 Usage example (store -> data_path -> Insight) — goes in the docstring + §3 demo

```python
from src.workflow.FileManager import FileManager
from openms_insight import Heatmap, StateManager
import polars as pl

fm = FileManager(workspace_dir, cache_path=workspace_dir / "cache")

# 1) store a (lazy) frame -> parquet, indexed by (dataset_id, name_tag)
fm.store_data("demo", "peaks", pl.scan_parquet("raw_peaks.parquet"))

# 2) hand the parquet PATH to an Insight component (subprocess preprocessing + disk cache)
sm = StateManager(session_key="demo_grid")
Heatmap(
    cache_id="demo_peaks_heatmap",
    data_path=fm.result_path("demo", "peaks"),     # <- the new path API
    x_column="rt", y_column="mass", intensity_column="intensity",
    cache_path=str(fm.cache_path / "insight"),     # keep Insight caches under the workspace
)(state_manager=sm)
```

> Note: store the layout dict (`{'layout': ..., 'side_by_side': ...}`) via `store_data` exactly as
> the oracle does — it's a plain dict, so it round-trips through the `.pkl.gz` branch unchanged.
> The `LayoutManager.Store` protocol (§1.2) is satisfied by this FileManager directly.

---

## 3. `streamlit-template/content/visualization_template.py` (NEW) — demo page

A self-contained demo registered in `app.py` that exercises the full stack on **small example
parquet** under `example-data/insight/`, so the template proves the grid + LayoutManager +
side-by-side + `Table<->LinePlot<->Heatmap<->SequenceView` linking end-to-end (and is the
`template:page` oracle for the FLASHApp viewers).

### 3.1 Example data to generate (committed under `example-data/insight/`)

Tiny, hand-built parquet (a one-off generator script `example-data/insight/_make_example.py`,
run once; commit the `.parquet`). Schemas chosen to match the Insight components' tidy contracts
(§4):
- `spectra.parquet` — master table: `scan_id:int, rt:float, ms_level:int, precursor_mz:float, n_peaks:int` (~20 rows).
- `peaks.parquet` — per-peak long format: `scan_id:int, peak_id:int, mass:float, intensity:float, is_annotated:int, ion_label:str` (~400 rows; `peak_id` globally unique).
- `heat.parquet` — peak map: `scan_id:int, rt:float, mass:float, intensity:float, peak_id:int` (a few thousand rows).
- `sequences.parquet` — `scan_id:int, sequence:str, precursor_charge:int` (one seq per a few scans).

### 3.2 Page body (the demo wiring)

```python
from pathlib import Path
import streamlit as st
from src.common.common import page_setup, save_params, show_linked_grid  # §4 below
from src.workflow.FileManager import FileManager
from src.view.grid import LayoutManager
from openms_insight import Table, LinePlot, Heatmap, SequenceView

params = page_setup()
DATA = Path("example-data/insight")
fm = FileManager(st.session_state.workspace, cache_path=Path(st.session_state.workspace, "cache"))
cache = str(Path(st.session_state.workspace, "cache", "insight"))

OPTIONS = ["Spectrum table", "Spectrum plot", "Peak map", "Sequence view"]
NAMES   = ["spectra_table", "spectrum_plot", "peak_map", "sequence_view"]

def builders():
    return {
        "spectra_table": lambda: Table(
            cache_id="demo_spectra", data_path=str(DATA/"spectra.parquet"),
            cache_path=cache, interactivity={"spectrum": "scan_id"},
            index_field="scan_id", default_row=0,
        ),
        "spectrum_plot": lambda: LinePlot(
            cache_id="demo_spectrum_plot", data_path=str(DATA/"peaks.parquet"),
            cache_path=cache, filters={"spectrum": "scan_id"},
            interactivity={"peak": "peak_id"}, x_column="mass", y_column="intensity",
            highlight_column="is_annotated", annotation_column="ion_label",
            title="MS/MS Spectrum",
        ),
        "peak_map": lambda: Heatmap(
            cache_id="demo_peak_map", data_path=str(DATA/"heat.parquet"),
            cache_path=cache, x_column="rt", y_column="mass", intensity_column="intensity",
            interactivity={"spectrum": "scan_id", "peak": "peak_id"}, title="Peak Map",
        ),
        "sequence_view": lambda: SequenceView(
            cache_id="demo_seq", sequence_data_path=str(DATA/"sequences.parquet"),
            peaks_data_path=str(DATA/"peaks.parquet"), cache_path=cache,
            filters={"spectrum": "scan_id"}, interactivity={"peak": "peak_id"},
            deconvolved=True, title="Fragment Coverage",
        ),
    }

DEFAULT_LAYOUT = [["spectra_table", "spectrum_plot"], ["peak_map", "sequence_view"]]

tab_view, tab_layout = st.tabs(["Viewer", "Layout Manager"])
lm = LayoutManager(OPTIONS, NAMES, store=fm, layout_id="demo_layout", session_prefix="demo")
with tab_layout:
    lm.render()
with tab_view:
    saved = lm.get_layout()
    layout, side_by_side = (saved if saved else ([DEFAULT_LAYOUT], False))
    show_linked_grid(layout, builders(), tool="demo", side_by_side=side_by_side)
save_params(params)
```

### 3.3 Register in `app.py`

Add to the `pages` dict (mirrors how the FLASHApp viewers are registered):

```python
"Visualization Template": [
    st.Page(Path("content", "visualization_template.py"),
            title="Linked Grid Demo", icon="🔗"),
],
```

---

## 4. `streamlit-template/src/common/common.py` — `show_linked_grid()` one-liner

Add a thin convenience over `render_linked_grid` that handles the **multi-experiment + side-by-side**
page concern (the part the oracle viewer pages hand-roll), so any template/FLASHApp viewer collapses
to one call. Keep `show_fig`/`show_table` untouched.

```python
# append to src/common/common.py
def show_linked_grid(layout, builders, *, tool, side_by_side=False,
                     grid_key="linked_grid", height=None, column_heights=None):
    """Render an N-experiment linked grid. `layout` is List[experiment]; each experiment is the
    nested rows list consumed by render_linked_grid. One independent StateManager per experiment
    (session_key f'{tool}__exp{i}') so experiments never cross-link. When exactly two experiments
    and side_by_side=True, render them in two st.columns; otherwise stack with st.divider()."""
    from src.view.grid import render_linked_grid
    import streamlit as st

    def _one(exp_idx, exp_layout, container):
        with container:
            render_linked_grid(
                exp_layout, builders, state_key=f"{tool}__exp{exp_idx}",
                grid_key=f"{grid_key}_{exp_idx}", height=height, column_heights=column_heights,
            )

    if len(layout) == 2 and side_by_side:
        c1, c2 = st.columns(2)
        _one(0, layout[0], c1); _one(1, layout[1], c2)
    else:
        for i, exp_layout in enumerate(layout):
            if i: st.divider()
            _one(i, exp_layout, st.container())
```

This is the "one-liner" the viewers call. Experiment selection (the `st.selectbox("choose
experiment", ...)` per experiment) stays in the viewer page because it is tool/data specific
(it needs the FileManager results list + display names); `show_linked_grid` only owns the
grid+side-by-side rendering. (The selectbox+grid pairing in the oracle is exactly this split.)

---

## 5. FLASHApp rebuild (from the frozen template)

### 5.1 `src/render/schema.py` (NEW) — FileManager caches -> Insight-ready tidy parquet

The oracle ships *wide, list-column, index-addressed* caches (one row per scan with array cells;
selection by positional `iloc`/`SignalPeaks[massIndex]`). Insight components want **tidy parquet
with stable value IDs** addressed by `filters`/`interactivity`. `schema.py` is the adapter: it
reads existing FileManager caches and writes derived tidy parquet (via `store_data`, so they live
in the same SQLite-indexed store and get a `result_path`). It is a **pure post-process** — it does
not touch `src/parse/*` producers.

Public API:

```python
# src/render/schema.py
def build_insight_caches(file_manager, dataset_id, tool, logger=None) -> None:
    """Read the oracle caches for (dataset_id, tool) and write the tidy parquet that the
    Insight builders (§5.2) consume via data_path=. Idempotent + cache-guarded: skip a target
    if its name_tag already exists (file_manager.result_exists) unless regenerate=True."""
```

Call site: append `build_insight_caches(file_manager, dataset_id, tool)` at the end of each parse
step (`parseDeconv`/`parseTnT`/`parseQuant` in `src/parse/*`, or right after them in `Workflow.py`),
OR lazily the first time a viewer loads a dataset (guarded by `result_exists`). Lazy-on-first-view
is recommended so re-processing isn't required for the migration.

**Stable IDs minted here** (deterministic, dataset-scoped): `scan_id` (= oracle scan-table `index`,
already 0..N), `mass_id` (per (scan, mass) — global running id), `peak_id` (per exploded signal/raw
peak — global running id), `protein_id` (= protein_df `index`), `tag_id` (per tag row), `feature_id`
(= FeatureGroupIndex). These become the `interactivity`/`filters` columns.

#### 5.1.1 Per-component tidy-parquet schemas (the data contract)

Mapping each oracle structure -> the parquet each Insight component consumes. Columns are the
*minimum* each component reads; carry extra display columns freely (render-time, uncached-hash).

**(a) Scan table** — oracle `scan_table` (already tidy). Component: `Table`.
`scans.parquet`: `scan_id:int(=index), Scan:int, MSLevel:int, RT:float, PrecursorMass:float, #Masses:int`.
Builder: `Table(interactivity={"scan": "scan_id"}, index_field="scan_id", default_row=0)`.
*Replaces oracle:* clicking a row set `scanIndex` (== the row's `index`); now sets selection
`scan` = `scan_id`.

**(b) Mass table** — oracle `mass_table` (one row/scan, list cells `MonoMass`,`SumIntensity`,
charges/isotopes/scores). Component: `Table`, filtered by scan. **Explode list cells to one row
per mass.** `masses.parquet`:
`scan_id:int, mass_id:int, mass_in_scan:int(0-based pos within scan), MonoMass:float,
SumIntensity:float, MinCharges:int, MaxCharges:int, MinIsotopes:int, MaxIsotopes:int,
CosineScore:float, SNR:float, QScore:float`.
Builder: `Table(filters={"scan": "scan_id"}, interactivity={"mass": "mass_id"}, index_field="mass_id")`.
*Replaces oracle:* `iloc[scanIndex]` row + frontend reading the list cells; `massIndex` ->
`mass_in_scan` is retained so 3D/spectrum overlays can still index a scan's mass arrays, and
`mass_id` is the cross-link value.

**(c) Deconvolved spectrum** — oracle `deconv_spectrum` (list `MonoMass`,`SumIntensity` per scan).
Component: `LinePlot` (default stick mode), filtered by scan. **Explode to one row per peak.**
`deconv_spectrum.parquet`: `scan_id:int, peak_id:int, MonoMass:float, SumIntensity:float`.
Builder: `LinePlot(filters={"scan": "scan_id"}, x_column="MonoMass", y_column="SumIntensity",
interactivity={"mass": "peak_id"})`.

**(d) Annotated / Augmented spectrum** — oracle `combined_spectrum` (deconv masses + `SignalPeaks`
nested cell + anno arrays). Two builders share this source:
- *Annotated Spectrum* = `LinePlot` over the **raw m/z** arrays (`MonoMass_Anno`/`SumIntensity_Anno`).
  `anno_spectrum.parquet`: `scan_id:int, peak_id:int, mz:float, intensity:float, is_signal:int`
  (explode `MonoMass_Anno`/`SumIntensity_Anno`; `is_signal` from membership in any `SignalPeaks`
  record's `peak_index` for that scan -> `highlight_column`).
- *Augmented Deconvolved Spectrum* = `LinePlot.tagger(...)` (top-down recipe; README §LinePlot
  modes). This mode consumes the **per-scan list-column frame as-is** (it does its own explode),
  so write `combined_tagger.parquet` = one row per scan with list columns:
  `scan_id:int, MonoMass:list<f64>, SumIntensity:list<f64>, SignalPeaks:list<list<list<f64>>>,
   Mzs:list<f64>, MzIntensities:list<f64>`.
  Builder: `LinePlot.tagger(filters={"spectrum":"scan_id"}, x_column="MonoMass",
  y_column="SumIntensity", signal_peaks_column="SignalPeaks", mz_column="Mzs",
  mz_intensity_column="MzIntensities", interactivity={"tagger_mass":"peak_id"},
  tag_identifier="tag")`. (`SignalPeaks[mass][peak] = [peak_index, mz, intensity, charge]` — exactly
  the inner record produced by `masstable._compute_peak_cells`, confirmed in oracle.)

**(e) 3D S/N plot ("Precursor Signals")** — oracle `threedim_SN_plot` (per scan: `SignalPeaks`,
`NoisyPeaks` nested cells; `update.py` picks `SignalPeaks[massIndex]` then renders points
`[peak_index, mz, intensity, charge]`). Component: `Plot3D`. **Explode the nested cells fully to
one row per point**, tagged Signal/Noise. `precursor_signals.parquet`:
`scan_id:int, mass_in_scan:int, peak_id:int, mz:float, charge:int, intensity:float, series:str("Signal"|"Noise")`.
Builder:
```python
Plot3D(filters={"scan": "scan_id", "mass": "mass_in_scan"},
       filter_defaults={"scan": -1},
       x_column="mz", y_column="charge", z_column="intensity",
       category_column="series", category_colors={"Signal":"#3366CC","Noise":"#DC3912"},
       title="Precursor Signals")
```
*Replaces oracle:* `scanIndex`+`massIndex` two-level positional filter -> value filters
`scan`(=scan_id) + `mass`(=mass_in_scan), exactly mirroring `update.py`'s
`SignalPeaks[mass_index]` slice but value-based. (README Plot3D example uses precisely this
`filters={'spectrum':'scan','mass':'mass_index'}` shape.)

**(f) Heatmaps (Raw/Deconv MS1/MS2)** — oracle builds a *full* `ms{1,2}_{deconv,raw}_heatmap`
plus precomputed compression levels and re-downsamples on zoom (`update.py::render_heatmap`,
`compression.downsample_heatmap`). Component: `Heatmap` (does its **own** multi-resolution
downsampling + zoom). So we **drop the precomputed `_<size>` levels and the bespoke
`render_heatmap`/`downsample_heatmap` zoom path entirely** and feed the full frame:
`ms{lvl}_{kind}_heatmap.parquet`: `rt:float, mass:float, intensity:float` (already the oracle's
full-resolution schema — `getMSSignalDF` aliases `mz_array->mass`, `intensity_array->intensity`).
Builder: `Heatmap(x_column="rt", y_column="mass", intensity_column="intensity", title=...)`.
*Replaces oracle:* `xRange/yRange` zoom + `render_heatmap` cache -> Insight's internal zoom +
multi-resolution cache. **No schema.py work needed for heatmaps** beyond pointing the builder at
the existing full-resolution `.pq` via `result_path` (these are already tidy). The `_<size>`
caches simply stop being produced (optional cleanup in `parse/deconv.py`, not required to delete).

**(g) Score Distribution / FDR plot** — oracle `density_target`/`density_decoy` (and the
`density_id_*` pair for tnt), each a `{x,y}` KDE DataFrame. Component: `LinePlot.density(...)`
(README density mode). **Concatenate the two into one long/tidy frame with a category column.**
`qscore_density.parquet`: `x:float (qscore/qvalue), y:float (density), group:str("target"|"decoy")`.
Builder: `LinePlot.density(x_column="x", y_column="y", category_column="group",
target_value="target", decoy_value="decoy", title="Score Distribution")`. (deconv uses
`density_{target,decoy}`; tnt uses `density_id_{target,decoy}` -> same tidy output.)

**(h) Protein table** — oracle `protein_dfs` (already tidy pandas). Component: `Table`.
`proteins.parquet`: `protein_id:int(=index), accession:str, description:str, sequence:str,
length:int, ProteoformMass:float, ProteoformLevelQvalue:float, Scan:int, ...`.
Builder: `Table(interactivity={"protein": "protein_id"}, index_field="protein_id", default_row=0)`.
*Replaces oracle:* row click set `proteinIndex` -> selection `protein` = `protein_id`.

**(i) Tag table** — oracle `tag_dfs` (one row per (tag,proteoform), sorted by `Scan`;
`update.py` resolves the selected `proteinIndex` -> scan via `proteoform_scan_map`, filters by
`Scan`, stamps `ProteinIndex`). Component: `Table` filtered by protein. Bake the
proteoform-scan resolution **into the parquet at build time** (no runtime `scan_map`):
`tags.parquet`: `tag_id:int, protein_id:int (resolved proteoform index, via
scan_resolution.build_proteoform_scan_map + tag_resolution mapping), scan_id:int, Scan:int,
TagSequence:str, StartPos:int, EndPos:int, Length:int, Score:float, mzs:str`.
Builder: `Table(filters={"protein": "protein_id"}, interactivity={"tag": "tag_id"},
index_field="tag_id")`. *Replaces oracle:* the entire `proteoform_scan_map` + `Scan`-pushdown +
`ProteinIndex`-stamp dance in `filter_data` collapses to a precomputed `protein_id` column +
a value filter.

**(j) Sequence view** — oracle: FLASHDeconv computes fragments at render time from a sequence in
the `('sequence','sequence')` cache (`update.py::get_sequence` + `render_sequence_data`);
FLASHTnT reads a per-proteoform `sequence_data` parquet (`sequence_data_store.py`, one row per
proteoform with fragment-mass list-of-lists, coverage, modifications) and `load_entry(pid)`.
Component: `SequenceView` (it does fragment matching itself from sequence + peaks). Two cases:
- *FLASHDeconv* (single global sequence): build `seq_deconv.parquet` with one row per scan:
  `scan_id:int, sequence:str, precursor_charge:int` (sequence is the global input sequence,
  charge from precursor). Peaks = the deconv spectrum long frame (`deconv_spectrum.parquet`,
  neutral masses -> `deconvolved=True`). Builder:
  `SequenceView(sequence_data_path="seq_deconv.parquet", peaks_data_path="deconv_spectrum.parquet",
  filters={"scan":"scan_id"}, interactivity={"mass":"peak_id"}, deconvolved=True)`.
- *FLASHTnT* (per-proteoform): build `seq_tnt.parquet` one row per proteoform:
  `protein_id:int, sequence:str, precursor_charge:int, coverage:list<f64>,
  proteoform_start:int, proteoform_end:int` (coverage/start/end straight from
  `sequence_data_store` entry; SequenceView's `coverage_column`/`proteoform_start_column`/
  `proteoform_end_column` opt-ins consume them). Peaks = per-scan deconv masses resolved by the
  proteoform's scan. Builder:
  `SequenceView(sequence_data_path="seq_tnt.parquet", peaks_data_path=..., filters={"protein":
  "protein_id"}, interactivity={"mass":"peak_id"}, deconvolved=True, coverage_column="coverage",
  proteoform_start_column="proteoform_start", proteoform_end_column="proteoform_end")`.
  *Note:* the rich theoretical-fragment list-of-lists the oracle precomputed
  (`getFragmentDataFromSeq`) is **no longer needed** — SequenceView enumerates + matches ion types
  itself from `sequence` + `annotation_config={"ion_types": settings["ion_types"], "tolerance":
  settings["tolerance"]}` (read from the oracle `settings` cache). The `sequence_data_store.py`
  table can stay as a coverage/modification source only, or be replaced by `seq_tnt.parquet`.
- *Internal Fragment Map* is **disabled** in the oracle TnT manager (commented out) and the deconv
  `internal_fragment_map` branch is dead code — do not rebuild it; if ever re-enabled, it maps to
  `SequenceView(internal_fragments=True)`.

**(k) FLASHQuant** — oracle `quant_dfs` (one row per FeatureGroup: scalar columns + list columns
`Charges/IsotopeIndices/CentroidMzs/RTs/MZs/Intensities`, each a list of comma-joined strings per
trace). Components: `Table` (feature list) `<->` `Plot3D` (the feature's traces in 3D). Build two:
- `quant_features.parquet` (tidy scalars): `feature_id:int(=FeatureGroupIndex), MonoisotopicMass,
  AverageMass, StartRT, EndRT, ApexRT, FeatureGroupQuantity, AllAUC, MinCharge, MaxCharge,
  MostAbundantFeatureCharge, IsotopeCosineScore`. Builder: `Table(interactivity={"feature":
  "feature_id"}, index_field="feature_id", default_row=0)`.
- `quant_traces.parquet` (long, the comma-split explode): for each feature, each trace, split the
  comma-joined `MZs`/`RTs`/`Intensities` strings to one row per point:
  `feature_id:int, charge:int, isotope:int, centroid_mz:float, rt:float, mz:float, intensity:float`.
  Builder: `Plot3D(filters={"feature":"feature_id"}, filter_defaults={"feature":-1},
  x_column="rt", y_column="mz", z_column="intensity", category_column="charge", title="Feature Traces")`.
  *Replaces oracle:* the bespoke `FLASHQuantView` Vue component -> `Table<->Plot3D` linked pair.

> **Explode/long-format helpers** in `schema.py`: `_explode_list_cols(df, by, list_cols, id_name)`
> (polars `explode` + running id), `_explode_nested_signal_peaks(df, col, series_label)` (two-level
> `explode` for `SignalPeaks`/`NoisyPeaks` -> `[peak_index,mz,intensity,charge]` rows),
> `_comma_split_long(df, cols)` (str.split("," ) + `explode` for quant traces),
> `_kde_to_long(target_df, decoy_df)` (concat with `group` col). All polars-lazy, written via
> `file_manager.store_data(..., row_group_size=...)` so Insight pushdown stays efficient.

### 5.2 The builders factory (`comp_name -> () -> BaseComponent(data_path=...)`)

`render.py` is **repurposed** from "grid render loop" to "FLASHApp's builder factory" (the grid
loop itself is deleted — §5.5 — and the page imports the frozen template grid). New `render.py`:

```python
# src/render/render.py  (post-migration: builders only; no grid loop)
from pathlib import Path
from openms_insight import (Table, LinePlot, Heatmap, Plot3D, SequenceView)

def make_builders(file_manager, dataset_id, tool, settings=None):
    """Return {comp_name: () -> BaseComponent} for one (tool, dataset). Each factory closes over
    dataset_id + file_manager + an Insight cache dir, and uses file_manager.result_path(...) to
    feed data_path=. cache_id is f'{tool}__{dataset_id}__{comp_name}' so caches are per-dataset
    (this is the oracle's 'dataset changed -> reset' guarantee, expressed via cache_id)."""
    p   = lambda tag: file_manager.result_path(dataset_id, tag)        # parquet path
    cid = lambda name: f"{tool}__{dataset_id}__{name}"
    cache = str(Path(file_manager.cache_path, "insight"))

    B = {
      "scan_table":     lambda: Table(cache_id=cid("scan_table"), data_path=p("scans"),
                            cache_path=cache, interactivity={"scan":"scan_id"},
                            index_field="scan_id", default_row=0, title="Scan Table"),
      "mass_table":     lambda: Table(cache_id=cid("mass_table"), data_path=p("masses"),
                            cache_path=cache, filters={"scan":"scan_id"},
                            interactivity={"mass":"mass_id"}, index_field="mass_id",
                            title="Mass Table"),
      "deconv_spectrum":lambda: LinePlot(cache_id=cid("deconv_spectrum"),
                            data_path=p("deconv_spectrum"), cache_path=cache,
                            filters={"scan":"scan_id"}, interactivity={"mass":"peak_id"},
                            x_column="MonoMass", y_column="SumIntensity",
                            title="Deconvolved Spectrum"),
      "anno_spectrum":  lambda: LinePlot(cache_id=cid("anno_spectrum"), data_path=p("anno_spectrum"),
                            cache_path=cache, filters={"scan":"scan_id"},
                            interactivity={"mass":"peak_id"}, x_column="mz", y_column="intensity",
                            highlight_column="is_signal", title="Annotated Spectrum"),
      "combined_spectrum": lambda: LinePlot.tagger(cache_id=cid("combined_spectrum"),
                            data_path=p("combined_tagger"), cache_path=cache,
                            filters={"spectrum":"scan_id"}, interactivity={"tagger_mass":"peak_id"},
                            x_column="MonoMass", y_column="SumIntensity",
                            signal_peaks_column="SignalPeaks", mz_column="Mzs",
                            mz_intensity_column="MzIntensities", tag_identifier="tag",
                            title="Augmented Deconvolved Spectrum"),
      "3D_SN_plot":     lambda: Plot3D(cache_id=cid("3D_SN_plot"), data_path=p("precursor_signals"),
                            cache_path=cache, filters={"scan":"scan_id","mass":"mass_in_scan"},
                            filter_defaults={"scan":-1}, x_column="mz", y_column="charge",
                            z_column="intensity", category_column="series",
                            category_colors={"Signal":"#3366CC","Noise":"#DC3912"},
                            title="Precursor Signals"),
      "ms1_deconv_heat_map": lambda: Heatmap(cache_id=cid("ms1_deconv_heat_map"),
                            data_path=p("ms1_deconv_heatmap"), cache_path=cache,
                            x_column="rt", y_column="mass", intensity_column="intensity",
                            interactivity={"scan":"scan_id"}, title="Deconvolved MS1 Heatmap"),
      "ms2_deconv_heat_map": lambda: Heatmap(cache_id=cid("ms2_deconv_heat_map"),
                            data_path=p("ms2_deconv_heatmap"), cache_path=cache, x_column="rt",
                            y_column="mass", intensity_column="intensity",
                            title="Deconvolved MS2 Heatmap"),
      "ms1_raw_heatmap":lambda: Heatmap(cache_id=cid("ms1_raw_heatmap"), data_path=p("ms1_raw_heatmap"),
                            cache_path=cache, x_column="rt", y_column="mass",
                            intensity_column="intensity", title="Raw MS1 Heatmap"),
      "ms2_raw_heatmap":lambda: Heatmap(cache_id=cid("ms2_raw_heatmap"), data_path=p("ms2_raw_heatmap"),
                            cache_path=cache, x_column="rt", y_column="mass",
                            intensity_column="intensity", title="Raw MS2 Heatmap"),
      "fdr_plot":       lambda: LinePlot.density(cache_id=cid("fdr_plot"), data_path=p("qscore_density"),
                            cache_path=cache, x_column="x", y_column="y", category_column="group",
                            target_value="target", decoy_value="decoy", title="Score Distribution"),
      "id_fdr_plot":    lambda: LinePlot.density(cache_id=cid("id_fdr_plot"),
                            data_path=p("qscore_density_id"), cache_path=cache, x_column="x",
                            y_column="y", category_column="group", target_value="target",
                            decoy_value="decoy", title="Score Distribution"),
      "protein_table":  lambda: Table(cache_id=cid("protein_table"), data_path=p("proteins"),
                            cache_path=cache, interactivity={"protein":"protein_id"},
                            index_field="protein_id", default_row=0, title="Protein Table"),
      "tag_table":      lambda: Table(cache_id=cid("tag_table"), data_path=p("tags"), cache_path=cache,
                            filters={"protein":"protein_id"}, interactivity={"tag":"tag_id"},
                            index_field="tag_id", title="Tag Table"),
      "sequence_view":  lambda: _sequence_view(file_manager, dataset_id, tool, cid, cache, p, settings),
      "quant_visualization": lambda: Table(cache_id=cid("quant_features"), data_path=p("quant_features"),
                            cache_path=cache, interactivity={"feature":"feature_id"},
                            index_field="feature_id", default_row=0, title="Features"),
      "quant_traces_3d": lambda: Plot3D(cache_id=cid("quant_traces"), data_path=p("quant_traces"),
                            cache_path=cache, filters={"feature":"feature_id"},
                            filter_defaults={"feature":-1}, x_column="rt", y_column="mz",
                            z_column="intensity", category_column="charge", title="Feature Traces"),
    }
    return B
```

`_sequence_view(...)` branches on `tool` to pick the deconv vs tnt SequenceView wiring described
in §5.1.1(j) (deconv: global sequence from `('sequence','sequence')`; tnt: per-proteoform
`seq_tnt.parquet` + coverage/proteoform columns + `annotation_config` from the `settings` cache).

**StateManager — one per (tool, experiment).** The grid creates it from `state_key`. The viewer
passes `state_key=f"{tool}__{experiment_id}"` (via `show_linked_grid`'s `tool=` -> `f'{tool}__exp{i}'`,
combined with the selected experiment id baked into builders' `cache_id`). Net effect: experiment A
and experiment B shown together have independent selections and independent component caches —
exactly the oracle's `state_tracker[tool][identifier]` isolation, now provided by Insight.

### 5.3 The OLD index-based selection -> Insight value-based interactivity (oracle map, cite `update.py`)

| Oracle (`update.py` / `filter_data`) | Insight (`filters`/`interactivity` + StateManager) |
|---|---|
| `selection_store['scanIndex']`; `per_scan_data.iloc[scanIndex:scanIndex+1]` | selection `scan` = `scan_id`; every per-scan panel `filters={"scan":"scan_id"}` |
| `selection_store['massIndex']`; `SignalPeaks[massIndex]`/`NoisyPeaks[massIndex]` | selection `mass` = `mass_in_scan` (Plot3D) / `mass_id` (Mass Table); `filters={"mass": ...}` |
| `proteinIndex` -> `proteoform_scan_map[proteinIndex]` -> filter `Scan`; stamp `ProteinIndex` (Tag/Seq, tnt) | selection `protein` = `protein_id`; `tags.parquet`/`seq_tnt.parquet` carry a precomputed `protein_id` column; `filters={"protein":"protein_id"}` (scan-map resolution moved to build time) |
| heatmap `selection_store['heatmap_*'] = {xRange,yRange}` -> `render_heatmap` re-downsample | Heatmap internal zoom + multi-resolution cache (no Python zoom path; per-instance `zoom_identifier`) |
| `get_sequence(selection_store)` + `render_sequence_data` (deconv) | `SequenceView(filters={"scan":"scan_id"}, deconvolved=True)` matches fragments itself |
| `load_entry(sequence_data_ds, proteinIndex)` (tnt) | `SequenceView(filters={"protein":"protein_id"}, coverage_column=..., proteoform_*_column=...)` |
| `StateTracker` (counter+id, per identifier) | `StateManager(session_key=state_key)` (identical counter+id reconcile, `update_from_vue`) |
| cleared selection echoed as `None` (render.py drop-None) | StateManager `clear_selection`/`set_selection(None)` semantics (already handled) |

The cross-link chains the deconv viewer must preserve (oracle): **scan -> mass -> spectrum -> 3D**
(Scan Table click sets `scan`; Mass Table + spectra + 3D filter by `scan`; Mass Table click sets
`mass`; 3D + spectrum highlight by `mass`); **protein -> tag -> sequence** (tnt: Protein Table sets
`protein`; Tag Table + Sequence View filter by `protein`; Tag/peak click sets `tag`/`mass`);
**heatmap zoom** (now component-internal). All expressible purely through the identifier vocabulary
above — no Python per-rerun filtering.

### 5.4 The three viewer pages — each shrinks to: pick experiment(s) -> load layout -> render

Reference (FLASHDeconvViewer.py, post-migration ~35 lines; TnT/Quant analogous):

```python
import streamlit as st
from pathlib import Path
from src.common.common import page_setup, save_params, show_linked_grid
from src.workflow.FileManager import FileManager
from src.render.render import make_builders
from src.render.schema import build_insight_caches

DEFAULT_LAYOUT = [["ms1_deconv_heat_map"], ["scan_table","mass_table"],
                  ["anno_spectrum","deconv_spectrum"], ["3D_SN_plot"]]

params = page_setup()
fm = FileManager(st.session_state.workspace, Path(st.session_state.workspace, "cache"))
results = fm.get_results_list(["threedim_SN_plot"])
if not results:
    st.error("No results to show yet. Please run a workflow first!"); st.stop()

names  = [fm.get_display_name(r) for r in results]
to_id  = {fm.get_display_name(r): r for r in results}

saved  = fm.get_results("layout","layout", partial=True).get("layout") if \
         fm.result_exists("layout","layout") else None
layout, side_by_side = (saved["layout"], saved["side_by_side"]) if saved else ([DEFAULT_LAYOUT], False)
# append sequence_view to default if a sequence is set (oracle parity)
if fm.result_exists("sequence","sequence") and not saved:
    layout = [DEFAULT_LAYOUT + [["sequence_view"]]]

# one experiment selector per layout slot (tool/data-specific -> stays in the page)
chosen = []
for i in range(len(layout)):
    label = "choose experiment" if i == 0 else None
    sel = st.selectbox(label or "choose experiment", names, key=f"deconv_exp_{i}")
    chosen.append(to_id[sel])

# lazily build Insight caches for chosen datasets (idempotent / cache-guarded)
for ds in set(chosen):
    build_insight_caches(fm, ds, "flashdeconv")

# builders for the (first) chosen dataset per experiment slot; multi-exp uses per-slot builders
def builders_for(ds): return make_builders(fm, ds, "flashdeconv",
                                           settings=None)
# render: show_linked_grid drives side-by-side / stacked + one StateManager per experiment
if len(layout) == 2 and side_by_side:
    show_linked_grid([layout[0]], builders_for(chosen[0]), tool=f"flashdeconv_{chosen[0]}",
                     side_by_side=False)
    show_linked_grid([layout[1]], builders_for(chosen[1]), tool=f"flashdeconv_{chosen[1]}",
                     side_by_side=False)
else:
    for i, exp_layout in enumerate(layout):
        if i: st.divider()
        show_linked_grid([exp_layout], builders_for(chosen[i]),
                         tool=f"flashdeconv_{chosen[i]}", side_by_side=False)
save_params(params)
```

- **FLASHTnTViewer.py**: same shape, `tool="flashtnt"`, `DEFAULT_LAYOUT =
  [["protein_table"],["sequence_view"],["tag_table"],["combined_spectrum"]]`, layout cache
  `("flashtnt_layout","layout")`, results gate `["protein_dfs"]`, `settings` passed to
  `make_builders` for SequenceView ion-types/tolerance.
- **FLASHQuantViewer.py**: simplest — gate `["quant_dfs"]`, FileManager rooted at
  `workspace/flashquant/cache` (oracle keeps this), default layout
  `[["quant_visualization","quant_traces_3d"]]` (feature Table `<->` 3D traces, the quant recipe),
  no LayoutManager needed.

Layout Manager pages become one-liners too:
```python
# content/FLASHDeconv/FLASHDeconvLayoutManager.py
from src.view.grid import LayoutManager
from src.workflow.FileManager import FileManager
from src.common.common import page_setup, save_params
# ... DECONV_OPTIONS / DECONV_NAMES constants live here (the only tool-specific bit) ...
params = page_setup()
fm = FileManager(st.session_state.workspace, Path(st.session_state.workspace,"cache"))
lm = LayoutManager(DECONV_OPTIONS, DECONV_NAMES, store=fm, layout_id="layout",
                   session_prefix="deconv", title="Layout Manager")
if fm.result_exists("sequence","sequence"):
    lm.add_options(["Sequence view (Mass table needed)"], ["sequence_view"])
lm.render(); save_params(params)
```

### 5.5 What gets DELETED / changed in FLASHApp

- **Delete** `src/render/components.py` (Vue declaration + all `FlashViewer*`/`PlotlyHeatmap`/
  `Tabulator`/`SequenceView`/`Plotly3Dplot`/... wrapper classes) — replaced by Insight components.
- **Delete** `src/render/initialize.py` (per-panel cache loading) — replaced by §5.1 schema +
  §5.2 builders feeding `data_path=`.
- **Delete** `src/render/update.py` (index-based `update_data`/`filter_data`/`render_heatmap`/
  `get_sequence`/`render_sequence_data`) — replaced by Insight `filters`/`interactivity` + each
  component's own preprocessing.
- **Delete** `src/render/StateTracker.py` — replaced by `openms_insight.StateManager`.
- **Delete** `src/render/render.py`'s grid loop (`render_grid`/`render_component`) — the grid now
  comes from the frozen template; `render.py` is repurposed to the builders factory (§5.2).
- **Optionally retire** `src/render/compression.py` zoom/`downsample_heatmap` and the producer's
  `_<size>` compression-level outputs (Heatmap downsamples itself). `compute_compression_levels`
  can go once initialize.py is gone. Safe to leave in `parse/deconv.py` until cleanup.
- **`util.py::hash_complex`** (used only by `render_component`) -> delete with the loop.
- **`js-component/`**: stop using it. Remove the `path=build_dir` declaration (in deleted
  components.py) and the submodule from build/CI (`Dockerfile*`, `.gitmodules` if present, the
  `js-component/dist` packaging in `run_app_temp.spec`). Insight ships its own Vue bundle.
- **`requirements.txt`**: add `openms-insight` (pin a version, e.g. `openms-insight==0.1.11`).
  Insight pulls polars/pyarrow; keep existing pins. Drop any js-build deps that existed only for
  the local component.
- **Keep**: `src/workflow/FileManager.py` (now mirrors the template's; see §6 note),
  `src/render/scan_resolution.py` + `tag_resolution.py` + `sequence.py` + `sequence_data_store.py`
  (now *build-time* helpers used by `schema.py` to mint `protein_id`/coverage), `src/parse/*`
  producers (unchanged; `schema.py` post-processes their output).

---

## 6. Non-divergence — FLASHApp uses the template's `grid.py` UNCHANGED

**Mechanism (recommended): git submodule of `streamlit-template` + thin re-export shim.**
`nondivergence.py` normalizes (strip trailing whitespace, drop blank lines and full-line
comments) then SHA-256-compares the two registered files. So the FLASHApp side must be the
*same source text* as the template's frozen `grid.py` (comments/blank-lines aside). The cleanest
way that the gate accepts and that avoids stale copies:

1. Add `OpenMS/streamlit-template` as a git submodule at `FLASHApp/streamlit-template/` (pinned to
   the frozen commit).
2. Create `FLASHApp/src/view/grid.py` as the registered FLASHApp path whose **content is byte-identical
   to the template's** `src/view/grid.py`. Two acceptable implementations:
   - **(preferred) vendored copy kept in sync by CI:** a tiny `make sync-grid` /
     pre-commit step copies `streamlit-template/src/view/grid.py` -> `src/view/grid.py`. The
     normalized-hash gate then trivially passes, and FLASHApp imports `from src.view.grid import
     render_linked_grid, LayoutManager` with no path gymnastics.
   - **(alt) symlink:** `src/view/grid.py -> ../../streamlit-template/src/view/grid.py`. Same bytes
     by construction; works on Linux/CI (the deployment target). The vendored copy is safer across
     Windows packaging (`run_app_temp.spec`), so prefer it.

Either way the **registered file pair is identical content**, so `_normalized_hash(a) ==
_normalized_hash(b)` and the gate is GREEN. The submodule guarantees the template source is
present locally for the hash comparison and pins the exact frozen version.

**Register the pair** in `migration/units.yaml -> meta.nondivergence_pairs` (uncomment + set):

```yaml
  nondivergence_pairs:
    - [/home/user/FLASHApp/src/view/grid.py, /home/user/FLASHApp/streamlit-template/src/view/grid.py]
```

(With the submodule, the template path resolves *inside FLASHApp*, so the gate is self-contained and
does not depend on a sibling checkout. If the submodule route is rejected, point the second element
at `/home/user/streamlit-template/src/view/grid.py` and keep `src/view/grid.py` a vendored copy.)

**Only `grid.py` is the frozen, non-divergent unit.** `common.py::show_linked_grid`,
`FileManager.py`, and `visualization_template.py` are template *features* FLASHApp may mirror but
are not byte-frozen (FLASHApp keeps its own richer FileManager; the template's is the ported
subset). The single source of truth that must never fork is `grid.py` (the grid loop + LayoutManager).

---

## 7. Build / implementation order (template first -> freeze -> FLASHApp rebuild)

1. **Template `src/workflow/FileManager.py`** — port the FLASHApp store API + add `as_path=` /
   `result_path` (§2). Unit-test store -> `result_path` -> file exists.
2. **Template `src/view/grid.py`** — `render_linked_grid` + `LayoutManager` (§1). Unit-test the
   render loop with stub builders (assert ≤3 columns, per-cell keys, shared StateManager) and the
   LayoutManager trim/expand/validate/JSON round-trip against the oracle's behavior.
3. **Template `src/common/common.py`** — add `show_linked_grid` (§4).
4. **Template `content/visualization_template.py`** + `example-data/insight/*.parquet` + register
   in `app.py` (§3). Smoke: page parses (the `template-app-smoke` gate) and renders the 4-panel
   linked grid + LayoutManager + side-by-side over example data.
5. **FREEZE `grid.py`**; set up the submodule + vendored copy/symlink in FLASHApp; fill
   `units.yaml meta.nondivergence_pairs`; confirm `python migration/nondivergence.py` is GREEN.
6. **FLASHApp `src/render/schema.py`** (§5.1) — adapters + the per-component tidy parquet; unit-test
   each explode against a golden (reuse `reconstruct_all` from `sequence_data_store` for seq parity;
   compare exploded peak counts to oracle `SignalPeaks[mass]` lengths).
7. **FLASHApp `src/render/render.py`** -> builders factory (§5.2); delete the grid loop.
8. **FLASHApp viewers + layout managers** rebuilt (§5.4): `FLASHDeconvViewer.py`,
   `FLASHTnTViewer.py`, `FLASHQuantViewer.py`, both `FLASH*LayoutManager.py`. Smoke gates:
   `flashapp-app-smoke` (FLASHDeconvViewer parses) + manual per-panel + cross-link check.
9. **Delete** `components.py`/`initialize.py`/`update.py`/`StateTracker.py` + grid loop + js-component
   usage (§5.5); add `openms-insight` to `requirements.txt`.
10. Run the Phase-3 gates (`nondivergence`, both app-smokes) + the three critics
    (template / original-parity / final) per `units.yaml`.

---

## Appendix A — Quick reference: oracle cache -> Insight component -> tidy parquet

| comp_name | oracle cache(s) | Insight component | tidy parquet (key cols) | filters / interactivity |
|---|---|---|---|---|
| scan_table | `scan_table` | Table | `scans` (scan_id,Scan,MSLevel,RT,PrecursorMass,#Masses) | — / `scan`=scan_id |
| mass_table | `mass_table` | Table | `masses` (scan_id,mass_id,mass_in_scan,MonoMass,SumIntensity,charges,scores) | `scan` / `mass`=mass_id |
| deconv_spectrum | `deconv_spectrum` | LinePlot | `deconv_spectrum` (scan_id,peak_id,MonoMass,SumIntensity) | `scan` / `mass`=peak_id |
| anno_spectrum | `combined_spectrum` | LinePlot | `anno_spectrum` (scan_id,peak_id,mz,intensity,is_signal) | `scan` / `mass`=peak_id |
| combined_spectrum | `combined_spectrum` | LinePlot.tagger | `combined_tagger` (scan_id, list:MonoMass/SumIntensity/SignalPeaks/Mzs/MzIntensities) | `spectrum` / `tagger_mass` |
| 3D_SN_plot | `threedim_SN_plot` | Plot3D | `precursor_signals` (scan_id,mass_in_scan,peak_id,mz,charge,intensity,series) | `scan`+`mass` / — |
| ms{1,2}_{deconv,raw}_heatmap | `ms*_{deconv,raw}_heatmap` (full) | Heatmap | reuse existing (rt,mass,intensity) | — / (optional `scan`) |
| fdr_plot / id_fdr_plot | `density_{target,decoy}` / `density_id_*` | LinePlot.density | `qscore_density` / `qscore_density_id` (x,y,group) | — / — |
| protein_table | `protein_dfs` | Table | `proteins` (protein_id,accession,description,sequence,length,...) | — / `protein`=protein_id |
| tag_table | `tag_dfs` (+`proteoform_scan_map`) | Table | `tags` (tag_id,protein_id,scan_id,Scan,TagSequence,StartPos,EndPos,Length,Score,mzs) | `protein` / `tag`=tag_id |
| sequence_view (deconv) | `sequence`,`sequence_view` | SequenceView | `seq_deconv` (scan_id,sequence,precursor_charge) + peaks=`deconv_spectrum` | `scan` / `mass`=peak_id |
| sequence_view (tnt) | `sequence_data`,`settings` | SequenceView | `seq_tnt` (protein_id,sequence,charge,coverage,proteoform_start/end) + peaks | `protein` / `mass`=peak_id |
| quant_visualization | `quant_dfs` | Table | `quant_features` (feature_id,masses,RTs,quant,charges,score) | — / `feature`=feature_id |
| quant_traces_3d | `quant_dfs` | Plot3D | `quant_traces` (feature_id,charge,isotope,centroid_mz,rt,mz,intensity) | `feature` / — |

## Appendix B — `SignalPeaks` inner record (confirmed from `masstable._compute_peak_cells`)

`SignalPeaks` / `NoisyPeaks` are `list[mass_idx] -> list[peak] -> [peak_index, mz, intensity, charge]`
(all float64). This is exactly the structure `LinePlot.tagger(signal_peaks_column=...)` documents and
`Plot3D` consumes after a two-level explode (mass_idx -> `mass_in_scan`, peak -> a row). The oracle's
`update.py` selected `SignalPeaks[massIndex]`; the migration explodes ahead of time and filters by
`mass_in_scan` value instead.
