# FLASHApp → OpenMS-Insight migration

Status doc + parity contract for moving every FLASHApp visualization off the
bespoke `flash_viewer_grid` (`src/render/*` + the `openms-streamlit-vue-component`
submodule) and onto reusable [OpenMS-Insight](https://github.com/t0mdavid-m/OpenMS-Insight)
components. Branch (both repos): `claude/peaceful-mayer-YqiXZ`.

## Why

`flash_viewer_grid` is a single mega-component that draws the whole synchronized
grid and manages cross-plot selection internally (`StateTracker.py`, `update.py`).
OpenMS-Insight already re-implements most of these as individual components
sharing a `StateManager` (a functional superset of `StateTracker`) and one Vue
bundle. We close the gap by porting the few missing visualizations *into*
OpenMS-Insight, then compose them with native Streamlit layout.

## Component mapping

| FLASHApp | OpenMS-Insight | State |
|---|---|---|
| `TabulatorScan/Mass/Protein/TagTable` | `Table` | **reuse** |
| `PlotlyLineplot` (deconv/annotated) | `LinePlot` | **reuse** |
| `PlotlyLineplotTagger` (augmented) | `LinePlot` + overlay | **done (overlay); tag-geometry pending**¹ |
| `PlotlyHeatmap` (MS1/2 deconv/raw) | `Heatmap` | **reuse** |
| `SequenceView` | `SequenceView` (+coverage/fixed-mods/ion-types) | **done** ✅ |
| `Plotly3Dplot` (Precursor S/N) | `Scatter3D` | **done** ✅ |
| `FDRPlotly` (target/decoy KDE) | `DensityPlot` | **done** ✅ |
| `FLASHQuantView` (traces) | `FeatureView` | **done** ✅ |
| `InternalFragmentMap` (disabled) | `InternalFragmentMap` | **done** ✅ |

¹ The reusable overlay-series primitive (deconv sticks + raw/annotated peaks) is
done and tested. The tag-annotation **geometry** (charge-state buttons +
inter-residue amino-acid arrows driven by `tagIndex`/`tagData`) is tightly
coupled to the per-scan-array data model and is best validated by driving the
app, so it lands with the Phase-2 TnT wiring rather than ported blind.

## State-key mapping (StateTracker → OpenMS-Insight identifier→column)

- `scanIndex` → scan-table `interactivity={'scanIndex':'index'}`; spectra / mass
  table / sequence `filters={'scanIndex':'index'}`
- `massIndex` → mass-table `interactivity={'massIndex':'mass_id'}`; 3D plot
  `filters={'massIndex':'mass_id'}`
- `proteinIndex` → protein-table `interactivity`; sequence / tag / spectrum
  `filters`
- `tagIndex`/`tagData`, `AApos`, `sequenceOut` → sequence/tag interactivity
- `heatmap_deconv/_deconv2/_raw/_raw2` → per-`Heatmap` `zoom_identifier` (one per
  heatmap) — the four bespoke range keys collapse into each heatmap's own zoom

## The critical data-model change (highest risk)

FLASHApp filters by **row index** (`iloc[scanIndex]`) over **arrays-per-scan**
(`MonoMass[]`, `SumIntensity[]`). OpenMS-Insight filters by **column value** over
**long format** (one row per peak). `src/parse/long_format.py` bridges this and
is fully unit-tested (`tests/test_long_format.py`):

- `explode_spectrum_long` — per-scan mass/intensity arrays → one row per peak
  with explicit `index` + `mass_id` (assigned **before** any intensity filter so
  `massIndex` maps to the original array position). Filtering `index == k`
  reproduces the old `iloc[k]`.
- `explode_combined_spectrum_long` — deconv + annotated → primary + overlay.
- `explode_signal_peaks_long` — nested `SignalPeaks`/`NoisyPeaks` → Scatter3D
  long format (`index, mass_id, mz, charge, intensity, kind`).
- `density_series_long` — precomputed target/decoy curves → DensityPlot long
  `{series, x, y}`.

These are **additive**: the existing render pipeline is untouched, so the old and
new paths coexist during the phased rollout.

## Phased rollout (remaining)

Each phase swaps one workflow's viewer to OpenMS-Insight, then retires the
corresponding `src/render/*` usage. **Do not delete `src/render/*` until a
workflow's no-feature-loss audit passes against the contract above.**

1. **FLASHDeconv** — `content/FLASHDeconv/FLASHDeconvViewer.py` /
   `FLASHDeconvLayoutManager.py`: heatmaps, scan/mass `Table`s, deconv/annotated
   `LinePlot`s, `Scatter3D`, `DensityPlot`, `SequenceView` + `InternalFragmentMap`
   (when a sequence is provided), one shared `StateManager`. Parse layer →
   long format via the adapters above.
2. **FLASHTnT** — protein `Table` → `SequenceView` → tag `Table` →
   combined-spectrum `LinePlot` (tagger overlay) → `DensityPlot` → heatmaps.
   Preserve proteoform→scan resolution
   (`scan_resolution.py:build_proteoform_scan_map`): resolve `proteinIndex`→scan,
   expose a `scan`/`deconv_index` column so value-filters reproduce the PyArrow
   pushdown.
3. **FLASHQuant** — `content/FLASHQuant/FLASHQuantViewer.py`: single
   `FeatureView` (use `FeatureView.explode_traces` to convert the per-group
   array format).

### Layout parity (all phases)

Rebuild the configurable layout on native Streamlit:
`COMPONENT_OPTIONS`/`COMPONENT_NAMES` pickers drive which components instantiate;
render `[experiment][row][col]` with `st.columns` per row (≤3 cols); multi-
experiment side-by-side via top-level `st.columns(n)` (≤5), **each column its own
`StateManager` with a distinct `session_key`** so selections don't leak across
panels; preserve save/load of layouts
(`file_manager.get_results('layout','layout')`, `side_by_side`).

## Verification

- **OpenMS-Insight**: per-component unit tests (preprocess→cache→
  `_prepare_vue_data`) + `npm run build`. Current: **455 passed**, build green.
- **FLASHApp**: parse adapters unit-tested. Current: **53 passed**.
- **Per phase (requires a real dataset + a running app — not possible in CI
  without sample data):** drive each workflow's `COMPONENT_OPTIONS`, confirm
  every component renders and every cross-link in the plan's interaction
  checklist round-trips (scan→spectra/mass/3D/sequence, mass→3D isolation,
  protein→sequence/tag/spectrum, heatmap zoom levels, selection clear → None,
  fragment overlay, multi-experiment isolation, layout add/remove + side-by-side).

## Dependency

`openms-insight` is declared in `requirements.txt` (git dependency). For local
development use an editable/path install: `pip install -e ../OpenMS-Insight`,
and either rebuild `js-component/dist` after Vue changes or run the bundle in dev
mode with `SVC_DEV_MODE=true`.
