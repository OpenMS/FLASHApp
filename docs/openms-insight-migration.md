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

## Phased rollout (status)

All three workflows render through the OpenMS-Insight engine (`src/render_oi/`)
**by default**. The `FLASHAPP_USE_OPENMS_INSIGHT` env flag is now an **opt-out**:
set it to `0`/`false`/`no`/`off` to fall back to the legacy `flash_viewer_grid`.
The layout managers are unchanged (they only emit component-name strings). Docker
images install OpenMS-Insight with its Vue bundle built from the pinned commit
(`openms-insight-build` stage in both Dockerfiles); see "Docker / packaging" below.

> **The default is on, but `src/render/*` is not yet deletable.** Everything up
> to the Vue render + click round-trip is verified on the bundled real data; the
> in-browser interaction is not (no headless browser in CI). The browser checklist
> below remains the gate for **removing** the legacy engine and the opt-out path —
> not for flipping the default, which the maintainer has chosen to enable now.

1. **FLASHDeconv** ✅ built (`src/render_oi/deconv_viewer.py`): heatmaps,
   scan/mass `Table`s, deconv/annotated `LinePlot`s, `Scatter3D` (optional
   `massIndex`), `DensityPlot`, `SequenceView` + `InternalFragmentMap` (when a
   sequence is provided), per-experiment `StateManager`. Parse → long format.
2. **FLASHTnT** ✅ built (`src/render_oi/tnt_viewer.py`): protein `Table` →
   `SequenceView` (coverage) → tag `Table` → combined-spectrum `LinePlot`
   (tagger overlay) → `DensityPlot` → heatmaps. Proteoform→scan resolution
   (`build_proteoform_scan_map`) resolves `proteinIndex`→`deconvIndex` before
   downstream panels render.
3. **FLASHQuant** ✅ built (`src/render_oi/quant_viewer.py`): feature-group
   `Table` + single `FeatureView`; `explode_quant_traces_long` converts the
   per-group array format to long trace-point format.

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
  `_prepare_vue_data`) + `npm run build`. Current: **464 passed**, build green.
- **FLASHApp**: parse adapters + per-workflow viewer engines tested against the
  **bundled real example workspaces** (`tests/test_{deconv,tnt,quant}_viewer_realdata.py`,
  `tests/test_long_format.py`). Current: **72 passed, 2 skipped**.
- **Browser audit (the deletion gate):** see the no-feature-loss audit and
  click-by-click checklist below. This is the one layer CI cannot cover.

## Docker / packaging

OpenMS-Insight's wheel is built by **hatchling**, which `force-include`s the
pre-built Vue bundle at `openms_insight/js-component/dist` **only if it exists on
disk** — and that bundle is gitignored. So a plain `pip install
git+https://…/OpenMS-Insight` (what the `requirements.txt` line would do) yields a
package with **no frontend**. Three install paths handle this correctly:

- **Docker** (`Dockerfile`, `Dockerfile.arm`): a dedicated `openms-insight-build`
  node stage clones the repo at the pinned commit (`ARG OPENMS_INSIGHT_REF`,
  default = the validated SHA), runs `npm run build`, mirrors `js-component/dist`
  → `openms_insight/js-component/dist`, and hands the populated checkout to the
  python stage, which `pip install`s it. The `requirements.txt` `openms-insight`
  line is stripped before `pip install -r` (alongside `pyopenms`) so it installs
  exactly once, with a working bundle. Bump the image by overriding
  `OPENMS_INSIGHT_REF` (or rely on the branch cache-bust `ADD`).
- **Local dev / CI / web sessions**: the SessionStart hook
  (`.claude/hooks/session-start.sh`) does the same (strip line → `npm run build` →
  mirror dist → `pip install -e ../OpenMS-Insight`). Or run the Vue bundle in dev
  mode with `SVC_DEV_MODE=true`.
- **`requirements.txt`**: keeps the `git+…@<branch>` ref as a declarative pointer
  for the migration branch; both install paths above strip and replace it, so it
  is never the thing that actually provides the frontend.

---

# No-feature-loss audit (deletion gate)

This is the per-workflow contract that must hold before `src/render/*` and the
`js-component/` submodule are deleted (Verification §4 of the plan). For each
workflow: (A) every legacy `StateTracker` key and `update.py` filter is mapped to
its `src/render_oi/*` replacement, (B) every component's columns/axes are
accounted for, and (C) a click-by-click browser checklist confirms the
interaction round-trips. Run the app with the new engine:

```bash
# The OpenMS-Insight engine is the default; this just makes it explicit.
FLASHAPP_USE_OPENMS_INSIGHT=1 streamlit run app.py local
# To audit the legacy engine instead: FLASHAPP_USE_OPENMS_INSIGHT=0 streamlit run app.py local
```

Legend: ✅ data-path verified on real data (test name in parentheses) · 👁 needs
browser confirmation.

## Workflow 1 — FLASHDeconv

Engine: `src/render_oi/deconv_viewer.py`. Tests:
`tests/test_deconv_viewer_realdata.py`.

### State keys (legacy `StateTracker`/`update.py` → new)

| legacy key | legacy filter (`update.py`) | new identifier→column | status |
|---|---|---|---|
| `scanIndex` | `per_scan_data.iloc[scanIndex]` (L131‑134) | scan `Table` `interactivity={scanIndex:index}`; spectra/mass/3D/seq `filters={scanIndex:index}` | ✅ row counts exact (`test_scan_click_cross_link_row_counts`) |
| `massIndex` | `SignalPeaks[massIndex]`/`NoisyPeaks[massIndex]` (L142‑146) | mass `Table` `interactivity={massIndex:mass_id}`; `Scatter3D` `optional_filters={massIndex:mass_id}` | ✅ 3D isolates mass 0 (same test) |
| `heatmap_deconv/_deconv2/_raw/_raw2` | four bespoke `render_heatmap` range keys (L149‑176) | one `zoom_identifier` per `Heatmap` (`{comp}_zoom`) | ✅ levels build; 👁 zoom round-trip |
| `sequenceOut` / sequence cache | `get_sequence()` (L13‑29) | `SequenceView(sequence_data=<string>)` when `result_exists('sequence','sequence')` | 👁 |

### Component coverage

| component | columns / axes preserved | status |
|---|---|---|
| Deconv/Raw MS1/MS2 `Heatmap` | x=`rt`, y=`mass`, color=`intensity` (log), zoom compression | ✅ build (incl. 608K→200K→20K levels); 👁 render |
| Scan `Table` | `index, Scan, MSLevel, RT, PrecursorMass, #Masses` | ✅ |
| Mass `Table` | exploded `MonoMass, SumIntensity, Min/MaxCharges, Min/MaxIsotopes, CosineScore, SNR, QScore` + `mass_id` | ✅ (`_explode_mass_table`) |
| Deconv `LinePlot` | x=`mass`, y=`intensity` | ✅ peak count exact |
| Annotated `LinePlot` | annotated peaks from `combined_spectrum` | ✅ |
| `Scatter3D` (Precursor S/N) | x=`mz*charge`, y=`charge`, z=`intensity`, signal/noise color | ✅ |
| `DensityPlot` (FDR) | precomputed `density_target`/`density_decoy` curves | ✅ (200‑pt target, empty decoy handled) |
| `SequenceView` (+`InternalFragmentMap`) | residues + fixed mods (C/M), deconvolved matching | 👁 (no sequence in bundled FD set) |

### 👁 Browser checklist (FLASHDeconv)

1. All `COMPONENT_OPTIONS` from `FLASHDeconvLayoutManager` render without error.
2. Click a scan row → deconv spectrum, annotated spectrum, mass table, **and**
   3D plot all update to that scan.
3. Click a mass row → 3D plot isolates that mass's signal/noisy peaks.
4. Zoom a heatmap → point density increases (correct compression level); zoom
   out → returns to overview.
5. Deselect (click empty) → dependent panels clear (selection round-trips to
   `None`), no stale data.
6. Submit a sequence → Sequence View + Internal Fragment Map appear and annotate
   the selected scan's peaks.
7. Layout editor: add/remove a component, save, reload → viewer reflects it.
8. Side-by-side (2 experiments): a selection in panel A does **not** move
   panel B (distinct `session_key` per panel — `oi_state_{panel_key}`).

## Workflow 2 — FLASHTnT

Engine: `src/render_oi/tnt_viewer.py`. Tests:
`tests/test_tnt_viewer_realdata.py`.

### State keys

| legacy behavior (`update.py`) | new | status |
|---|---|---|
| `proteinIndex` → `proteoform_scan_map` → `deconv_index`, pushdown `field('index')==deconv_index` (L122‑134) | protein `Table` `interactivity={proteinIndex:index}`; `render_experiment_tnt` resolves `proteinIndex`→`deconvIndex`; spectra/mass `filters={deconvIndex:index}` | ✅ resolution + tag count exact (`test_protein_click_cross_links`) |
| tag table: pushdown `field('Scan')==scan`, stamp `ProteinIndex` (L177‑192) | tag `Table` `filters={proteinIndex:ProteinIndex}` | ✅ 15 tags for protein 0 |
| sequence: `load_entry(sequence_data_ds, proteinIndex)` (L194‑213) | `SequenceView` `filters={proteinIndex:proteoform_index}`, per-proteoform `coverage[]`/`maxCoverage` | ✅ 450 residues + coverage |

### Component coverage

| component | preserved | status |
|---|---|---|
| Protein `Table` | `index, accession, description, ProteoformMass, Coverage(%), TagCount, ProteoformLevelQvalue` | ✅ |
| Tag `Table` | `TagIndex, TagSequence, StartPos, EndPos, Length, Score, DeltaMass` | ✅ |
| Combined `LinePlot` (tagger) | deconv primary (564) + annotated overlay (65 657) | ✅ overlay (`test_combined_spectrum_overlay`) |
| `SequenceView` | coverage shading, fixed mods, `settings.{tolerance,ion_types}` | ✅ data; 👁 render |
| `DensityPlot` (id-FDR) | `density_id_target/decoy` (empty in antibody set) | ✅ both-empty handled |
| MS1 raw/deconv `Heatmap` | as Deconv | ✅ |

### 👁 Browser checklist (FLASHTnT)

1. Click a protein row → sequence view (coverage colored), tag table, **and**
   combined spectrum all update to that proteoform.
2. Tag table shows only the selected protein's tags; on-spectrum tag overlay (if
   enabled) matches.
3. Combined spectrum shows both series (deconv sticks + annotated overlay).
4. Sequence view: per-residue coverage shading + fixed-mod (C/M) styling +
   correct `ion_types`/`tolerance` from settings.
5. Selection clear → dependent panels reset.
6. Side-by-side isolation (distinct `oi_tnt_state_{panel_key}`).
7. **Tagger tag-geometry** (charge buttons + inter-residue amino-acid arrows
   driven by `tagIndex`/`tagData`): confirm parity with legacy or log as the one
   known gap (see footnote ¹ in the component table) before deletion.

## Workflow 3 — FLASHQuant

Engine: `src/render_oi/quant_viewer.py`. Tests:
`tests/test_quant_viewer_realdata.py`. No cross-component linking, no layout
manager (fixed single view).

| component | preserved | status |
|---|---|---|
| Feature-group `Table` | scalar summary columns; click sets `featureGroup` | ✅ 1437 rows |
| `FeatureView` | per-charge 3D traces (x=`mz`, y=`rt`, z=`intensity`), isotope hover; `filters={featureGroup:feature_group}` | ✅ exact per-group point counts (0/5/100 → 1384/239/208) |

### 👁 Browser checklist (FLASHQuant)

1. Feature-group table renders all groups; selecting a group updates the
   FeatureView to that group's traces.
2. Traces colored per charge; isotope info in hover; conflict-resolution
   highlighting if a `conflict_resolution_dfs` cache is present.

## After all three checklists pass

1. Remove the `FLASHAPP_USE_OPENMS_INSIGHT` opt-out (it already defaults on),
   making the new engine unconditional in the three viewers.
2. Delete `src/render/{components,render,update,StateTracker,initialize,
   compression}.py` and the `js-component/` bundle + `openms-streamlit-vue-component`
   submodule. Keep `src/render/{sequence,sequence_data_store,scan_resolution}.py`
   if still imported by the parse layer (they are — `scan_resolution` is used by
   `tnt_viewer`).
3. Drop the now-dead `quant_visualization` / `flash_viewer_grid` code paths.
4. Re-run the full suite + `npm run build`; commit the retirement per workflow.
