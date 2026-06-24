# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

**FLASHApp** is a Streamlit web application for visualizing **top-down proteomics** results from the OpenMS FLASH* tool family. It is built on the [OpenMS streamlit-template](https://github.com/OpenMS/streamlit-template) and bundles three independent sub-applications, each registered as a section in `app.py`:

- **FLASHDeconv** (⚡️) — spectral deconvolution: raw MS ion peaks → neutral monoisotopic masses (proteoforms), with isotope/charge resolution and FDR scoring.
- **FLASHTnT** (🧨) — tag-and-track top-down identification: runs FLASHDeconv, then matches short sequence "tags" against a protein FASTA database to identify proteins (PrSMs), with target/decoy FDR.
- **FLASHQuant** (📊) — proteoform quantification from FLASHDeconv mass traces (view-only; no run step).

The heavy lifting is done by **TOPP command-line tools** (`FLASHDeconv`, `FLASHTnT`, `DecoyDatabase`) shipped in the Docker image; the app drives them, parses their output into pandas DataFrames, caches them per workspace, and renders them through a custom **Vue.js Streamlit component** (`flash_viewer_grid`).

## Commands

```bash
# Run locally (online_deployment=false in settings.json → always "local" mode)
streamlit run app.py local

# Unit tests (pytest; uses fakeredis, needs pyopenms importable)
pytest tests/ -v
pytest tests/test_selection_clear.py -v                 # single file
pytest tests/test_selection_clear.py::test_name -v      # single test

# Lint (errors-only; mirrors .github/workflows/pylint.yml, which runs on `main`)
pylint $(git ls-files '*.py') --errors-only --disable=C0103,C0114,C0301,C0411,W0212,W0631,W0602,W1514,W2402,E0401,E1101,F0001,R1732

# Docker (full image with OpenMS + TOPP tools + Vue build)
docker build -f Dockerfile --no-cache -t flashapp:latest --build-arg GITHUB_TOKEN=<gh-token> .
docker run -p 8501:8501 flashapp:latest          # → http://localhost:8501
# Dockerfile.arm is the linux/arm64 variant (swaps miniforge installer to aarch64).
```

Python is pinned to **3.11** (matches the Docker runtime). `GITHUB_TOKEN` is required at build time to fetch the private `openms-streamlit-vue-component` submodule and OpenMS resources.

### Prerequisites before any production / Docker build

These are already the defaults on the `main`/release branches; verify them when building or debugging a blank viewer:

1. **Submodule present:** `git submodule init && git submodule update` (update to latest: `git submodule update --remote`).
2. **Vue component built & copied:** the bundle in `js-component/dist/` is produced from the `openms-streamlit-vue-component/` submodule (a Vite/Vue project). **Always prefer building the bundle via Docker, never a local Node toolchain** — the repo `Dockerfile` `js-build` stage (`node:21` → `npm install && npm run build`) is the canonical, reproducible build. To rebuild the committed bundle from *local* submodule source (e.g. after editing a `.vue` file), use a small Docker stage that `COPY`s the local submodule and runs `npm run build`, then export and copy `dist/` over `js-component/dist/`:
   ```dockerfile
   FROM node:21 AS build
   WORKDIR /openms-streamlit-vue-component
   COPY . .
   RUN npm install && npm run build
   FROM scratch AS export
   COPY --from=build /openms-streamlit-vue-component/dist /
   ```
   ```bash
   docker build -f vue-build.Dockerfile --target export \
     --output type=local,dest=./vue-dist openms-streamlit-vue-component
   # then replace js-component/dist/ with ./vue-dist/
   ```
   Only the prebuilt `js-component/dist/` is committed; the submodule source is fetched separately. (A bare local `npm install && npm run build` also works but is **not** preferred — Docker guarantees the toolchain.)
3. **`src/render/components.py` → `_RELEASE = True`** (loads the bundle from `js-component/dist/`). When `False`, the component is loaded from the Vite dev server at `http://localhost:5173` for live Vue development.
4. **`.streamlit/config.toml` → `developmentMode = false`**.

> Build order matters: the OpenMS/TOPP build must precede the Vue build in the Dockerfiles (see recent commits reordering this). `--no-cache` is recommended for the full image.

## Architecture

### The data pipeline (the core mental model)

Every sub-app follows the same path; understanding it requires reading `src/Workflow.py`, `src/parse/`, `src/workflow/FileManager.py`, and `src/render/`:

```
mzML upload ──► WorkflowManager.execution()
                   └─ executor.run_topp('FLASHDeconv' / 'FLASHTnT' / 'DecoyDatabase')
                        └─ writes out_deconv.mzML, annotated.mzML, *.tsv, *.feature, *.msalign
                            └─ src/parse/* turns those into pandas DataFrames
                                └─ FileManager.store_data(dataset_id, key, df)  ──► workspace cache
                                                                                       │
Viewer page ◄── render_grid() (src/render/render.py) ◄── reads cached dfs + layout ◄──┘
                   └─ get_component_function()  ──► Vue `flash_viewer_grid` component
```

- **Workflows:** `src/Workflow.py` defines `DeconvWorkflow`, `TagWorkflow` (FLASHTnT), and `QuantWorkflow`, all subclasses of `WorkflowManager` (`src/workflow/WorkflowManager.py`). Each implements `upload()`, `configure()`, `execution()`, `results()`. `execution()` runs TOPP tools, stores every output file via `FileManager.store_file()`, then calls the parsers. `TagWorkflow` chains `DecoyDatabase` → `FLASHDeconv` → `FLASHTnT`.
- **Parsers (`src/parse/`):** `deconv.py::parseDeconv`, `tnt.py::parseTnT`, `quant.py::parseQuant` / `flashquant.py::parseFLASHQuantOutput` are the entry points. The mzML→DataFrame heavy lifting lives in `masstable.py` (`parseFLASHDeconvOutput`, `parseFLASHTaggerOutput`) using a **multiprocessing pool**; `tag_resolution.py` maps tags ↔ proteoforms.
- **Renderers (`src/render/`):** `render.py` (`render_grid`, `render_component`) pushes state into the Vue component and reads user selections back out. `components.py` declares the component and defines the per-cell component classes: `PlotlyHeatmap`, `Tabulator` (Scan/Mass/Protein/Tag tables), `PlotlyLineplot`, `Plotly3Dplot`, `FDRPlotly`, `SequenceView`, `InternalFragmentMap`, `FLASHQuant`. `compression.py` compresses payloads sent to the browser; `StateTracker.py` tracks selection state across cells (cross-component linking via shared scan/mass identifiers).

### Pages and the Layout Manager

Each sub-app's pages (in `content/<SubApp>/`) are registered in `app.py`. The distinctive FLASHApp concept is the **Layout Manager** (`FLASHDeconvLayoutManager.py`, `FLASHTnTLayoutManager.py`): a grid editor where the user composes which visualization components appear in which cells (≤5 experiments, ≤3 columns/row). The chosen layout is persisted to the workspace cache and the **Viewer** page renders it (falling back to a built-in default if none is saved).

- Layout is stored via `FileManager.store_data` under dataset key **`'layout'`** for FLASHDeconv and **`'flashtnt_layout'`** for FLASHTnT (separate namespaces — they share the underlying `deconv_dfs`/`anno_dfs` data but keep independent layouts). It is JSON-importable/exportable.
- **Sequence Input** (`FLASHDeconvSequenceInput.py`) saves a proteoform sequence + fixed modifications to the `'sequence'` dataset; doing so unlocks the `Sequence view` and `Internal fragment map` components in the Layout Manager.
- **FLASHQuant** is simpler: File Upload + a single fixed-layout Viewer, no Layout Manager, and uses a separate cache subdirectory.

### Workspaces & FileManager

State lives in per-session **workspaces** (`enable_workspaces: true`, `workspaces_dir: ".."` → `../workspaces-FLASHApp/`). `FileManager` (`src/workflow/FileManager.py`) is the single gateway to the workspace's `cache/`: `store_file`, `store_data`, `get_results`, `result_exists`, `get_results_list`, `get_files`, keyed by a `dataset_id` (typically `<filename>_<timestamp>`). Demo workspaces are seeded from `example-data/workspaces/` (`demo_workspaces` in `settings.json`).

### Parameters

`ParameterManager` (`src/workflow/ParameterManager.py`) persists widget state to JSON and generates TOPP `.ini` files. `configure()` exposes TOPP tool parameters via `self.ui.input_TOPP('FLASHDeconv', exclude_parameters=[...], custom_defaults={...})`. **Widget keys must match keys in `default-parameters.json`.** `presets.json` holds named parameter bundles (`test_parameter_presets.py` guards this).

### Deployment & runtime (`entrypoint.sh`, `k8s/`)

The container entrypoint starts **Redis** + one or more **RQ workers** (queue `openms-workflows`) + Streamlit. When `STREAMLIT_SERVER_COUNT > 1`, it runs N Streamlit instances behind an **nginx** load balancer with sticky-cookie session routing. `QueueManager` (`src/workflow/QueueManager.py`) offloads `execution()` to RQ when `online_deployment` is set. The entrypoint is written to tolerate **Apptainer/Singularity** read-only rootfs (all runtime state goes under `$RUNTIME_DIR`, default `/tmp/opendiakiosk`). `k8s/` is a kustomize base + `overlays/prod` deploying to namespace `openms` as `ghcr.io/openms/flashapp:latest` behind nginx/traefik ingress; `clean-up-workspaces.py` runs via cron for periodic GC.

### CI (`.github/workflows/`)

- `build-and-test.yml` — multi-arch (amd64 + arm64) Docker build → merged manifest, kustomize/kubeconform lint, and **container smoke tests** under apptainer / nginx-on-kind / traefik-on-kind, then publishes images + an ORAS SIF to GHCR. (Its "test" jobs are deployment smoke tests, **not** pytest.)
- `unit-tests.yml` — the pytest suite. `pylint.yml` — lint. `build-windows-executable-app.yaml` + `test-win-exe-w-embed-py.yaml` — the PyInstaller desktop build. `ghcr-cleanup.yml` — registry GC.

## Conventions & gotchas

- **`app.py` sets multiprocessing start method to `spawn`** (polars + Unix fork are incompatible) and imports `pyopenms` early (required for the Windows build). Don't remove these.
- **Running workflows locally needs the TOPP binaries** (`FLASHDeconv`, `FLASHTnT`) on `PATH` — they ship only in the full Docker image. The upload/viewer/download paths work on pre-computed result files without them.
- **Windows packaged build:** `run_app.py` + `run_app_temp.spec` (PyInstaller). A `windows` arg in `sys.argv` triggers a working-directory `chdir` in `page_setup()`.
- Pages start with `page_setup()` from `src/common/common.py`, which initializes the workspace, sidebar, and params; call `save_params(params)` at the end. Use `show_fig()` / `show_table()` for consistent display.
- Decorate `configure()` / page sections with `@st.fragment` for partial reruns.
- Workflow display names map to lowercase-hyphenated keys ("FLASHTnT" → workflow dir / preset keys); TOPP params use colon paths (`tag:min_length`).
