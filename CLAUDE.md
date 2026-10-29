# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

**FLASHApp** is a Streamlit web application for visualizing **top-down proteomics** results from the OpenMS FLASH* tool family. It is built on the [OpenMS streamlit-template](https://github.com/OpenMS/streamlit-template) and bundles three independent sub-applications, each registered as a section in `app.py`:

- **FLASHDeconv** (⚡️) — spectral deconvolution: raw MS ion peaks → neutral monoisotopic masses (proteoforms), with isotope/charge resolution and FDR scoring.
- **FLASHTnT** (🧨) — tag-and-track top-down identification: runs FLASHDeconv, then matches short sequence "tags" against a protein FASTA database to identify proteins (PrSMs), with target/decoy FDR.
- **FLASHQuant** (📊) — proteoform quantification from FLASHDeconv mass traces (view-only; no run step).

The heavy lifting is done by **TOPP command-line tools** (`FLASHDeconv`, `FLASHTnT`, `DecoyDatabase`) shipped in the Docker image; the app drives them, parses their output into pandas DataFrames, caches them per workspace, and renders them with **openms-insight** components (`Table`, `Heatmap`, `DeconvolvedSpectrum`, `Plot3D`, `SeriesPlot`, `ProteinSequenceView`).

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

# Docker (full image with OpenMS + TOPP tools + openms-insight built from source)
docker build -f Dockerfile --no-cache -t flashapp:latest --build-arg GITHUB_TOKEN=<gh-token> .
docker run -p 8501:8501 flashapp:latest          # → http://localhost:8501
# Dockerfile.arm is the linux/arm64 variant (swaps miniforge installer to aarch64).
```

Python is pinned to **3.11** (matches the Docker runtime). `GITHUB_TOKEN` is required at build time to fetch OpenMS resources.

### openms-insight

The viewers need `openms_insight` importable. Until the release with the components the viewers use (`DeconvolvedSpectrum`, `Plot3D`, `SeriesPlot`, `ProteinSequenceView`, the Table options they use) is on PyPI, it is **built from source**: the Dockerfiles' `insight-build` stage (`INSIGHT_REF` build arg) and the Windows workflow's `build-openms-insight` job clone `t0mdavid-m/openms-insight`, build its Vue bundle and install it with `--no-deps` (its dependencies come from `requirements.txt`). For local work, see README → Visualizations. When the release is out, add `openms-insight>=<version>` to `requirements.txt` and drop those build steps. The unit tests do not import it.

`.streamlit/config.toml` → `developmentMode = false` for production builds.

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
Viewer page ◄── render_grid() (src/render/render.py) ◄── reads cached results + layout ◄──┘
                   └─ one openms-insight component per layout cell (shared StateManager)
```

- **Workflows:** `src/Workflow.py` defines `DeconvWorkflow`, `TagWorkflow` (FLASHTnT), and `QuantWorkflow`, all subclasses of `WorkflowManager` (`src/workflow/WorkflowManager.py`). Each implements `upload()`, `configure()`, `execution()`, `results()`. `execution()` runs TOPP tools, stores every output file via `FileManager.store_file()`, then calls the parsers. `TagWorkflow` chains `DecoyDatabase` → `FLASHDeconv` → `FLASHTnT`.
- **Parsers (`src/parse/`):** `deconv.py::parseDeconv`, `tnt.py::parseTnT`, `quant.py::parseQuant` / `flashquant.py::parseFLASHQuantOutput` are the entry points. The mzML→DataFrame heavy lifting lives in `masstable.py` (`parseFLASHDeconvOutput`, `parseFLASHTaggerOutput`) using a **multiprocessing pool**; `tag_resolution.py` maps tags ↔ proteoforms.
- **Renderers (`src/render/`):** `insight_data.py` turns the cached results into the lazy frames each component is built from (one row per mass, signal peaks per mass, best hit per spectrum, …; Streamlit-free, tested in `tests/test_insight_data.py`). `render.py::render_grid` builds one openms-insight component per layout cell; the components of a grid share a `StateManager` (identifiers `scan`, `mass`, `protein`, `pscan`, `tag`, `aa`, `sequence`, `feature_group`), cleared when the grid switches dataset. Each component's preprocessed cache lives in `<cache>/insight/<dataset>/`; its cache id hashes the source results (size, mtime) and the code of `render.py`/`insight_data.py`, so editing either rebuilds the caches. `compression.py` is still used by the parsers.

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
