# FLASHApp

FLASHApp's own guidance, imported by the template-owned `CLAUDE.md` (which describes the framework).

## What This Is

**FLASHApp** is a Streamlit web application for **top-down proteomics** with the OpenMS FLASH* tool family. It is built on the [OpenMS streamlit-template](https://github.com/OpenMS/streamlit-template). `app.py` registers one analysis section (⚡️ **Top-down analysis**) and the developer-only FLASHIda section:

- **Top-down analysis** (⚡️, one workflow, `content/TopDown/`) — spectral deconvolution of every run (FLASHDeconv: raw MS ion peaks → neutral monoisotopic masses, with isotope/charge resolution and FDR scoring) and, with **Identify proteoforms** on and a protein database, tag-and-track identification (FLASHTnT matches short sequence "tags" against the database to identify proteoforms (PrSMs), with target/decoy FDR; FLASHAnnotator annotates their fragments). Off, a run is deconvolved only (MS1-only or intact-mass runs, samples without a fitting database).
- **FLASHIda** (🧠) — developer feature of the local app only: registered in `app.py` when `src/flashida.py::available()` (local app + developer mode unlocked; the unlock is remembered in `<local workspaces>/.developer-mode` until Lock). Acquisition (`IdaWorkflow` starts Flash.exe with `<method>.json` for every new `*FLASHIda_<method>.raw`; it remembers the method's `runtime.log_dir`) and a live Viewer (`src/render/ida_live.py`, an `st.fragment(run_every=…)`) that re-reads a run folder's `scan_commands.tsv`, `scan_results.tsv`, `identification.tsv` and `pooled_identification.tsv` (`src/parse/ida_log.py`) and shows the survey masses with the chosen precursors, the precursors and the scans commanded for each.

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
docker build -f Dockerfile --no-cache -t flashapp:latest --secret id=github_token,env=GITHUB_TOKEN .
docker run -p 8501:8501 flashapp:latest          # → http://localhost:8501
# Dockerfile.arm is the linux/arm64 variant (swaps miniforge installer to aarch64).
```

Python is pinned to **3.11** (matches the Docker runtime). `GITHUB_TOKEN` is required at build time to fetch OpenMS resources.

### openms-insight

The viewers need `openms_insight` importable. Until a release after 0.4.0 includes t0mdavid-m/openms-insight#21 (MirrorPlot shared selections, SeriesPlot `y_range_mode`/`hover_column`, NaN-safe table formatters, used by the validation, QC and run comparison views) #23 (per-scan fragmentation method in the sequence view) and #34 (Heatmap `hover_column`, `highlight` ring and `on_select`, used by the FLASHTnT proteoform heatmap), merged into the branch `claude/project-thread-hxjbbs`, it is **built from source**: the Dockerfiles' `insight-build` stage (`INSIGHT_REF` build arg) and `build-desktop-app.yaml`'s `build-openms-insight` job (a wheel for the Windows, macOS and Linux apps) clone `t0mdavid-m/openms-insight`, build its Vue bundle and install it with `--no-deps` (its dependencies come from `requirements.txt`). For local work, see README → Visualizations. When the release is out, add `openms-insight>=<version>` to `requirements.txt` and drop those build steps. The unit tests do not import it.

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

- **Workflows:** `src/Workflow.py` defines `TopDownWorkflow` (the former FLASHDeconv and FLASHTnT workflows, merged; `TagWorkflow`/`DeconvWorkflow` are aliases for queued jobs of older builds) and `QuantWorkflow`, all subclasses of `WorkflowManager` (`src/workflow/WorkflowManager.py`). Each implements `upload()`, `configure()`, `execution()`, `results()`. `execution()` runs TOPP tools, stores every output file via `FileManager.store_file()`, then calls the parsers. `TopDownWorkflow` runs `FLASHDeconv` per run and, with `identify` on, `DecoyDatabase` → `FLASHTnT` → `FLASHAnnotator`. It keeps the FLASHTnT workflow's name and folder (`flashtnt/`, its params and the `flashtnt` presets); the MS data of the former `flashdeconv/` workflow is listed by reference (`external_files.txt`) on first use (`adopt_legacy_inputs`).
- **Parsers (`src/parse/`):** `deconv.py::parseDeconv`, `tnt.py::parseTnT` are the entry points. The mzML→DataFrame heavy lifting lives in `masstable.py` (`parseFLASHDeconvOutput`, `parseFLASHTaggerOutput`) using a **multiprocessing pool**; `tag_resolution.py` maps tags ↔ proteoforms.
- **Renderers (`src/render/`):** `insight_data.py` turns the cached results into the lazy frames each component is built from (one row per mass, signal peaks per mass, best hit per spectrum, …; Streamlit-free, tested in `tests/test_insight_data.py`). `render.py::render_grid` builds one openms-insight component per layout cell; the components of a grid share a `StateManager` (identifiers `scan`, `mass`, `protein`, `pscan`, `tag`, `aa`, `sequence`, `feature_group`), cleared when the grid switches dataset. Each component's preprocessed cache lives in `<cache>/insight/<dataset>/`; its cache id hashes the source results (size, mtime) and the code of `render.py`/`insight_data.py`, so editing either rebuilds the caches. `compression.py` is still used by the parsers.

### Pages and the Layout Manager

Each sub-app's pages (in `content/<SubApp>/`, `content/TopDown/` for the merged workflow) are registered in `app.py`. The distinctive FLASHApp concept is the **Layout Manager**, opened as a dialog by **Edit layout** on each Viewer page (`src/render/layout_editor/dialogs.py`; the editor itself is `src/render/layout_editor/frontend/index.html`, plain JS with a searchable palette of collapsible groups from `catalog.py`): a grid editor where the user composes which visualization components appear in which cells (≤5 experiments, ≤3 columns/row). The chosen layout is persisted to the workspace cache and the **Viewer** page renders it (falling back to a built-in default if none is saved).

- One Viewer (`TopDownViewer.py`) shows every run in one combined view (`tool='topdown'` in `render_grid`): the catalog is `TOPDOWN_COMPONENTS`, both kinds together. Identification tables set the `scan` selection too (beside `pscan`), so spectra, the mass table and the tag table follow either the scan table or an identification; the sequence view shows the proteoform of an identified run, else the sequence set in the view; identification components on a run deconvolved only say so. One layout, **`'topdown_layout'`** (`src/render/topdown_layout.py`). Until one is saved, each experiment starts from the default of its run: the layout a former Viewer saved for its kind (**`'layout'`** deconvolution, **`'flashtnt_layout'`** identification), else the deconvolution default, or the Native/Lysate default for an identified run by its `-native` setting (`sample_type.of_result`). It is JSON-importable/exportable.
- **The sequence of a run that was not identified** is set in the sequence view's own **Change sequence** dialog (there is no separate page). `render.py::_deconv_sequence` saves it with its fixed modifications to the `'sequence'` dataset; a saved sequence unlocks the `Internal fragment map` component in the Layout Manager. Without one the view starts empty.
- **The analysis choice** (`src/analysis.py`, Streamlit-free): `identify` (on/off), the sample type (Native/Lysate, the `flashtnt` presets, applied by the **Sample** control; `sample_type.is_on` reads the flag, whose checkbox holds the ini text `'true'`/`'false'`) and `quantification` (`none` | `label-free` | `isobaric`, falling back to the former `compare-runs` switch) with `isobaric-plex`. Six analyses result. `TopDownWorkflow.analysis_section()` draws them in the **Your analysis** box above the parameters (the page calls it before `show_parameter_section()`); `configure()` only draws the tabs that follow from them. The Native/Lysate preset buttons are hidden (`_no_preset_buttons`, module-level: the workflow object is pickled for its run).
- **Isobaric labels** (`src/quant/isobaric.py`, Streamlit-free): FLASHDeconv runs with `-iq:type <plex>` and `-out_quant`; `read_quant` reads the reporter table (stored per run as `isobaric_dfs`, with `isobaric_plex`), the channel table (channel, sample, condition) lives in `<workflow_dir>/isobaric-channels.tsv`, and `build` turns all runs of a start (fractions of one labelled set) into the comparison contract per level (feature: spectra within 10 ppm; proteoform/protein: FLASHTnT PrSMs at q ≤ 0.01; a native PrSM without any mass is a feature of its own). `src/quant/run.py::compare_isobaric` stores it with `kind: 'isobaric'` in the manifest, so `QuantSession.isobaric` makes the comparison speak of samples instead of runs. The Viewer adds `reporter_table`/`reporter_ions` (`ISOBARIC_COMPONENTS`) to the default of a labelled run (`topdown_layout.with_reporters`).
- **User Guide** (`content/user_guide.py` + `.md`): a line `@video name.mp4 | caption` embeds `static/videos/name.mp4`. The four instruction videos are Playwright recordings of the app with caption overlays.
- **Workflow pages** follow OpenDDA's three steps: `<Tool>Upload.py` (Upload data, with the manual result upload as a second tab), `<Tool>Configure.py` (Describe experiment; `configure()` draws numbered tabs, the first being **Samples and groups**) and `<Tool>Run.py` (Run analysis). Each opens with `page_intro` and ends with `next_steps` from `src/common/guide.py`.
- **Comparing runs** (`src/quant/`, `src/render/quant_view.py`, `src/render/quant_cells.py`): with **Quantification: Label-free** on, FLASHDeconv runs with `-quant` and, after the per-file loop, `src/quant/run.py` runs FLASHQuantConsensus over the stored runs (FLASHTnT adds the PrSMs as `-in_id` and, when the tool has it, `-out_id`, the proteoform of each consensus feature). `src/quant/contract.py` reads the output into the contract shared with the quant unification work (a long table per level: feature, proteoform, protein); `src/quant/session.py` stores it as a dataset `comparison_<time>` in the same cache (`quant_*` names, `tool` in `quant_manifest`) and runs the comparison chain of `src/parse/quant_compare.py` per level, cached by a stamp of the design and settings. The Viewer lists it as **Comparison: N runs (date)**; selecting it replaces the grid with the comparison grid (`QUANT_COMPONENTS`, its own layout under `quant_layout`; `quant_view.migrate_layout` moves the former `flashtnt_quant_layout` there) and adds **Edit comparison** (conditions, level, analysis settings). The conditions start from the workflow's Samples table (`experimental-design.tsv`).

### Workspaces & FileManager

State lives in per-session **workspaces** (`enable_workspaces: true`, `workspaces_dir: ".."` → `../workspaces-FLASHApp/`). `FileManager` (`src/workflow/FileManager.py`) is the single gateway to the workspace's `cache/`: `store_file`, `store_data`, `get_results`, `result_exists`, `get_results_list`, `get_files`, keyed by a `dataset_id` (typically `<filename>_<timestamp>`). Demo workspaces are seeded from `example-data/workspaces/` (`demo_workspaces` in `settings.json`).

### Parameters

`ParameterManager` (`src/workflow/ParameterManager.py`) persists widget state to JSON and generates TOPP `.ini` files. `configure()` exposes TOPP tool parameters via `self.ui.input_TOPP('FLASHDeconv', exclude_parameters=[...], custom_defaults={...})`. **Widget keys must match keys in `default-parameters.json`.** `presets.json` holds named parameter bundles (`test_parameter_presets.py` guards this).

### Deployment & runtime (`entrypoint.sh`, `k8s/`)

The container entrypoint starts **Redis** + one or more **RQ workers** (queue `openms-workflows`) + Streamlit. When `STREAMLIT_SERVER_COUNT > 1`, it runs N Streamlit instances behind an **nginx** load balancer with sticky-cookie session routing. `QueueManager` (`src/workflow/QueueManager.py`) offloads `execution()` to RQ when `online_deployment` is set. The entrypoint is written to tolerate **Apptainer/Singularity** read-only rootfs (all runtime state goes under `$RUNTIME_DIR`, default `/tmp/opendiakiosk`). `k8s/` is a kustomize base + `overlays/prod` deploying to namespace `openms` as `ghcr.io/openms/flashapp:latest` behind nginx/traefik ingress; `clean-up-workspaces.py` runs via cron for periodic GC.

### CI (`.github/workflows/`)

- `build-and-test.yml` — multi-arch (amd64 + arm64) Docker build → merged manifest, kustomize/kubeconform lint, and **container smoke tests** under apptainer / nginx-on-kind / traefik-on-kind, then publishes images + an ORAS SIF to GHCR. (Its "test" jobs are deployment smoke tests, **not** pytest.)
- `unit-tests.yml` — the pytest suite. `pylint.yml` — lint. `build-desktop-app.yaml` — one run builds the Electron app for all three platforms, from the template's unified build. The `openms` job (`openms.yml`) builds OpenMS `t0mdavid-m/OpenMS@FDmerge` (it has FLASHAnnotator) from source on Windows, macOS and Linux, once per commit (cached by SHA; a released version such as `3.5.0` would be downloaded instead). One `build` job, a matrix over Windows / macOS (arm64) / Linux, then runs the same bash steps everywhere: `electron/scripts/stage-openms.sh` (the built package as `OPENMS_PACKAGE`), `stage-python.sh` plus the openms-insight wheel and `polars-lts-cpu`, copy the app, `package.sh` (NSIS / signed and notarized dmg (when the `MACOS_CERTIFICATE`/`APPLE_*` secrets exist) / AppImage, named `<INSTALLER_NAME>-Setup.exe`, `-macOS-arm64.dmg`, `-Linux-x86_64.AppImage` while the app is `APP_NAME` FLASHApp), then `install-app.sh`, `smoke-payload.js`, `smoke-app.js` and `e2e.js --scenario tests/e2e/electron-workflows.js`. On a release every installer is attached as an asset, then `publish-web-image` runs `build-and-test.yml`. The Quickstart's download buttons (`src/common/installers.py`) link the installers of the latest release. `ghcr-cleanup.yml` — registry GC.
- `template-sync.yml` — weekly merge of the streamlit-template (`template` in `app.yaml`); conflicts in files listed under `sync.app_owned` keep FLASHApp's side. Merge those PRs with a merge commit, never squash, or the next sync conflicts again. `auto-merge.yml` and `.github/CODEOWNERS` form the review queue (`docs/review-queue.md`).

## Conventions & gotchas

- **`app.py` sets multiprocessing start method to `spawn`** (polars + Unix fork are incompatible) and imports `pyopenms` early (required for the Windows build). Don't remove these.
- **Vendor input:** FLASHDeconv and FLASHTnT accept Thermo `.raw` and Bruker timsTOF `.d` (`MS_DATA_TYPES` in `src/Workflow.py`). `execution()` converts them with `self.executor.convert_to_mzml()` (from the streamlit-template): FileConverter + ThermoRawFileParser for `.raw` (mono in the Docker images, .NET on Windows, `FileConverter` in the Windows `TOPP_TOOLS`), `src/python-tools/bruker_to_mzml.py` (timsrust + pyOpenMS) for `.d`. Converted mzML is cached in `<workflow_dir>/converted-mzML/`. Online, a `.d` is uploaded zipped.
- **Import from a link** (`src/data_import/`, from streamlit-template #428): the upload page gets a section that imports files from a link (HTTPS incl. presigned, public S3/Swift folders such as CSC Allas, OneDrive/SharePoint, Zenodo, Figshare; tar/zip listed in place) into the matching `upload_widget()` folder, as a background process reporting through `<workflow_dir>/link-import/`. Online, an SSRF guard refuses non-public addresses. Toggled by `data_import` in `settings.json`.
- **Running workflows locally needs the TOPP binaries** (`FLASHDeconv`, `FLASHTnT`) on `PATH` — they ship only in the full Docker image. The upload/viewer/download paths work on pre-computed result files without them.
- **Windows packaged build:** `run_app.py` + `run_app_temp.spec` (PyInstaller). A `windows` arg in `sys.argv` triggers a working-directory `chdir` in `page_setup()`.
- Pages start with `page_setup()` from `src/common/common.py`, which initializes the workspace, sidebar, and params; call `save_params(params)` at the end. Use `show_fig()` / `show_table()` for consistent display.
- **Light and dark theme**: `.streamlit/config.toml` defines `[theme.light]` and `[theme.dark]` without `theme.base`, so the app follows the OS color scheme (switch under ⋮ → Settings). `show_fig()` runs figures through `follow_app_theme()`; call it yourself before a direct `st.plotly_chart`. The Vue component reads the active theme from Streamlit's render args. Logos: `*-light.svg` (white tile, dark text) for the light theme, `flashapp-icon.svg` / `*-white.svg` (dark tile, white text) for the dark theme, picked by `st.context.theme.type`.
- Decorate `configure()` / page sections with `@st.fragment` for partial reruns.
- Workflow display names map to lowercase-hyphenated keys ("FLASHTnT" → workflow dir / preset keys); TOPP params use colon paths (`tag:min_length`).
