# Random "white screen" in FLASHApp: root cause analysis

*Investigated 2026-09-10 on the local build (Windows 11, Edge driven by Playwright, Streamlit 1.49.1, polars 1.44.2, pyarrow 19.0.1, Python 3.11).*

## Summary

The "white screen that needs a refresh" is the **Streamlit server process dying with a native access violation**, not a front-end glitch. It was reproduced three times in the FLASHDeconv Viewer, each time within about a second of an action that makes the grid re-render (selecting an experiment, navigating back to the Viewer, clicking a Scan Table row).

Causal chain:

1. Every FLASH viewer grid cell (Vue component `flash_viewer_grid`) calls `Streamlit.setComponentValue` on mount and on every echoed render (`App.vue` watches the selection store with `{ deep: true, immediate: true }`; `updateRenderData` copies the `selection_store` that Python sends back into that store).
2. Streamlit's `WidgetStateManager.setJsonValue` has no equality check, so every echo becomes a full `rerun_script` request. Six cells produce six requests per render.
3. With the default `runner.fastReruns=true`, `AppSession.request_rerun` calls `request_stop()` on the running ScriptRunner and immediately starts a new script thread. The old thread only notices at its next Streamlit call; inside `initialize_data` / `filter_data` it spends long stretches in polars and pyarrow without one, so it keeps executing. The server log shows four "Beginning script thread" lines within 330 ms.
4. The concurrent threads share the same `st.session_state.plot_data` LazyFrames / pyarrow datasets and the same `st.cache_data` entries (polars DataFrames stored as pickles and unpickled on every read). polars' native deserializer, invoked by `st.cache_data` from `render_heatmap` (`src/render/update.py:158`), crashes with an access violation while a second thread of the same page is in `filter_data` (`update.py:141`, a pyarrow `to_table` filter). No Python exception is possible; the whole process exits.

User-visible result: the status bar shows CONNECTING, then the "Connection error - Streamlit server is not responding" dialog. If the crash lands right after a sidebar click the main area has already been cleared for the new page, so the page is empty and white. Refreshing starts a fresh session (run status, selections and layout state are gone). On Kubernetes the Streamlit process is PID 1 of its container, so a crash restarts the pod and takes down any workflow that fell back to local execution; the RQ path survives.

## Evidence

Faulthandler output from crash #2 (10:18:21):

```
Windows fatal exception: access violation

Current thread 0x0000adf8 (most recent call first):
  File "polars\dataframe\frame.py", line 562 in deserialize
  File "polars\dataframe\frame.py", line 1230 in __setstate__
  File "streamlit\runtime\caching\cache_data_api.py", line 651 in read_result
  File "streamlit\runtime\caching\cache_utils.py", line 251 in _get_or_create_cached_value
  File "streamlit\runtime\caching\cache_utils.py", line 227 in __call__
  File "src\render\update.py", line 158 in filter_data            # render_heatmap(...)
  File "src\render\render.py", line 38 in render_component
  File "src\render\render.py", line 105 in render_grid
  File "content\FLASHDeconv\FLASHDeconvViewer.py", line 117 in <module>
  File "streamlit\runtime\scriptrunner\script_runner.py", line 685 in _run_script

Thread 0x0000ad04 (most recent call first):                          # second runner, same page, still executing
  File "src\render\update.py", line 141 in filter_data            # pyarrow to_table(filter=index == -1)
  File "src\render\render.py", line 38 in render_component
  File "src\render\render.py", line 105 in render_grid
  File "content\FLASHDeconv\FLASHDeconvViewer.py", line 117 in <module>
  File "streamlit\runtime\scriptrunner\script_runner.py", line 685 in _run_script
```

Server timeline before crash #2:

| Time | Event | Meaning |
|---|---|---|
| 10:18:20.519 | Running script (Viewer) | sidebar click on Viewer |
| 10:18:21.005 | `rerun_script` with component value `{scanIndex:0, massIndex:0, counter:2}` | first grid cell echoes its state |
| 10:18:21.006 | Beginning script thread | runner #2 starts while runner #1 still executes |
| 10:18:21.286 | `rerun_script` (two component values) | more cells echo |
| 10:18:21.288 | Beginning script thread | runner #3 |
| 10:18:21.334 | `rerun_script` (three component values) | |
| 10:18:21.335 | Beginning script thread | runner #4 |
| 10:18:21.349 | Memory cache HIT, unpickle polars frame | runner #4 reads the heatmap cache |
| ~10:18:21.36 | access violation | process gone; browser shows CONNECTING |

Crash #3 (default configuration, separate server on another port) has the identical stack after a click on a Scan Table row, with four script threads started between 10:26:44.046 and 10:26:45.074. Crash #1 (10:14:20, before faulthandler was enabled) ended the log inside the same `render_heatmap` cache path.

## What was tested

| Scenario (Playwright, Edge) | Clicks | Result |
|---|---|---|
| Random sidebar walk, no workflow running (5 s sidebar auto-rerun ticking) | 30 | no white screen, 0 page errors |
| Viewer <-> Workflow while the Run tab's 1 s `st.rerun()` loop runs a FLASHDeconv job | 24 | no white screen |
| User-guide flow (Configure widget, Start Workflow, Viewer / Layout / Download), fresh `__pycache__` | 70 | no white screen |
| Viewer with an experiment selected (grid visible), navigate / click rows | 2, 24 | server crash twice, identical stack (plus crash #1 earlier) |
| Same, with `--runner.fastReruns false` | 40 | no crash, but the Viewer rendered 1 of 6 cells and every later navigation click was ignored |
| Standalone: 4 threads unpickling polars frames + 2 streaming collects + 2 pyarrow readers, 30 s | ~120k ops | no crash; the failure needs Streamlit's stop-and-restart interleaving, not just threads |

The purely front-end white screen known in Streamlit >= 1.44 (a delta-path mismatch such as `Bad 'setIn' index` thrown inside a React `setState` updater with no root error boundary, which unmounts the whole app) did **not** occur in about 250 navigations. FLASHApp uses the upstream-implicated patterns (a `run_every` fragment on every page, `st.spinner` / `st.status` with `time.sleep` + `st.rerun`), so the risk exists, but it is not what the reproduction found.

## Secondary findings

- **Source watcher purges modules mid-run.** `src/`, `src/render`, `src/parse`, `src/common` and the `content` folders have no `__init__.py`, so Streamlit registers the directories as watched paths. Every `.pyc` written under them (first imports, the spawned workflow child) counts as a source change: 67 change signals in one seven-minute run. Each one shows the "File change - Rerun / Always rerun" prompt and deletes all `src.*` modules from `sys.modules` from the watchdog thread. One run hit `KeyError: 'src.render.compression'` inside the import machinery because of that race.
- **Table filter dialog leaks into the page.** `TabulatorTable.vue` teleports its filter dialog into the parent Streamlit document. Navigating with it open leaves a fixed, full-viewport, 50 % black backdrop (z-index 9999) behind; only a refresh removes it reliably.
- **Run tab reruns once per second while a job runs** (`time.sleep(1); st.rerun()`), on every tab of the Workflow page. Not a crash cause on its own, but it multiplies the runs that can overlap with user actions.
- **No safety net for front-end mismatches** since Streamlit 1.44: a delta-tree mismatch blanks the page instead of showing "Bad message format".

## Fixes and whether they work alone

The crash needs two things at once: overlapping script runs, and shared native state that is unsafe under overlap. Removing either one stops the crash.

- **Fix 1: stop the echo cascade in the Vue component** (drop `immediate: true`; do not echo Python's own state back; compare with the last payload sent). Alone it removes most overlapping runs during grid render, so crashes become rare but not zero: any rerun that lands mid-run (a row click, a sidebar click during a Viewer run, the 5-second auto-rerun) can still start a second thread. Crash #3 was triggered by a single real row click.
- **Fix 2a: stop pickling polars DataFrames through `st.cache_data`** (`update.py:34`; return pandas or Arrow, or cache bytes yourself). Alone it removes the exact crash site seen in all three stacks. Residual risk remains because concurrent threads still share the LazyFrames and pyarrow datasets in session state; the standalone test did not crash those, so this is likely sufficient for the observed crash but not provably for every path.
- **Fix 2b: per-session lock around `render_grid`.** Alone it removes the overlap itself for the Viewer, so it fixes the crash class rather than one crash site. Cost: a new run waits until the old one reaches its next Streamlit call, which can be a second or two inside `initialize_data`.
- **`runner.fastReruns=false`** also removes the overlap, but alone it stalls the Viewer (tested: one of six cells rendered, navigation ignored). Possibly viable together with Fix 1; untested.
- **Fix 3: `__init__.py` files, or `--server.fileWatcherType none` in production** (Windows `.bat`, Docker entrypoint, Kubernetes command). Independent of the crash; fixes the "File change" prompts and the mid-run module purge.
- **Fix 4: clean up the teleported filter dialog** (remove the backdrop on the iframe's `pagehide` / `beforeunload` and on RENDER events with the dialog closed, or render the dialog inside the iframe). Independent of the crash.
- **Fix 5: `PYTHONFAULTHANDLER=1` in production, pin polars** (currently unpinned as `polars>=1.0.0`). Diagnostics and hygiene only; consider reporting the `DataFrame.deserialize` access violation to polars with the stack above.

Recommended combination: Fix 1 plus Fix 2b is the smallest set that eliminates the crash class and keeps the Viewer responsive; add Fix 2a as defence in depth. Fixes 3, 4 and 5 address separate problems and can ship independently.

## Reproduction

Start the app with faulthandler so a native crash leaves a stack in the log:

```
.venv\Scripts\python.exe -X faulthandler -m streamlit run app.py local --server.headless true --logger.level debug
```

Open the FLASHDeconv Viewer, pick an experiment, then click a table row or a sidebar link while the grid re-renders. The crash appeared after 2 and 24 clicks in the two faulthandler runs. The Playwright harness used for the investigation (`repro.py`, scenario `viewer_grid`) decodes Streamlit's WebSocket protobuf frames and captures console and page errors; `polars_thread_test.py` is the standalone concurrency check.

## Measured impact of Fix 1 + Fix 2b (2026-09-10)

Setup: branch `fix/viewer-rerun-cascade`, same machine, server started with `-X faulthandler`, `--server.headless true`, `--logger.level debug`, default `runner.fastReruns`; Edge driven by Playwright (`measure.py`, fresh browser context per trial, experiment `example_fd_20260910-095658`, default 6-cell layout). "Settled" = no WebSocket frame in either direction for 1.5 s and no running indicator; "script runs" = `NewSession` messages received in that window; "rerun requests" = `rerun_script` messages the browser sent. Values are median (min-max) over 5 trials.

Fix 1 = Vue component only sends `setComponentValue` for state that did not just arrive from Python and differs from the last payload sent (patch: `vue_fix1.patch`, applied to `openms-streamlit-vue-component@FVdeploy`, bundle rebuilt into `js-component/dist/`). Fix 2b = per-session lock around `render_grid` (`src/render/render.py`).

### Initial load of the Viewer (T0 = clicking the experiment in the dropdown)

| Metric | Before | After (Fix 1 + 2b) |
|---|---|---|
| first grid cell visible (s) | 0.16 (0.14-0.17) | 0.15 (0.15-0.16) |
| last of 6 cells appears (s) | 0.37 (0-1.08) | 0.26 (0-0.73) |
| page settled (s) | 2.47 (2.24-2.92) | 2.06 (2.03-2.46) |
| script runs caused | 11 (8-12) | 5 (5-6) |
| rerun requests sent by browser | 20 (20-24) | 4 (4-4) |

### Navigation on the Viewer (grid loaded and settled)

| Metric | Before | After (Fix 1 + 2b) |
|---|---|---|
| Scan Table row click: settled (s) | 0.45 (0.45-0.53) | 0.27 (0.26-0.30) |
| row click: script runs | 3 (3-4) | 2 (2-2) |
| row click: rerun requests | 7 (7-7) | 1 (1-1) |
| row click updated the Mass Table | 5 / 5 trials | 5 / 5 trials |
| leave to Layout Manager: content visible (s) | 0.18 (0.16-0.18) | 0.18 (0.17-0.21) |
| leave: settled (s) | 0.16 (0.15-0.56) | 0.16 (0.14-1.60) |
| return to Viewer: all 6 cells back (s) | 0.19 (0.18-0.26) | 0.17 (0.17-0.23) |
| return to Viewer: settled (s) | 1.33 (1.30-1.37) | 0.19 (0.17-0.24) |
| return: script runs / rerun requests | 4 (4-5) / 7 (7-7) | 1 (1-1) / 1 (1-1) |

### Crash check (`repro.py viewer_grid`, 60 navigation clicks per run, faulthandler on)

| Run | Before | After (Fix 1 + 2b) |
|---|---|---|
| measurement session | server crashed in trial 3 (row click), same polars `deserialize` stack | 5 / 5 trials clean |
| seed 3 | crashed after 4 clicks | 60 clicks, no crash |
| seed 5 | crashed after 6 clicks (blank main area) | 60 clicks, no crash |
| seed 3 with `--runner.fastReruns false` | not run | 60 clicks, no crash, but the Viewer still stalls after returning to it (1 of 6 cells, later navigation ignored) |

No `fatal exception` appears in any AFTER server log.

### Interpretation

- The echo cascade was the main source of overlapping runs: rerun requests per Viewer load drop from 20-24 to 4, per row click from 7 to 1, per return to the Viewer from 7 to 1. The remaining 4 requests on load are the experiment selectbox plus the default selections the tables and heatmap make on first render (genuine state changes, sent once each).
- Responsiveness improves or stays equal on every axis; the biggest user-visible gains are the return to the Viewer (settled 1.33 s -> 0.19 s) and the row click (0.45 s -> 0.27 s). Initial load settles about 0.4 s earlier; the first cell appears at the same time.
- The per-session lock (Fix 2b) added no measurable latency in these runs because after Fix 1 the runs rarely overlap; it is the safety net for the interactions that still do (row clicks, sidebar clicks during a run).
- `runner.fastReruns=false` remains unusable even with both fixes: it stops the crash but the Viewer stalls; keep the default.
- Pre-existing browser-side messages unchanged by the fix: one 404 for a static resource on load, and a Tabulator `Scroll Error - No matching row found` page error when the selected row is not in the table (both present before and after).
