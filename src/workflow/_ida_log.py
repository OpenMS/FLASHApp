"""
Pure helpers for the optional IDA-log -> mzML passthrough to FLASHDeconv.

Kept streamlit- and pyopenms-free so the log discovery and filename auto-match
logic can be unit-tested without booting a Streamlit runtime (mirrors the
src/workflow/_log_status.py pattern).
"""

from pathlib import Path
from os.path import basename, splitext, exists

# Upload widget key + "no log" sentinel, shared by the UI and execution.
IDA_LOG_KEY = "ida-log-files"
IDA_NONE = "(none)"


def available_ida_logs(workflow_dir) -> list:
    """Return all uploaded IDA ``.log`` file paths for a workflow directory.

    Covers both copy-mode uploads (files placed directly in
    ``input-files/ida-log-files/``) and local mode (existing paths listed in
    ``external_files.txt``).
    """
    log_dir = Path(workflow_dir, "input-files", IDA_LOG_KEY)
    if not log_dir.exists():
        return []
    logs = [
        str(f) for f in log_dir.iterdir()
        if f.name.endswith(".log") and f.name != "external_files.txt"
    ]
    external_files = Path(log_dir, "external_files.txt")
    if external_files.exists():
        logs += [
            line.strip() for line in external_files.read_text().splitlines()
            if line.strip().endswith(".log") and exists(line.strip())
        ]
    return logs


def auto_match_log(mzml_name: str, logs: list) -> str:
    """Return the log whose file-name stem matches ``mzml_name``, else ``IDA_NONE``.

    Matching is by file-name stem (extension stripped), so ``sample1.mzML`` maps
    to ``.../sample1.log`` regardless of directory.
    """
    target = splitext(basename(mzml_name))[0]
    for log in logs:
        if splitext(basename(log))[0] == target:
            return log
    return IDA_NONE


def selected_mzml_files(session_state, params, param_prefix) -> list:
    """Return the *live* mzML selection, preferring the Streamlit widget value.

    ``params`` is a snapshot loaded once at workflow construction and is **not**
    refreshed on Streamlit fragment reruns, so it lags the user's current
    multiselect choice (e.g. right after selecting files and toggling a reactive
    checkbox). The multiselect widget always mirrors its current value into
    ``session_state[f"{param_prefix}mzML-files"]``, so prefer that and fall back
    to the persisted ``params`` value only when session_state has nothing.
    """
    return session_state.get(f"{param_prefix}mzML-files") or params.get("mzML-files", [])
