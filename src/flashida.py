"""FLASHIda: intelligent data acquisition for top-down proteomics.

FLASHIda runs on the instrument PC next to a Thermo Orbitrap. The acquisition
page holds a table of runs, each a raw file name with its method
(src/flashida_method.py). While it runs, the newest Xcalibur sequence (.sld) in
the raw file folder shows which runs are scheduled, and FLASHIda is started
with a run's method when its raw file appears (src/Workflow.py:IdaWorkflow). The
FLASHIda Viewer follows the logs it writes (src/parse/ida_log.py). Both are
shown only in the local app, and only once developer mode is unlocked
(``available()``).

Everything but ``available()`` is Streamlit-free, so the workflow process can
use it.
"""

import json
import re
import time
from pathlib import Path


def raw_names(folder) -> set:
    """Names (without .raw) of the raw files in a folder."""
    try:
        return {p.stem for p in Path(folder).iterdir()
                if p.suffix.lower() == '.raw' and p.is_file()}
    except OSError:
        return set()


def newest_sequence(folder):
    """The newest Xcalibur sequence (.sld) in a folder, or None."""
    try:
        slds = [p for p in Path(folder).iterdir() if p.suffix.lower() == '.sld' and p.is_file()]
    except OSError:
        return None
    return max(slds, key=lambda p: p.stat().st_mtime, default=None)


# Every row of an Xcalibur sequence keeps its raw file name in a JSON blob
# (UTF-16 in the .sld): "OriginalRawFileNameWithoutExtension":"<name>","DisplayRowNumber":<n>
_SLD_ROW = re.compile(
    r'"OriginalRawFileNameWithoutExtension":"((?:[^"\\]|\\.)*)","DisplayRowNumber":(\d+)')


def sequence_rows(path) -> list:
    """[{name, row}] of the raw files an Xcalibur sequence (.sld) writes, in order."""
    return sequence_rows_from_bytes(Path(path).read_bytes())


def sequence_rows_from_bytes(data: bytes) -> list:
    """``sequence_rows`` of an .sld already read, e.g. one the user uploaded."""
    text = data.decode('utf-16-le', errors='ignore')
    rows = []
    for m in _SLD_ROW.finditer(text):
        try:
            name = json.loads('"' + m.group(1) + '"')
        except ValueError:
            name = m.group(1)
        rows.append({'name': name, 'row': int(m.group(2))})
    return rows


def scheduled(rows: list, names) -> dict:
    """{run name: sequence row} of the runs the sequence lists (names compared
    as Windows does, ignoring case)."""
    by_name = {r['name'].casefold(): r['row'] for r in rows}
    return {n: by_name[n.casefold()] for n in names if n.casefold() in by_name}


# What the running acquisition reports to the page.
STATUS_FILE = 'ida-status.json'


def read_status(workflow_dir) -> dict:
    try:
        return json.loads(Path(workflow_dir, STATUS_FILE).read_text(encoding='utf-8'))
    except (OSError, ValueError):
        return {}


def write_status(workflow_dir, status: dict) -> None:
    path = Path(workflow_dir, STATUS_FILE)
    status = dict(status, updated=time.time())
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(status, indent=2), encoding='utf-8')
    tmp.replace(path)


def available() -> bool:
    """FLASHIda pages are shown in the local app once developer mode is on."""
    import streamlit as st
    from src.common.developer import is_developer

    settings = st.session_state.get('settings') or {}
    return not settings.get('online_deployment', True) and is_developer()


# The log folders of the methods FLASHIda was started with, newest last, in
# the workflow directory. The FLASHIda Viewer offers them.
LOG_FOLDERS = 'log-folders.txt'


def method_log_dir(method_path, working_dir) -> Path:
    """Folder FLASHIda writes the run folders of a method to: its
    runtime.log_dir, relative to FLASHIda's working directory (the default)."""
    try:
        method = json.loads(Path(method_path).read_text(encoding='utf-8-sig'))
        log_dir = str((method.get('runtime') or {}).get('log_dir') or '')
    except (OSError, ValueError, AttributeError):
        log_dir = ''
    path = Path(log_dir) if log_dir else Path('.')
    return (path if path.is_absolute() else Path(working_dir, path)).resolve()


def remember_log_folder(workflow_dir, folder) -> None:
    file = Path(workflow_dir, LOG_FOLDERS)
    try:
        folders = [f for f in file.read_text(encoding='utf-8').splitlines() if f.strip()]
    except OSError:
        folders = []
    folders = [f for f in folders if f != str(folder)] + [str(folder)]
    file.write_text('\n'.join(folders[-20:]) + '\n', encoding='utf-8')


def log_folders(workflow_dir) -> list:
    """Log folders FLASHIda was started with, newest first."""
    try:
        lines = Path(workflow_dir, LOG_FOLDERS).read_text(encoding='utf-8').splitlines()
    except OSError:
        return []
    return [line for line in reversed(lines) if line.strip()]
