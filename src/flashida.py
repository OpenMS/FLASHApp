"""FLASHIda: intelligent data acquisition for top-down proteomics.

FLASHIda runs on the instrument PC next to a Thermo Orbitrap. The app starts it
for every raw file the acquisition writes (src/Workflow.py:IdaWorkflow), and the
FLASHIda Viewer follows the logs it writes (src/parse/ida_log.py). Both are
shown only in the local app, and only once developer mode is unlocked
(``available()``).

Everything but ``available()`` is Streamlit-free, so the workflow process can
use it.
"""

import re
from pathlib import Path

# Raw files the acquisition writes: <anything>FLASHIda_<method>[_<anything>].raw.
# FLASHIda is started with <method>.json from the method folder.
_RAW = re.compile(r'.*FLASHIda_([^.]+)\.raw$', re.IGNORECASE)


def parse_raw_name(file_name: str):
    """The method named in a FLASHIda raw file name, or None."""
    match = _RAW.match(file_name)
    if match is None:
        return None
    return match.group(1).partition('_')[0] or None


def find_raws(folder) -> dict:
    """{raw file name: method} of the FLASHIda raw files in a folder."""
    raws = {}
    for path in Path(folder).iterdir():
        if path.is_file():
            parsed = parse_raw_name(path.name)
            if parsed is not None:
                raws[path.name] = parsed
    return raws


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
    import json

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
