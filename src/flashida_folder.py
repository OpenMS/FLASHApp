"""The method folder of a FLASHIda acquisition.

Methods live as JSON files in a folder the user sets (one file per method, the
file name without .json is the method's name). A run of the table names its
method file (``method_file``) and gets the file as it is when FLASHIda starts
for it. The built-in presets are written into the folder once (an existing
file is never overwritten), and ``.favorites.json`` in the folder lists the
methods shown first in every method choice.

Streamlit-free, so the workflow process and the tests use it.
"""

import json
import re
import time
from pathlib import Path

from src import flashida_method as fm

FAVORITES_FILE = '.favorites.json'
_BAD = re.compile(r'[\\/:*?"<>|]')


def usable(folder) -> bool:
    """True when ``folder`` names an existing folder."""
    text = str(folder or '').strip()
    return bool(text) and Path(text).is_dir()


def stem(name: str) -> str:
    """The file name (without .json) a method of this name is saved as."""
    name = (name or '').strip()
    if name.lower().endswith('.json'):
        name = name[:-5]
    return _BAD.sub('_', name).strip(' .')


def name_problem(name: str) -> str:
    """Why a method name cannot be a file name ('' when it can)."""
    name = (name or '').strip()
    if not name:
        return 'Enter a name for the method.'
    if _BAD.search(name):
        return 'A method name cannot hold \\ / : * ? " < > |.'
    if name.startswith('.'):
        return 'A method name cannot start with a dot.'
    return ''


def path(folder, method: str) -> Path:
    return Path(str(folder).strip(), f'{method}.json')


def list_methods(folder) -> list:
    """Names of the methods in the folder, sorted as Windows sorts (ignoring case)."""
    if not usable(folder):
        return []
    try:
        files = [p for p in Path(str(folder).strip()).iterdir()
                 if p.suffix.lower() == '.json' and p.is_file() and not p.name.startswith('.')]
    except OSError:
        return []
    return sorted((p.stem for p in files), key=str.casefold)


def exists(folder, method) -> bool:
    return bool(method) and usable(folder) and path(folder, method).is_file()


def read(folder, method: str) -> dict:
    """The method in ``<folder>/<method>.json``. Raises ValueError with a
    message for the user when it cannot be read."""
    p = path(folder, method)
    try:
        data = json.loads(p.read_text(encoding='utf-8-sig'))
    except OSError as e:
        raise ValueError(f"The method file '{p.name}' cannot be read: {e.strerror or e}") from e
    except ValueError as e:
        raise ValueError(f"The method file '{p.name}' is not valid JSON: {e}") from e
    if not isinstance(data, dict):
        raise ValueError(f"The method file '{p.name}' does not hold a method ({{...}}).")
    return data


def write(folder, method: str, data: dict, overwrite: bool = False) -> Path:
    """Saves a method as ``<folder>/<method>.json``. Raises FileExistsError
    when the file is there and ``overwrite`` is off."""
    p = path(folder, method)
    if p.exists() and not overwrite:
        raise FileExistsError(p.name)
    tmp = p.with_name(p.name + '.tmp')
    tmp.write_text(json.dumps(data, indent=2), encoding='utf-8')
    tmp.replace(p)
    return p


def seed_presets(folder, presets: dict) -> list:
    """Writes each preset whose file is not in the folder yet; the names written."""
    if not usable(folder):
        return []
    written = []
    have = {m.casefold() for m in list_methods(folder)}
    for name, data in presets.items():
        s = stem(name)
        if s and s.casefold() not in have:
            write(folder, s, data)
            have.add(s.casefold())
            written.append(s)
    return written


# --- favorites --------------------------------------------------------------------

def favorites(folder) -> list:
    if not usable(folder):
        return []
    try:
        data = json.loads(Path(str(folder).strip(), FAVORITES_FILE).read_text(encoding='utf-8'))
    except (OSError, ValueError):
        return []
    return [str(m) for m in data if isinstance(m, str)] if isinstance(data, list) else []


def toggle_favorite(folder, method: str) -> bool:
    """Adds or removes a favorite; True when it is a favorite now."""
    favs = favorites(folder)
    on = method not in favs
    favs = favs + [method] if on else [m for m in favs if m != method]
    Path(str(folder).strip(), FAVORITES_FILE).write_text(json.dumps(favs, indent=2), encoding='utf-8')
    return on


def ordered(methods: list, favs: list) -> list:
    """Favorites first (in the order of ``methods``), then the others."""
    fav = set(favs)
    return [m for m in methods if m in fav] + [m for m in methods if m not in fav]


# --- runs -------------------------------------------------------------------------

def run_method(run: dict, folder) -> dict | None:
    """The method a run gets: its method file, or (a table saved before the
    method folder) the method kept in the run. None when there is none."""
    name = run.get('method_file')
    if name:
        try:
            return read(folder, name) if usable(folder) else None
        except ValueError:
            return None
    return run.get('method') if isinstance(run.get('method'), dict) else None


def adopt_embedded(runs: list, folder) -> tuple:
    """(runs, written): runs that still keep their method in the table get a
    method file instead. A method equal to a file in the folder uses that
    file; any other is saved under its preset's name (or the run's name), with
    a number added when that name is taken."""
    if not usable(folder):
        return runs, []
    contents = {}
    for m in list_methods(folder):
        try:
            contents[m] = read(folder, m)
        except ValueError:
            pass
    written, out = [], []
    for run in runs:
        method = run.get('method')
        if run.get('method_file') or not isinstance(method, dict):
            out.append(run)
            continue
        match = next((m for m, data in contents.items() if data == method), None)
        if match is None:
            base = stem(run.get('preset') or run.get('name') or 'method') or 'method'
            match, n = base, 2
            while match.casefold() in {m.casefold() for m in contents}:
                match, n = f'{base} ({n})', n + 1
            write(folder, match, method)
            contents[match] = method
            written.append(match)
        out.append({'name': run.get('name', ''), 'method_file': match})
    return out, written


def method_problems(runs: list, folder) -> list:
    """[(index, text)] of runs whose method cannot be used."""
    out = []
    for i, run in enumerate(runs):
        name = run.get('method_file')
        if not name:
            if not isinstance(run.get('method'), dict):
                out.append((i, 'Pick a method.'))
        elif not usable(folder):
            out.append((i, 'Set the method folder.'))
        elif not exists(folder, name):
            out.append((i, f"The method file '{name}.json' is not in the method folder."))
        else:
            try:
                read(folder, name)
            except ValueError as e:
                out.append((i, str(e)))
    return out


# --- moving and duplicating ---------------------------------------------------------

_TRAILING = re.compile(r'^(.*?)(\d+)$')


def next_name(name: str, taken) -> str:
    """A name for a copy of a run: its trailing number counted up (R1 → R2,
    rep01 → rep02) until it is free, else ``_2``, ``_3``, … added."""
    used = {t.casefold() for t in taken}
    m = _TRAILING.match(name or '')
    if m:
        head, digits = m.group(1), m.group(2)
        n = int(digits)
        while True:
            n += 1
            cand = f'{head}{n:0{len(digits)}d}'
            if cand.casefold() not in used:
                return cand
    n = 2
    while f'{name}_{n}'.casefold() in used:
        n += 1
    return f'{name}_{n}'


def duplicate(runs: list, index: int) -> list:
    """The runs with a copy of run ``index`` right after it, under the next free name."""
    out = [dict(r) for r in runs]
    copy = dict(out[index])
    copy['name'] = next_name(copy.get('name', ''), [r.get('name', '') for r in out])
    out.insert(index + 1, copy)
    return out


def move(runs: list, index: int, step: int) -> list:
    """The runs with run ``index`` moved ``step`` places (-1 up, +1 down)."""
    out = list(runs)
    j = index + step
    if 0 <= index < len(out) and 0 <= j < len(out):
        out[index], out[j] = out[j], out[index]
    return out


# --- duration ---------------------------------------------------------------------

def duration(method: dict | None):
    """A method's ``global.duration`` in minutes, or None."""
    value = fm.get(method or {}, 'global.duration')
    return float(value) if isinstance(value, (int, float)) and not isinstance(value, bool) else None


def total_minutes(runs: list, folder) -> tuple:
    """(minutes, unknown): the summed duration of the runs' methods and how
    many runs have none."""
    total, unknown = 0.0, 0
    for run in runs:
        d = duration(run_method(run, folder))
        if d is None:
            unknown += 1
        else:
            total += d
    return total, unknown


def format_minutes(minutes: float) -> str:
    m = int(round(minutes))
    return f'{m // 60} h {m % 60:02d} min' if m >= 60 else f'{m} min'


def format_end(minutes: float, now: float | None = None) -> str:
    """When something that takes ``minutes`` from ``now`` ends, as 'Wed 07:30'."""
    end = (time.time() if now is None else now) + minutes * 60
    return time.strftime('%a %H:%M', time.localtime(end))
