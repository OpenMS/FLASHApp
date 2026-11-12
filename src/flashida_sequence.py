"""Sequences of FLASHIda runs: export and import, names from an Xcalibur
sequence, and names made from a pattern.

A run is ``{name, method_file}``: its raw file name and the name of its method
in the method folder (src/flashida_folder.py). Tables saved before the method
folder kept the method in the run (``{name, preset, method}``); those are read
too. Streamlit-free, so it is tested on its own.
"""

import csv
import datetime
import io
import json
import re

from src import flashida
from src import flashida_folder as ff

FORMAT = 'flashida-runs'
VERSION = 2

# Built-in preset that the pattern dialog proposes for wash runs.
WASH_PRESET = 'Wash'


# --- export and import ------------------------------------------------------------

def export_runs(runs: list, folder=None) -> str:
    """The runs table as a JSON sequence file. Each run names its method file
    and carries the method too, so the file can be imported on a PC whose
    method folder does not have it."""
    out = []
    for r in runs:
        entry = {'name': r.get('name', '')}
        if r.get('method_file'):
            entry['method_file'] = r['method_file']
        method = ff.run_method(r, folder)
        if method is not None:
            entry['method'] = method
        out.append(entry)
    return json.dumps({'format': FORMAT, 'version': VERSION, 'runs': out}, indent=2)


def export_csv(runs: list) -> str:
    """The run names (and their methods) as CSV, in table order."""
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator='\r\n')
    w.writerow(['Row', 'Raw file name', 'Method'])
    for i, r in enumerate(runs, 1):
        w.writerow([i, r.get('name', ''), r.get('method_file') or r.get('preset') or ''])
    return buf.getvalue()


# Xcalibur's sequence as CSV (Sequence Setup: File > Export / Import): a first
# line "Bracket Type=4," and then a header row and one row per sample.
XCALIBUR_COLUMNS = ['Sample Type', 'File Name', 'Sample ID', 'Path', 'Instrument Method', 'Process Method',
                    'Calibration File', 'Position', 'Inj Vol', 'Level', 'Sample Wt', 'Sample Vol', 'ISTD Amt',
                    'Dil Factor', 'L1 Study', 'L2 Client', 'L3 Laboratory', 'L4 Company', 'L5 Phone',
                    'Comment', 'Sample Name']


def export_xcalibur_csv(runs: list, path: str = '', instrument_method: str = '') -> str:
    """The runs as an Xcalibur sequence CSV: one Unknown sample per run, its
    raw file name as File Name and Sample Name, writing into ``path`` with the
    instrument method ``instrument_method`` (both may be left for Xcalibur)."""
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator='\r\n')
    w.writerow(['Bracket Type=4', ''])
    w.writerow(XCALIBUR_COLUMNS)
    for r in runs:
        row = dict.fromkeys(XCALIBUR_COLUMNS, '')
        row.update({'Sample Type': 'Unknown', 'File Name': r.get('name', ''), 'Path': path,
                    'Instrument Method': instrument_method, 'Sample Wt': '0', 'Sample Vol': '0',
                    'ISTD Amt': '0', 'Dil Factor': '1', 'Sample Name': r.get('name', '')})
        w.writerow([row[c] for c in XCALIBUR_COLUMNS])
    return buf.getvalue()


def csv_runs(data) -> tuple:
    """(runs, has methods) of a CSV sequence: an Xcalibur sequence CSV (its
    File Name column, no methods) or the table Export sequence writes (Raw file
    name and Method). Raises ValueError for any other CSV."""
    if isinstance(data, bytes):
        try:
            data = data.decode('utf-8-sig')
        except UnicodeDecodeError:
            data = data.decode('cp1252')
    rows = list(csv.reader(io.StringIO(data)))
    for i, row in enumerate(rows):
        cells = [c.strip() for c in row]
        for column, method in (('File Name', None), ('Raw file name', 'Method')):
            if column in cells:
                at = cells.index(column)
                m_at = cells.index(method) if method in cells else None
                runs = [{'name': r[at].strip(),
                         'method_file': (r[m_at].strip() or None) if m_at is not None and len(r) > m_at else None}
                        for r in rows[i + 1:] if len(r) > at and r[at].strip()]
                return runs, m_at is not None
    raise ValueError('This CSV is neither an Xcalibur sequence (no "File Name" column) nor a table '
                     'exported here (no "Raw file name" column).')


def parse_runs(text) -> list:
    """The runs of a JSON sequence file (or a bare list of runs), as
    ``{name, method_file?, method?, preset?}``.

    Raises ValueError with a message for the user when the file is not one."""
    if isinstance(text, bytes):
        text = text.decode('utf-8-sig')
    try:
        data = json.loads(text)
    except ValueError as e:
        raise ValueError(f'Not valid JSON: {e}') from e
    if isinstance(data, dict):
        if data.get('format') not in (None, FORMAT):
            raise ValueError(f"This is not a FLASHIda sequence file (format {data.get('format')!r}).")
        data = data.get('runs')
    if not isinstance(data, list):
        raise ValueError('A sequence file holds a list of runs.')
    runs = []
    for i, r in enumerate(data, 1):
        if not isinstance(r, dict):
            raise ValueError(f'Run {i} is not a run.')
        method_file = r.get('method_file')
        method = r.get('method')
        if not (isinstance(method_file, str) and method_file.strip()) and not isinstance(method, dict):
            raise ValueError(f'Run {i} has no method.')
        run = {'name': str(r.get('name') or '').strip()}
        if isinstance(method_file, str) and method_file.strip():
            run['method_file'] = ff.stem(method_file)
        if isinstance(method, dict):
            run['method'] = method
        if isinstance(r.get('preset'), str):
            run['preset'] = r['preset']
        runs.append(run)
    return runs


def into_folder(runs: list, folder) -> tuple:
    """(runs, written, missing): imported runs as ``{name, method_file}``.
    A method file the folder lacks is written from the method the file
    carries; a run whose method is neither there nor carried is listed in
    ``missing`` (it keeps its method file name, which the table marks)."""
    written, missing, out, legacy = [], [], [], []
    for r in runs:
        name = r.get('method_file')
        if name:
            if not ff.exists(folder, name):
                if isinstance(r.get('method'), dict) and ff.usable(folder):
                    ff.write(folder, name, r['method'])
                    written.append(name)
                else:
                    missing.append(name)
            out.append({'name': r['name'], 'method_file': name})
        else:
            legacy.append(len(out))
            out.append(r)
    if legacy:
        adopted, more = ff.adopt_embedded([out[i] for i in legacy], folder)
        for i, run in zip(legacy, adopted):
            out[i] = run
        written += more
    return out, written, sorted(set(missing), key=str.casefold)


def sld_names(data: bytes) -> list:
    """The raw file names an Xcalibur sequence (.sld, as bytes) writes, in order."""
    return [r['name'] for r in flashida.sequence_rows_from_bytes(data)]


def runs_from_names(names, method_file: str | None) -> list:
    """One run per name, all with the same method file."""
    return [{'name': n, 'method_file': method_file} for n in names]


def combine(existing: list, new: list, mode: str) -> tuple:
    """(runs, duplicates): ``new`` appended to or replacing ``existing``
    (``mode`` 'append' or 'replace'). ``duplicates`` are the names that occur
    more than once afterwards (ignoring case, as Windows does); they are kept,
    and the table marks them."""
    if mode not in ('append', 'replace'):
        raise ValueError(f'Unknown mode {mode!r}')
    runs = (list(existing) if mode == 'append' else []) + list(new)
    return runs, duplicate_names(runs)


def duplicate_names(runs: list) -> list:
    """Names (as first written) that more than one run has, in table order."""
    seen, dup = {}, []
    for r in runs:
        name = (r.get('name') or '').strip()
        if not name:
            continue
        key = name.casefold()
        if key in seen and seen[key] not in dup:
            dup.append(seen[key])
        seen.setdefault(key, name)
    return dup


# --- names from a pattern ---------------------------------------------------------

REPLICATE_FORMATS = {'R{n}': 'R1', 'R{nn}': 'R01', 'rep{n}': 'rep1', 'rep{nn}': 'rep01'}
ORDERS = {'sample': 'All replicates of a sample, then the next',
          'replicate': 'Alternate: replicate 1 of every sample, then 2, …'}
DEFAULT_WASH_NAME = '{prefix}_{date}_wash_{n}'
SEPARATOR = '_'


def today() -> str:
    return datetime.date.today().strftime('%Y%m%d')


def replicate_label(n: int, fmt: str) -> str:
    """``fmt`` with {nn} as two digits and {n} as the plain number."""
    return fmt.replace('{nn}', f'{n:02d}').replace('{n}', str(n))


def join(parts) -> str:
    """The non-empty parts joined by '_' (an empty placeholder leaves no '__')."""
    return SEPARATOR.join(str(p).strip() for p in parts if str(p or '').strip())


def parse_samples(text: str, comment: str = '', replicates: int = 1) -> list:
    """[{sample, comment, replicates}] from one sample per line; a line may
    set its own values: ``Ecoli; comment=top3; replicates=3``.

    Raises ValueError naming the line when a value cannot be read."""
    samples = []
    for no, line in enumerate((text or '').splitlines(), 1):
        parts = [p.strip() for p in line.split(';')]
        if not parts[0]:
            if any(parts):
                raise ValueError(f'Line {no}: the sample name comes first.')
            continue
        s = {'sample': parts[0], 'comment': comment, 'replicates': replicates}
        for part in parts[1:]:
            if not part:
                continue
            key, eq, value = (x.strip() for x in part.partition('='))
            key = key.lower()
            if not eq or key not in ('comment', 'replicates'):
                raise ValueError(f'Line {no}: write comment=… or replicates=… after a ";".')
            if key == 'comment':
                s['comment'] = value
            else:
                try:
                    s['replicates'] = int(value)
                except ValueError:
                    raise ValueError(f'Line {no}: replicates must be a whole number.') from None
                if s['replicates'] < 1:
                    raise ValueError(f'Line {no}: at least 1 replicate.')
        samples.append(s)
    return samples


def pattern_names(samples: list, prefix: str = '', date: str = '', suffix: str = '',
                  replicate_format: str = 'R{n}', order: str = 'sample',
                  wash_every: int = 0, wash_name: str = DEFAULT_WASH_NAME) -> list:
    """[{name, wash}] of the runs a pattern makes:
    ``prefix_date_sample_comment_replicate_suffix``, empty parts left out.

    With ``wash_every`` > 0 a wash run follows every that many sample runs
    (none after the last); ``wash_name`` takes {prefix}, {date}, {suffix} and
    {n}, the wash number."""
    if order not in ORDERS:
        raise ValueError(f'Unknown order {order!r}')
    if order == 'sample':
        runs = [(s, k) for s in samples for k in range(1, s['replicates'] + 1)]
    else:
        most = max((s['replicates'] for s in samples), default=0)
        runs = [(s, k) for k in range(1, most + 1) for s in samples if k <= s['replicates']]

    out, washes = [], 0
    for i, (s, k) in enumerate(runs, 1):
        out.append({'name': join([prefix, date, s['sample'], s['comment'],
                                  replicate_label(k, replicate_format), suffix]),
                    'wash': False})
        if wash_every and wash_every > 0 and i % wash_every == 0 and i < len(runs):
            washes += 1
            out.append({'name': wash_run_name(wash_name, prefix, date, suffix, washes), 'wash': True})
    return out


def wash_run_name(pattern: str, prefix: str, date: str, suffix: str, n: int) -> str:
    """A wash run's name: the placeholders filled, the parts an empty one
    leaves behind joined again without empty pieces."""
    values = {'prefix': prefix, 'date': date, 'suffix': suffix, 'n': str(n)}
    filled = _PLACEHOLDER.sub(lambda m: values.get(m.group(1), m.group(0)), pattern or DEFAULT_WASH_NAME)
    return join(filled.split(SEPARATOR))


_PLACEHOLDER = re.compile(r'\{(\w+)\}')


def unknown_placeholders(pattern: str) -> list:
    """Placeholders of a wash name that are not filled in."""
    return [p for p in _PLACEHOLDER.findall(pattern or '') if p not in ('prefix', 'date', 'suffix', 'n')]
