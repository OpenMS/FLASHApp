"""Run the example analyses with the TOPP tools on PATH, or compare two such runs.

    python .github/scripts/compare_openms_builds.py run OUT
        Runs every example of example-data/build_examples.py exactly as the Run
        page does and keeps the stored results (cache/files/<example id>/) in
        OUT. The example inputs are used as they are in the repository.

    python .github/scripts/compare_openms_builds.py compare OLD NEW [--md FILE]
        Compares the stored results of two runs file by file: tables column by
        column (numbers as numbers), mzML spectra peak by peak, everything else
        line by line, and writes a Markdown report.
"""

import argparse
import gzip
import importlib.util
import json
import math
import pickle
import shutil
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

REL_TOL = 1e-6   # differences below this (relative) count as printing digits
SAMPLE = 3       # example differences shown per column


def load_examples():
    spec = importlib.util.spec_from_file_location('build_examples', REPO / 'example-data' / 'build_examples.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run(out):
    ex = load_examples()
    from src.workflow.FileManager import FileManager
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        workspace = Path(tmp, 'workspace')
        cache = workspace / 'cache'
        kept = []
        for example in ex.EXAMPLES:
            print(f"== {example['files']}", flush=True)
            run_ids, comparison = ex.run(example, workspace)
            shutil.copytree(workspace / 'flashtnt' / 'logs', out / 'logs' / Path(example['files'][0]).stem,
                            dirs_exist_ok=True)
            for old_id, entry in zip(run_ids, example['runs']):
                if entry is None:
                    FileManager(workspace / 'flashtnt', cache).remove_results(old_id)
                    continue
                ex.rename(cache, old_id, *entry)
                kept.append(entry[0])
            if comparison is not None:
                ex.rename(cache, comparison, *example['comparison'])
                if 'comparison_settings' in example:
                    from src.quant.session import QuantSession
                    session = QuantSession(FileManager(workspace / 'flashtnt', cache), example['comparison'][0])
                    session.save_settings({**session.settings, **example['comparison_settings']})
                    # the comparison chain runs again with the settings
                    session = QuantSession(FileManager(workspace / 'flashtnt', cache), example['comparison'][0])
                kept.append(example['comparison'][0])
        for dataset_id in kept:
            shutil.copytree(cache / 'files' / dataset_id, out / 'files' / dataset_id, dirs_exist_ok=True)
    print(f'kept {kept}')


# ---------------------------------------------------------------- comparison

class Diff:
    """Differences of one column (or of a file's values)."""

    def __init__(self):
        self.values = 0
        self.differ = 0        # beyond REL_TOL
        self.digits = 0        # numbers that differ below REL_TOL
        self.max_rel = 0.0
        self.max_abs = 0.0
        self.samples = []

    def number(self, a, b):
        self.values += 1
        if a == b or (isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b)):
            return
        if isinstance(a, float) and math.isnan(a) or isinstance(b, float) and math.isnan(b):
            self.other(a, b)
            return
        d = abs(a - b)
        rel = d / max(abs(a), abs(b), 1e-300)
        self.max_abs = max(self.max_abs, d)
        self.max_rel = max(self.max_rel, rel)
        if rel <= REL_TOL:
            self.digits += 1
        else:
            self.differ += 1
            if len(self.samples) < SAMPLE:
                self.samples.append((a, b))

    def other(self, a, b):
        self.values += 1
        if a != b:
            self.differ += 1
            if len(self.samples) < SAMPLE:
                self.samples.append((a, b))

    def value(self, a, b):
        """Compare two python values recursively."""
        if isinstance(a, bool) or isinstance(b, bool):
            self.other(a, b)
        elif isinstance(a, (int, float)) and isinstance(b, (int, float)):
            self.number(float(a), float(b))
        elif isinstance(a, str) and isinstance(b, str):
            fa, fb = as_float(a), as_float(b)
            if fa is not None and fb is not None:
                self.number(fa, fb)
            else:
                self.other(a, b)
        elif isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
            if len(a) != len(b):
                self.other(f'<{len(a)} items>', f'<{len(b)} items>')
            for x, y in zip(a, b):
                self.value(x, y)
        elif isinstance(a, dict) and isinstance(b, dict):
            if set(a) != set(b):
                self.other(sorted(map(str, a)), sorted(map(str, b)))
            for k in set(a) & set(b):
                self.value(a[k], b[k])
        else:
            try:
                import numpy as np
                if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
                    self.value(np.asarray(a).tolist(), np.asarray(b).tolist())
                    return
            except ImportError:
                pass
            self.other(a, b)

    def text(self):
        if not self.differ and not self.digits:
            return ''
        out = []
        if self.digits:
            out.append(f'{self.digits} last-digit')
        if self.differ:
            s = '; '.join(f'`{short(a)}` → `{short(b)}`' for a, b in self.samples)
            out.append(f'**{self.differ} of {self.values} differ** (max rel {self.max_rel:.2g}) e.g. {s}')
        return ', '.join(out)


def short(v, n=60):
    s = repr(v) if not isinstance(v, str) else v
    return s if len(s) <= n else s[:n] + '…'


def as_float(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def compare_frames(a, b):
    """Column-by-column differences of two polars frames, aligned by row."""
    import polars as pl
    notes = []
    if a.height != b.height:
        notes.append(f'**rows {a.height} → {b.height}**')
    ca, cb = list(a.columns), list(b.columns)
    if ca != cb:
        gone, new = [c for c in ca if c not in cb], [c for c in cb if c not in ca]
        if gone or new:
            notes.append(f'**columns removed {gone}, added {new}**')
        else:
            notes.append('columns reordered')
    n = min(a.height, b.height)
    cols = {}
    for c in [c for c in ca if c in cb]:
        sa, sb = a[c].head(n), b[c].head(n)
        if sa.dtype != sb.dtype:
            notes.append(f'`{c}` type {sa.dtype} → {sb.dtype}')
        d = Diff()
        if sa.dtype.is_numeric() and sb.dtype.is_numeric():
            import numpy as np
            x = sa.cast(pl.Float64).to_numpy()
            y = sb.cast(pl.Float64).to_numpy()
            both_nan = np.isnan(x) & np.isnan(y)
            eq = (x == y) | both_nan
            d.values = n
            idx = np.where(~eq)[0]
            for i in idx:
                d.values -= 1
                d.number(float(x[i]), float(y[i]))
        else:
            if sa.dtype == sb.dtype and sa.dtype in (pl.Utf8, pl.Boolean, pl.Categorical) or sa.dtype.is_temporal():
                eq = (sa == sb).fill_null(False) | (sa.is_null() & sb.is_null())
                d.values = n - int((~eq).sum())
                idx = (~eq).arg_true().to_list()
            else:
                idx = range(min(n, 20000))
            la, lb = sa.to_list(), sb.to_list()
            for i in idx:
                d.value(la[i], lb[i])
        t = d.text()
        if t:
            cols[c] = t
    return notes, cols


def read_table(path):
    import polars as pl
    if path.suffix == '.pq':
        return pl.read_parquet(path)
    # FLASH* tsv: header in the first line; keep everything as text first
    return pl.read_csv(path, separator='\t', infer_schema_length=0, quote_char=None,
                       truncate_ragged_lines=True, comment_prefix=None)


def compare_mzml(pa, pb):
    import numpy as np
    import pyopenms as oms
    ea, eb = oms.MSExperiment(), oms.MSExperiment()
    oms.MzMLFile().load(str(pa), ea)
    oms.MzMLFile().load(str(pb), eb)
    notes = []
    if ea.size() != eb.size():
        notes.append(f'**spectra {ea.size()} → {eb.size()}**')
    peaks, mz, inten, meta = Diff(), Diff(), Diff(), Diff()
    for sa, sb in zip(ea, eb):
        ma, ia = sa.get_peaks()
        mb, ib = sb.get_peaks()
        if len(ma) != len(mb):
            peaks.other(len(ma), len(mb))
            continue
        peaks.values += 1
        for x, y in ((ma, mb),):
            for i in np.where(x != y)[0]:
                mz.number(float(x[i]), float(y[i]))
        for i in np.where(ia != ib)[0]:
            inten.number(float(ia[i]), float(ib[i]))
        meta.number(sa.getRT(), sb.getRT())
        ka, kb = [], []
        sa.getKeys(ka)
        sb.getKeys(kb)
        for k in set(ka) | set(kb):
            meta.value(str(sa.getMetaValue(k)) if sa.metaValueExists(k) else None,
                       str(sb.getMetaValue(k)) if sb.metaValueExists(k) else None)
    cols = {}
    for name, d in (('peaks per spectrum', peaks), ('m/z', mz), ('intensity', inten), ('RT, meta values', meta)):
        if d.text():
            cols[name] = d.text()
    return notes, cols


def compare_lines(pa, pb):
    la = pa.read_text(errors='replace').splitlines()
    lb = pb.read_text(errors='replace').splitlines()
    notes = []
    if len(la) != len(lb):
        notes.append(f'**lines {len(la)} → {len(lb)}**')
    d = Diff()
    import re
    split = re.compile(r'[\t =,;:]+')
    for x, y in zip(la, lb):
        if x == y:
            d.values += 1
            continue
        ta, tb = split.split(x), split.split(y)
        if len(ta) != len(tb):
            d.other(x, y)
            continue
        for u, v in zip(ta, tb):
            d.value(u, v)
    cols = {'lines': d.text()} if d.text() else {}
    return notes, cols


def load_pickle(path):
    with gzip.open(path, 'rb') as f:
        return pickle.load(f)


def compare_file(pa, pb):
    name = pa.name
    try:
        if name.endswith('.pq') or name.endswith('.tsv') and not name.startswith('toppic'):
            return compare_frames(read_table(pa), read_table(pb))
        if name.endswith('.mzML'):
            return compare_mzml(pa, pb)
        if name.endswith('.pkl.gz'):
            d = Diff()
            d.value(normalize(load_pickle(pa)), normalize(load_pickle(pb)))
            return [], ({'values': d.text()} if d.text() else {})
        if name.endswith('.json'):
            d = Diff()
            d.value(json.loads(pa.read_text()), json.loads(pb.read_text()))
            return [], ({'values': d.text()} if d.text() else {})
        if pa.read_bytes() == pb.read_bytes():
            return [], {}
        return compare_lines(pa, pb)
    except Exception as error:  # report, never stop the comparison
        return [f'could not compare: {type(error).__name__}: {error}'], {}


def normalize(obj):
    """Frames and other objects as plain python values."""
    try:
        import pandas as pd
        import polars as pl
        if isinstance(obj, pl.LazyFrame):
            obj = obj.collect()
        if isinstance(obj, pl.DataFrame):
            return obj.to_dict(as_series=False)
        if isinstance(obj, pd.DataFrame):
            return {c: obj[c].tolist() for c in obj.columns}
    except ImportError:
        pass
    if isinstance(obj, dict):
        return {str(k): normalize(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [normalize(v) for v in obj]
    return obj


def compare(old, new, md):
    old, new = Path(old) / 'files', Path(new) / 'files'
    lines = ['# FLASHApp example results: old build → new build', '',
             f'Numbers that differ by less than {REL_TOL:g} (relative) are counted as last-digit differences.', '']
    total_diff = 0
    for dataset in sorted({p.name for p in old.iterdir()} | {p.name for p in new.iterdir()}):
        da, db = old / dataset, new / dataset
        lines.append(f'## {dataset}')
        if not da.exists() or not db.exists():
            lines += [f'**only in {"new" if db.exists() else "old"} build**', '']
            continue
        fa = {p.relative_to(da).as_posix() for p in da.rglob('*') if p.is_file()}
        fb = {p.relative_to(db).as_posix() for p in db.rglob('*') if p.is_file()}
        if fa - fb:
            lines.append(f'- **files only in old:** {sorted(fa - fb)}')
        if fb - fa:
            lines.append(f'- **files only in new:** {sorted(fb - fa)}')
        identical = []
        for rel in sorted(fa & fb):
            pa, pb = da / rel, db / rel
            if pa.read_bytes() == pb.read_bytes():
                identical.append(rel)
                continue
            notes, cols = compare_file(pa, pb)
            if not notes and not cols:
                identical.append(rel + ' (same values)')
                continue
            if any('**' in t for t in notes + list(cols.values())):
                total_diff += 1
            lines.append(f'- `{rel}`' + (': ' + '; '.join(notes) if notes else ''))
            for c, t in cols.items():
                lines.append(f'  - `{c}`: {t}')
        lines.append(f'- identical: {", ".join(identical) if identical else "none"}')
        lines.append('')
    lines.insert(4, f'**{total_diff} files with differences beyond the last digits.**\n')
    text = '\n'.join(lines)
    print(text)
    if md:
        Path(md).write_text(text)


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest='cmd', required=True)
    r = sub.add_parser('run')
    r.add_argument('out')
    c = sub.add_parser('compare')
    c.add_argument('old')
    c.add_argument('new')
    c.add_argument('--md')
    args = parser.parse_args()
    if args.cmd == 'run':
        run(args.out)
    else:
        compare(args.old, args.new, args.md)


if __name__ == '__main__':
    main()
