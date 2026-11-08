"""The AI assistant's view of FLASHApp's results (Streamlit-free).

FLASHApp keeps every result in the workspace's cache (``cache/cache.db`` and
``cache/files/<dataset>/``), not in run folders. ``CacheRuns`` lists those
datasets as the assistant's runs, so it can read their tables and show them
side by side in a view of the Viewer.
"""

from __future__ import annotations

import sqlite3
from datetime import datetime
from pathlib import Path

from src.render.layout_editor.catalog import (
    DECONV_DEFAULT_LAYOUT, ISOBARIC_COMPONENTS, TNT_DEFAULT_LAYOUTS, TOPDOWN_COMPONENTS,
    TOPDOWN_SEQUENCE_COMPONENTS,
)

# What the Viewer offers (the sequence components once a sequence is set)
COMPONENTS = TOPDOWN_COMPONENTS + TOPDOWN_SEQUENCE_COMPONENTS + ISOBARIC_COMPONENTS
# The ready-made views, as topdown_layout.presets names them
PRESETS = [
    {'name': 'deconvolution', 'label': 'Deconvolution', 'layout': DECONV_DEFAULT_LAYOUT},
    {'name': 'native', 'label': 'Native', 'layout': TNT_DEFAULT_LAYOUTS['Native']},
    {'name': 'lysate', 'label': 'Lysate', 'layout': TNT_DEFAULT_LAYOUTS['Lysate']},
]
# The Viewer's views file, relative to the workflow folder
VIEWS = '../layouts/viewer.json'

# Files of a run that say what it holds
_DECONVOLVED = ('spec1_tsv', 'spec2_tsv', 'out_deconv_mzML')
_IDENTIFIED = ('prsms_tsv', 'protein_tsv')


class CacheRuns:
    """The datasets of a workspace's cache, with the interface of RunManager
    the assistant uses: ``list_runs``, ``get``, ``latest``, ``results_dir``."""

    def __init__(self, workflow_dir):
        self.cache = Path(workflow_dir).parent / 'cache'

    def _rows(self) -> list[dict]:
        db = self.cache / 'cache.db'
        if not db.is_file():
            return []
        con = sqlite3.connect(f'file:{db}?mode=ro', uri=True)
        try:
            con.row_factory = sqlite3.Row
            tables = {r[0] for r in con.execute("select name from sqlite_master where type='table'")}
            rows = {}
            for table in ('stored_files', 'stored_data'):
                if table in tables:
                    for r in con.execute(f'select * from {table}'):
                        rows.setdefault(r['id'], {}).update(
                            {k: r[k] for k in r.keys() if r[k] is not None})
            return list(rows.values())
        finally:
            con.close()

    def list_runs(self) -> list[dict]:
        runs = []
        for row in self._rows():
            folder = self.cache / 'files' / row['id']
            if any(k in row for k in _DECONVOLVED):
                kind = 'identified' if any(k in row for k in _IDENTIFIED) else 'deconvolved'
            elif 'quant_manifest' in row:
                kind = 'comparison'
            else:
                continue
            created = folder.stat().st_mtime if folder.is_dir() else 0
            runs.append({
                'id': row['id'], 'name': row.get('display_name') or row['id'], 'status': 'finished',
                'kind': kind, 'isobaric': 'isobaric_dfs' in row,
                'created': datetime.fromtimestamp(created).isoformat(timespec='seconds') if created else None,
                '_mtime': created,
            })
        runs.sort(key=lambda r: r['_mtime'])
        for r in runs:
            del r['_mtime']
        return runs

    def get(self, run_id) -> dict | None:
        return next((r for r in self.list_runs() if r['id'] == run_id), None)

    def latest(self) -> str | None:
        runs = [r for r in self.list_runs() if r['kind'] != 'comparison']
        return runs[-1]['id'] if runs else None

    def results_dir(self, run_id) -> Path:
        return self.cache / 'files' / run_id
