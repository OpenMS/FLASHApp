"""The comparison table shows what changed by default and keeps the rest in
the column menu (src/render/quant_cells.py::_table_columns)."""

import polars as pl
import pytest

pytest.importorskip('openms_insight')

from src.quant.session import QuantSession  # noqa: E402
from src.render.quant_cells import _table_columns  # noqa: E402


class _Session:
    sample_noun = 'sample'

    def __init__(self, level, annotations, design):
        self.level = level
        self._annotations = annotations
        self._design = design
        self.manifest = {'runs': list(design['sample'])}

    def annotations(self, level=None):
        return self._annotations.lazy()

    def design(self):
        return self._design

    hidden = []
    all_runs = QuantSession.all_runs
    runs = QuantSession.runs
    conditions = QuantSession.conditions
    compared = QuantSession.compared


DESIGN = pl.DataFrame({
    'sample': ['1x_a', '1x_b', '128N', '4x_a', '4x_b'],
    'condition': ['1x', '1x', '', '4x', '4x'],
    'replicate': [1, 2, 1, 1, 2],
})


def _visible(columns, hidden):
    return [c['title'] for c in columns if c['field'] not in hidden]


def test_feature_table_shows_what_changed():
    annotations = pl.DataFrame({'id': ['1', '2'], 'mass': [1.0, 2.0], 'rt': [1.0, 2.0],
                                'parent': [None, None], 'PrSMs': [0, 0]})
    session = _Session('feature', annotations, DESIGN)
    columns, hidden = _table_columns(session, session.compared)
    assert _visible(columns, hidden) == ['Feature', 'Change', 'log2 fold change', 'Adjusted p-value', 'Detected in']
    # Everything else is one click away in the column menu
    assert {'mass', 'rt', 'parent', 'PrSMs', '1x_a', '4x_b'} <= set(hidden)
    assert set(hidden) <= {c['field'] for c in columns}


def test_identified_features_show_the_proteoform():
    annotations = pl.DataFrame({'id': ['1', '2'], 'mass': [1.0, 2.0], 'rt': [1.0, 2.0], 'parent': ['P1', None]})
    session = _Session('feature', annotations, DESIGN)
    columns, hidden = _table_columns(session, session.compared)
    assert 'Identified as' in _visible(columns, hidden)


def test_samples_without_condition_are_not_compared():
    session = _Session('feature', pl.DataFrame({'id': ['1']}), DESIGN)
    assert session.compared == ['1x_a', '1x_b', '4x_a', '4x_b']
    columns, _ = _table_columns(session, session.compared)
    fields = [c['field'] for c in columns]
    assert '128N' not in fields
    detected = next(c for c in columns if c['field'] == 'Detected')
    assert detected['formatterParams']['symbol'] == '/4'


def test_every_run_compared_while_none_has_a_condition():
    design = DESIGN.with_columns(pl.lit('').alias('condition'))
    session = _Session('feature', pl.DataFrame({'id': ['1']}), design)
    assert session.compared == list(DESIGN['sample'])
