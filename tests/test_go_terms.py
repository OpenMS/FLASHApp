"""GO enrichment of the proteins selected in the Viewer (src/render/go_terms.py)."""

import types

import pandas as pd
import polars as pl
import pytest

pytest.importorskip('openms_insight.analysis.enrichment')

from src.render import go_terms, hierarchy_data  # noqa: E402


TERMS = {
    'P1': {'BP': ['glycolysis', 'transport'], 'CC': ['cytosol'], 'MF': []},
    'P2': {'BP': ['glycolysis'], 'CC': ['cytosol'], 'MF': ['kinase']},
    'P3': {'BP': ['transport'], 'CC': ['membrane'], 'MF': []},
    'P4': {'BP': ['translation'], 'CC': ['ribosome'], 'MF': ['rRNA binding']},
}


class FakeFetch:
    def __init__(self, fail=False):
        self.calls = []
        self.fail = fail

    def __call__(self, accessions):
        self.calls.append(list(accessions))
        if self.fail:
            raise ConnectionError('offline')
        rows = [
            {'query': a, **{f'{t}_terms': TERMS[a][t] for t in ('BP', 'CC', 'MF')}}
            for a in accessions if a in TERMS
        ]
        return pd.DataFrame(rows, columns=['query', 'BP_terms', 'CC_terms', 'MF_terms'])


def _fm(tmp_path):
    return types.SimpleNamespace(cache_path=str(tmp_path))


def test_accession_helpers():
    assert go_terms.clean_accession('sp|P12345|ALDOA_HUMAN') == 'P12345'
    assert go_terms.clean_accession('P12345') == 'P12345'
    assert go_terms.protein_label('sp|P12345|ALDOA_HUMAN') == 'ALDOA_HUMAN'
    assert go_terms.protein_label('P12345') == 'P12345'


def test_annotations_are_cached_and_not_found_kept(tmp_path):
    fm = _fm(tmp_path)
    fetch = FakeFetch()
    frame, error = go_terms.annotations(fm, ['sp|P1|A_HUMAN', 'P2', 'Q9'], fetch=fetch, now=0)
    assert error is None
    assert sorted(frame['query']) == ['P1', 'P2']
    assert fetch.calls == [['P1', 'P2', 'Q9']]
    # Looked up once: the cache answers, also for the unknown Q9
    go_terms.annotations(fm, ['P1', 'P2', 'Q9'], fetch=fetch, now=10)
    assert len(fetch.calls) == 1
    # An unknown protein is asked again after a month
    go_terms.annotations(fm, ['Q9'], fetch=fetch, now=go_terms.RETRY_NOT_FOUND_S + 1)
    assert fetch.calls[-1] == ['Q9']


def test_annotations_offline_reports_and_keeps_cache(tmp_path):
    fm = _fm(tmp_path)
    go_terms.annotations(fm, ['P1'], fetch=FakeFetch(), now=0)
    frame, error = go_terms.annotations(fm, ['P1', 'P2'], fetch=FakeFetch(fail=True), now=1)
    assert 'MyGene.info' in error
    assert frame['query'].tolist() == ['P1']


def test_enrichment_of_a_selection(tmp_path):
    fm = _fm(tmp_path)
    selected = ['sp|P1|A_HUMAN', 'sp|P2|B_HUMAN']
    background = selected + ['sp|P3|C_HUMAN', 'sp|P4|D_HUMAN']
    result, error, labels = go_terms.enrichment(fm, selected, background, fetch=FakeFetch())
    assert error is None
    assert (result['fg_count'], result['bg_count']) == (2, 4)
    bp = result['categories']['BP']
    top = bp.iloc[0]
    assert top['GO_Term'] == 'glycolysis' and top['Count'] == 2
    frame = go_terms.term_frame(bp, labels)
    assert frame['rank'].to_list()[0] == 1
    assert frame.filter(pl.col('term') == 'glycolysis')['proteins'].item() == 'A_HUMAN, B_HUMAN'
    assert go_terms.term_frame(None, labels).height == 0


def test_run_background_skips_empty_accessions(monkeypatch):
    hits = pl.DataFrame({'accession': ['P2', 'P1', '', None, 'P1']})
    monkeypatch.setattr(hierarchy_data, 'prsms', lambda fm, ds: hits)
    assert go_terms.run_background(None, 'ds') == ['P1', 'P2']


class FakeSession:
    """A comparison at one level, with the tables comparison_selection reads."""

    def __init__(self, level):
        self.level = level
        self.levels = ['feature', 'proteoform', 'protein']
        self._stats = pl.DataFrame({
            'Proteoform': ['a', 'b', 'c', 'd'],
            'log2FC': [1.0, -2.0, 0.1, None],
            'p-adj': [0.01, 0.02, 0.5, None],
        })

    def table(self, level=None):
        base = self._stats.select('Proteoform')
        if self.level == 'protein':
            return base.with_columns(pl.Series('id', ['P1', 'P2', 'P3', 'P4']))
        if self.level == 'proteoform':
            return base.with_columns(pl.Series('accession', ['P1', 'P1', 'P3', 'P4']))
        return base.with_columns(pl.Series('parent', ['f1', 'f2', None, 'f3']))

    def annotations(self, level=None):
        return pl.LazyFrame({'id': ['f1', 'f2', 'f3'], 'accession': ['P1', 'P2', 'P4']})

    def result(self, name, level=None):
        return self._stats.lazy()


@pytest.mark.parametrize('level, picked, selected, background, without', [
    ('protein', ['a', 'b'], ['P1', 'P2'], ['P1', 'P2', 'P3'], 0),
    ('proteoform', ['a', 'b'], ['P1'], ['P1', 'P3'], 0),
    ('feature', ['a', 'c'], ['P1'], ['P1', 'P2'], 1),
])
def test_comparison_selection_maps_entries_to_proteins(level, picked, selected, background, without):
    got = go_terms.comparison_selection(FakeSession(level), picked)
    assert got == (selected, background, len(picked), without)


def test_comparison_selection_nothing_picked():
    assert go_terms.comparison_selection(FakeSession('protein'), None)[2:] == (0, 0)
