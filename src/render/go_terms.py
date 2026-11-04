"""GO term enrichment of the proteins selected in the Viewer.

The user selects several proteins or proteoforms at once: a box or lasso on
the volcano plot or a heatmap, the volcano plot's All changed / Up / Down
buttons, or the check boxes of a protein, proteoform or comparison table. The
components write the selection into one identifier of the grid's
StateManager (``RUN_SET`` in a run, holding UniProt accessions; ``COMPARISON_SET``
in a comparison, holding the entries' labels), and the GO enrichment cell
tests which Gene Ontology terms are over-represented in it against all
proteins of the run or comparison (Fisher's exact test, Benjamini-Hochberg),
with openms_insight.analysis.enrichment.

GO annotations come from MyGene.info by UniProt accession and are kept in the
workspace cache, so each protein is looked up once. Everything that does not
depend on the selection is prepared once and stored in the cache: the
workflow looks up the GO terms of every identified protein after a run and
writes the run's term index (``precompute_run``), a table with one row per
protein and term. A selection then only counts its proteins in that index and
runs the tests (``enrich``). Results of earlier builds get their index the
first time the cell shows them. Free of Streamlit.
"""

import hashlib
import json
import threading
import time
from collections import OrderedDict
from pathlib import Path

import numpy as np
import polars as pl

# Identifiers of the selected set
RUN_SET = 'go_proteins'
COMPARISON_SET = 'go_entries'

CATEGORIES = {'BP': 'Biological process', 'CC': 'Cellular component', 'MF': 'Molecular function'}

# A protein MyGene.info did not know is asked again after this long
RETRY_NOT_FOUND_S = 30 * 24 * 3600

_CACHE_LOCK = threading.Lock()

# Bump when the index files change shape
INDEX_VERSION = 1
# Term indexes kept in memory per process
_INDEXES_KEPT = 16
_INDEXES = OrderedDict()
_INDEX_LOCK = threading.Lock()


def clean_accession(name):
    """'sp|P12345|NAME_HUMAN' -> 'P12345'; a bare accession unchanged."""
    if name is None:
        return None
    parts = str(name).split('|')
    return parts[1] if len(parts) >= 2 else parts[0]


def protein_label(name):
    """'sp|P12345|NAME_HUMAN' -> 'NAME_HUMAN' (the entry name), else the name."""
    if name is None:
        return None
    parts = str(name).split('|')
    return parts[-1] if len(parts) >= 3 and parts[-1] else str(name)


# -- backgrounds and selections ------------------------------------------------

def run_background(file_manager, dataset_id):
    """Accessions of the proteins identified in a run (the PrSMs counted in
    the protein, proteoform and PrSM tables)."""
    from src.render import hierarchy_data
    hits = hierarchy_data.prsms(file_manager, dataset_id)
    if not hits.height:
        return []
    return sorted(a for a in hits["accession"].drop_nulls().unique().to_list() if a)


def comparison_entries(session):
    """The entries of a comparison's level with the protein each belongs to:
    columns Proteoform (the label the components select by) and accession
    (null for a feature without an identification)."""
    level = session.level
    table = session.table(level)
    if level == 'protein':
        frame = table.select('Proteoform', pl.col('id').alias('accession'))
    elif level == 'proteoform':
        frame = table.select('Proteoform', 'accession')
    else:
        # A feature: the protein of the proteoform its PrSMs identify
        if 'parent' in table.columns and 'proteoform' in session.levels:
            forms = session.annotations('proteoform').collect().select(
                pl.col('id').alias('parent'), 'accession',
            )
            frame = table.select('Proteoform', 'parent').join(forms, on='parent', how='left').select(
                'Proteoform', 'accession',
            )
        else:
            frame = table.select('Proteoform', pl.lit(None, pl.Utf8).alias('accession'))
    return frame.unique(maintain_order=True)


def comparison_universe(session):
    """Accessions of every entry of a comparison's level: the term index of
    a comparison is built for these, so the tested background and any
    selection come from the same index."""
    return sorted({a for a in comparison_entries(session)['accession'] if a})


def comparison_selection(session, labels):
    """(selected accessions, background accessions, selected entries, entries
    of the selection without a protein) for the labels selected in a
    comparison. The background is every tested entry's protein."""
    entries = comparison_entries(session)
    tested = set(
        session.result('statistics').drop_nulls(['log2FC', 'p-adj']).select('Proteoform').collect()['Proteoform']
    )
    background = sorted({a for a in entries.filter(pl.col('Proteoform').is_in(list(tested)))['accession'] if a})
    chosen = entries.filter(pl.col('Proteoform').is_in(list(labels or [])))
    selected = sorted({a for a in chosen['accession'] if a})
    without = chosen.filter(pl.col('accession').is_null()).height
    return selected, background, len(set(labels or [])), without


# -- annotation cache ------------------------------------------------------------

def _cache_file(file_manager):
    return Path(file_manager.cache_path, 'go', 'annotations.json')


def _read_cache(path):
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return {}


def annotations(file_manager, accessions, fetch=None, now=None):
    """GO terms of the accessions, as the frame openms_insight's
    go_enrichment takes (query, BP_terms, CC_terms, MF_terms), from the
    workspace cache, looking up the ones it does not hold.

    `fetch(clean_accessions)` returns a frame like fetch_go_annotations (the
    default; a network call to MyGene.info) and may raise; the cache keeps
    what it got. Returns (frame, error message or None).
    """
    import pandas as pd

    now = time.time() if now is None else now
    wanted = sorted({clean_accession(a) for a in accessions if a})
    path = _cache_file(file_manager)
    error = None
    with _CACHE_LOCK:
        cache = _read_cache(path)
        missing = [
            a for a in wanted
            if a not in cache or (
                cache[a].get('found') is False and now - cache[a].get('at', 0) > RETRY_NOT_FOUND_S
            )
        ]
        if missing:
            if fetch is None:
                from openms_insight.analysis.enrichment import fetch_go_annotations as fetch
            try:
                fetched = fetch(missing)
            except Exception as e:  # network or service error: keep what we have
                fetched = None
                error = f'GO terms could not be fetched from MyGene.info ({type(e).__name__}).'
            if fetched is not None:
                found = {}
                for _, row in fetched.iterrows():
                    found[str(row['query'])] = {t: list(row[f'{t}_terms']) for t in CATEGORIES}
                for accession in missing:
                    terms = found.get(accession)
                    cache[accession] = (
                        {'found': True, 'at': now, **terms} if terms is not None
                        else {'found': False, 'at': now}
                    )
                path.parent.mkdir(parents=True, exist_ok=True)
                tmp = path.with_suffix('.tmp')
                tmp.write_text(json.dumps(cache))
                tmp.replace(path)
    rows = [
        {'query': a, **{f'{t}_terms': cache[a].get(t, []) for t in CATEGORIES}}
        for a in wanted if a in cache and cache[a].get('found')
    ]
    frame = pd.DataFrame(rows, columns=['query'] + [f'{t}_terms' for t in CATEGORIES])
    return frame, error


# -- term index -------------------------------------------------------------------

def _digest(accessions):
    return hashlib.md5(json.dumps(sorted(accessions)).encode()).hexdigest()[:16]


def _index_file(file_manager, universe):
    return Path(file_manager.cache_path, 'go', 'index', f'v{INDEX_VERSION}-{_digest(universe)}.parquet')


def index_frame(annotation_frame, universe):
    """The term index of the proteins `universe` (clean accessions): columns
    accession, category, term, one row per protein and term of a category; a
    protein without any term has one row with null category and term, so
    the index also lists every protein it was built for."""
    rows_acc, rows_cat, rows_term = [], [], []
    seen = set()
    for _, row in annotation_frame.iterrows():
        accession = str(row['query'])
        if accession not in universe:
            continue
        for category in CATEGORIES:
            for term in row[f'{category}_terms']:
                rows_acc.append(accession)
                rows_cat.append(category)
                rows_term.append(str(term))
                seen.add(accession)
    for accession in sorted(set(universe) - seen):
        rows_acc.append(accession)
        rows_cat.append(None)
        rows_term.append(None)
    return pl.DataFrame(
        {'accession': rows_acc, 'category': rows_cat, 'term': rows_term},
        schema={'accession': pl.Utf8, 'category': pl.Utf8, 'term': pl.Utf8},
    ).unique(maintain_order=True)


def index_ready(file_manager, accessions):
    """Whether the term index of these proteins is stored."""
    universe = {clean_accession(a) for a in accessions if a}
    return _index_file(file_manager, universe).exists()


def term_index(file_manager, accessions, fetch=None):
    """The term index (see index_frame) of the proteins `accessions`, from
    memory or the workspace cache, else built from their annotations (which
    may ask MyGene.info for the ones not looked up yet). Returns (index,
    error message or None). An index built while MyGene.info could not be
    reached is not stored, so the next call asks again."""
    universe = {clean_accession(a) for a in accessions if a}
    path = _index_file(file_manager, universe)
    with _INDEX_LOCK:
        cached = _INDEXES.get(str(path))
        if cached is not None:
            _INDEXES.move_to_end(str(path))
            return cached, None
    error = None
    try:
        index = pl.read_parquet(path)
    except (OSError, pl.exceptions.PolarsError):
        frame, error = annotations(file_manager, sorted(universe), fetch=fetch)
        index = index_frame(frame, universe)
        if error is None:
            path.parent.mkdir(parents=True, exist_ok=True)
            tmp = path.with_suffix(f'.{time.time_ns()}.tmp')
            index.write_parquet(tmp)
            tmp.replace(path)
    if error is None:
        with _INDEX_LOCK:
            _INDEXES[str(path)] = index
            while len(_INDEXES) > _INDEXES_KEPT:
                _INDEXES.popitem(last=False)
    return index, error


def _run_file(file_manager, dataset_id):
    return Path(file_manager.cache_path, 'go', 'runs', f'{dataset_id}.json')


def _run_stamp(file_manager, dataset_id):
    from src.render import insight_data
    path = insight_data.result_path(file_manager, dataset_id, 'protein_dfs')
    if path is None:
        return None
    stat = path.stat()
    return [path.name, stat.st_size, stat.st_mtime_ns]


def stored_run_background(file_manager, dataset_id):
    """run_background, kept in the workspace cache next to the run's term
    index while the run's protein table is unchanged."""
    path = _run_file(file_manager, dataset_id)
    stamp = _run_stamp(file_manager, dataset_id)
    try:
        stored = json.loads(path.read_text())
        if stored.get('stamp') == stamp and stored.get('version') == INDEX_VERSION:
            return stored['background']
    except (OSError, ValueError, KeyError):
        pass
    background = run_background(file_manager, dataset_id) if stamp is not None else []
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(f'.{time.time_ns()}.tmp')
    tmp.write_text(json.dumps({'version': INDEX_VERSION, 'stamp': stamp, 'background': background}))
    tmp.replace(path)
    return background


def precompute_run(file_manager, dataset_id, fetch=None):
    """Prepare the GO enrichment cell of a run: its background, the GO terms
    of its proteins and its term index. Returns (number of proteins, error
    message or None). Called by the workflow after a run."""
    background = stored_run_background(file_manager, dataset_id)
    if not background:
        return 0, None
    _, error = term_index(file_manager, background, fetch=fetch)
    return len(background), error


# -- enrichment ---------------------------------------------------------------------

_TABLE_COLUMNS = [
    'GO_Term', 'Count', 'Background', 'GeneRatio', 'BgRatio', 'FoldEnrichment',
    'p_value', '-log10(p)', 'p_adj', 'Proteins',
]


def _category_table(terms, fg, n_fg, n_bg):
    """Fisher's exact test (one-sided) of every term of one category that a
    selected protein has, as openms_insight's go_term_table writes it."""
    import pandas as pd
    from scipy.stats import hypergeom
    from openms_insight.analysis.enrichment import benjamini_hochberg

    picked = terms.filter(pl.col('accession').is_in(list(fg)))
    if not picked.height:
        return pd.DataFrame(columns=_TABLE_COLUMNS)
    counts = (
        picked.group_by('term').agg(
            pl.col('accession').n_unique().alias('Count'),
            pl.col('accession').unique().sort().str.join(';').alias('Proteins'),
        )
        .join(terms.group_by('term').agg(pl.col('accession').n_unique().alias('Background')), on='term')
        .sort('term')
    )
    a = counts['Count'].to_numpy().astype(np.int64)
    k = counts['Background'].to_numpy().astype(np.int64)
    # P(X >= a) for X ~ hypergeometric: the one-sided Fisher's exact test
    # (many terms share their counts, so each pair is computed once)
    pairs, inverse = np.unique(np.stack([a, k]), axis=1, return_inverse=True)
    p = np.clip(hypergeom.sf(pairs[0] - 1, n_bg, n_fg, pairs[1]), 0.0, 1.0)[np.ravel(inverse)]
    frame = pd.DataFrame({
        'GO_Term': counts['term'].to_list(),
        'Count': a,
        'Background': k,
        'GeneRatio': [f'{x}/{n_fg}' for x in a],
        'BgRatio': [f'{x}/{n_bg}' for x in k],
        'FoldEnrichment': (a / n_fg) / (k / n_bg),
        'p_value': p,
        '-log10(p)': -np.log10(np.clip(p, 1e-300, None)),
        'p_adj': benjamini_hochberg(p),
        'Proteins': counts['Proteins'].to_list(),
    })
    return frame.sort_values(['p_value', 'GO_Term']).reset_index(drop=True)[_TABLE_COLUMNS]


def enrich(index, foreground, background, categories=tuple(CATEGORIES)):
    """GO enrichment of the proteins `foreground` against `background`
    (accessions) from a term index, with the result of openms_insight's
    go_enrichment: fg_total, bg_total, fg_count, bg_count (proteins with any
    term) and categories (category -> table of go_term_table)."""
    fg = {clean_accession(a) for a in foreground if a}
    bg = {clean_accession(a) for a in background if a} | fg
    terms = index.filter(
        pl.col('accession').is_in(list(bg)) & pl.col('category').is_in(list(categories))
    )
    annotated = set(terms['accession'].unique().to_list())
    fg_annotated = fg & annotated
    n_fg, n_bg = len(fg_annotated), len(annotated)
    return {
        'fg_total': len(fg),
        'bg_total': len(bg),
        'fg_count': n_fg,
        'bg_count': n_bg,
        'categories': {
            t: _category_table(terms.filter(pl.col('category') == t), fg_annotated, n_fg, n_bg)
            for t in categories
        },
    }


def enrichment(file_manager, selected, background, fetch=None, universe=None):
    """GO enrichment of `selected` against `background` (accessions), from
    the term index of `universe` (default: the background; it has to hold
    both). Returns (result of enrich, error message or None, accession ->
    label of the proteins)."""
    pool = list(universe if universe is not None else background) + list(selected)
    index, error = term_index(file_manager, pool, fetch=fetch)
    result = enrich(index, selected, background)
    labels = {clean_accession(a): protein_label(a) for a in list(selected) + list(pool) if a}
    return result, error, labels


def term_frame(table, labels):
    """The terms of one category for the Viewer's table: one row per term,
    the proteins by their entry names."""
    if table is None or len(table) == 0:
        return pl.DataFrame(schema={
            'rank': pl.Int64, 'term': pl.Utf8, 'count': pl.Int64, 'background': pl.Int64,
            'ratio': pl.Utf8, 'bg_ratio': pl.Utf8, 'fold': pl.Float64, 'pvalue': pl.Float64,
            'padj': pl.Float64, 'proteins': pl.Utf8,
        })
    proteins = [
        ', '.join(labels.get(a, a) for a in str(p).split(';') if a) for p in table['Proteins']
    ]
    return pl.DataFrame({
        'rank': list(range(1, len(table) + 1)),
        'term': table['GO_Term'].astype(str).tolist(),
        'count': table['Count'].astype(int).tolist(),
        'background': table['Background'].astype(int).tolist(),
        'ratio': table['GeneRatio'].astype(str).tolist(),
        'bg_ratio': table['BgRatio'].astype(str).tolist(),
        'fold': table['FoldEnrichment'].astype(float).tolist(),
        'pvalue': table['p_value'].astype(float).tolist(),
        'padj': table['p_adj'].astype(float).tolist(),
        'proteins': proteins,
    })
