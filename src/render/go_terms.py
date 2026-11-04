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
workspace cache, so each protein is looked up once. Free of Streamlit.
"""

import json
import threading
import time
from pathlib import Path

import polars as pl

# Identifiers of the selected set
RUN_SET = 'go_proteins'
COMPARISON_SET = 'go_entries'

CATEGORIES = {'BP': 'Biological process', 'CC': 'Cellular component', 'MF': 'Molecular function'}

# A protein MyGene.info did not know is asked again after this long
RETRY_NOT_FOUND_S = 30 * 24 * 3600

_CACHE_LOCK = threading.Lock()


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


def enrichment(file_manager, selected, background, fetch=None):
    """GO enrichment of `selected` against `background` (accessions).
    Returns (result of go_enrichment, error message or None, accession ->
    label of the proteins)."""
    from openms_insight.analysis.enrichment import go_enrichment
    frame, error = annotations(file_manager, list(selected) + list(background), fetch=fetch)
    result = go_enrichment(frame, selected, background, categories=tuple(CATEGORIES))
    labels = {clean_accession(a): protein_label(a) for a in list(selected) + list(background) if a}
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
