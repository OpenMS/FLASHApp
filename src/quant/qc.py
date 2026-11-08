"""What the identification and quality-control views of a comparison show.

Streamlit-free: every function takes the long table of a level (one row per
entity and run: ``id``, ``sample``, ``intensity``, ``status``, ``score``,
``qvalue``; src/quant/contract.py) and the runs (samples) to show, and returns
a polars frame the cells in src/render/quant_cells.py draw.

"Detected" means the run found the entity itself; "matched" that it was
measured where another run found it (match between runs), with or without an
identification of its own.
"""

import math

import polars as pl

FOUND = ('detected', 'transferred_identified', 'transferred')
STATUS_NAMES = {
    'detected': 'detected',
    'transferred_identified': 'matched, identified',
    'transferred': 'matched',
}


def identification_level(levels, level):
    """The level the identification views count: the level shown, but
    proteoforms instead of features when the runs were identified."""
    if level == 'feature' and 'proteoform' in levels:
        return 'proteoform'
    return level


def _runs(long, samples):
    return long.filter(pl.col('sample').is_in(list(samples)))


def counts_per_run(long, samples):
    """Entities found in each run, by how: ``sample``, ``status`` (named as
    in STATUS_NAMES), ``count``; one row per run and way of finding,
    in run order."""
    samples = list(samples)
    found = _runs(long, samples).filter(pl.col('status').is_in(FOUND))
    counts = found.group_by('sample', 'status').agg(pl.len().cast(pl.Int64).alias('count'))
    grid = pl.DataFrame({'sample': samples}).join(pl.DataFrame({'status': list(FOUND)}), how='cross')
    return (
        grid.join(counts, on=['sample', 'status'], how='left')
        .with_columns(pl.col('count').fill_null(0))
        .with_columns(
            pl.col('sample').replace_strict({s: i for i, s in enumerate(samples)}, return_dtype=pl.Int64).alias('_run'),
            pl.col('status').replace_strict({s: i for i, s in enumerate(FOUND)}, return_dtype=pl.Int64).alias('_how'),
        )
        .sort('_run', '_how').drop('_run', '_how')
        .with_columns(pl.col('status').replace_strict(STATUS_NAMES))
    )


def frequency(long, samples):
    """How many entities were found in 1, 2, ... n runs: ``runs``,
    ``detected`` (found by the runs themselves) and ``found`` (detected or
    matched)."""
    samples = list(samples)
    per_entity = _runs(long, samples).group_by('id').agg(
        (pl.col('status') == 'detected').sum().cast(pl.Int64).alias('detected'),
        pl.col('status').is_in(FOUND).sum().cast(pl.Int64).alias('found'),
    )
    rows = []
    for k in range(1, len(samples) + 1):
        rows.append({
            'runs': k,
            'detected': per_entity.filter(pl.col('detected') == k).height,
            'found': per_entity.filter(pl.col('found') == k).height,
        })
    return pl.DataFrame(rows, schema={'runs': pl.Int64, 'detected': pl.Int64, 'found': pl.Int64})


def overlap(long, samples):
    """Shared entities of every pair of runs, in percent of those either run
    detected (Jaccard index): ``Sample`` and one column per run."""
    samples = list(samples)
    detected = _runs(long, samples).filter(pl.col('status') == 'detected')
    sets = {s: set(detected.filter(pl.col('sample') == s)['id'].to_list()) for s in samples}
    rows = {'Sample': samples}
    for other in samples:
        column = []
        for sample in samples:
            union = sets[sample] | sets[other]
            column.append(100.0 * len(sets[sample] & sets[other]) / len(union) if union else None)
        rows[other] = column
    return pl.DataFrame(rows, schema={'Sample': pl.Utf8, **{s: pl.Float64 for s in samples}})


def identification_table(long, annotations, samples):
    """One row per entity found in any run: its annotations, in how many
    runs it was detected and found (``Detected``, ``Found``) and its status
    in every run (columns named after the runs), most often found first."""
    samples = list(samples)
    runs = _runs(long, samples)
    status = (
        runs.select('id', 'sample', pl.col('status').replace_strict(STATUS_NAMES, default='missing'))
        .pivot(on='sample', index='id', values='status', aggregate_function='first')
    )
    for sample in samples:
        if sample not in status.columns:
            status = status.with_columns(pl.lit('missing').alias(sample))
    counts = runs.group_by('id').agg(
        (pl.col('status') == 'detected').sum().cast(pl.Int64).alias('Detected'),
        pl.col('status').is_in(FOUND).sum().cast(pl.Int64).alias('Found'),
    )
    keep = [c for c in ('label', 'accession', 'description', 'proforma', 'mass', 'rt', 'PrSMs', 'IdentifiedRuns',
                        'proteoforms') if c in annotations.columns]
    table = (
        annotations.select('id', *keep).join(counts, on='id', how='inner').join(status, on='id', how='left')
        .filter(pl.col('Found') > 0)
        .sort(['Found', 'Detected', 'label'], descending=[True, True, False])
        .with_row_index('RowIndex')
        .with_columns(pl.col('RowIndex').cast(pl.Int64), pl.col(samples).fill_null('missing'))
    )
    if 'accession' in table.columns:
        table = table.with_columns(
            pl.col('accession').str.split('|').list.get(-1, null_on_oob=True).alias('protein_label'),
        )
    return table


def totals(long, samples):
    """Summed intensity of each run (sample): ``sample``, ``total``,
    ``quantified`` (entities with a value) and ``relative`` (total over the
    median total), in run order."""
    samples = list(samples)
    values = _runs(long, samples).filter(pl.col('intensity').is_not_null() & (pl.col('intensity') > 0))
    sums = values.group_by('sample').agg(
        pl.col('intensity').sum().alias('total'), pl.len().cast(pl.Int64).alias('quantified'),
    )
    frame = pl.DataFrame({'sample': samples}).join(sums, on='sample', how='left').with_columns(
        pl.col('total').fill_null(0.0), pl.col('quantified').fill_null(0),
    )
    positive = [t for t in frame['total'].to_list() if t > 0]
    median = sorted(positive)[len(positive) // 2] if positive else None
    if positive and len(positive) % 2 == 0:
        median = (sorted(positive)[len(positive) // 2 - 1] + sorted(positive)[len(positive) // 2]) / 2
    return frame.with_columns(
        (pl.col('total') / median if median else pl.lit(None, pl.Float64)).alias('relative'),
    )


def unlike_replicates(totals_frame, conditions, factor=2.0):
    """The runs whose total is under 1/``factor`` or over ``factor`` times
    the median total of their condition (conditions of two or more runs)."""
    by_condition = {}
    for row in totals_frame.iter_rows(named=True):
        condition = conditions.get(row['sample'])
        if condition and row['total'] > 0:
            by_condition.setdefault(condition, []).append((row['sample'], row['total']))
    unlike = []
    for members in by_condition.values():
        if len(members) < 2:
            continue
        values = sorted(t for _, t in members)
        n = len(values)
        median = values[n // 2] if n % 2 else (values[n // 2 - 1] + values[n // 2]) / 2
        unlike += [s for s, t in members if t < median / factor or t > median * factor]
    return unlike


def correlation(long, samples, min_shared=3):
    """Pearson correlation of the log2 intensities of every pair of runs,
    over the entities both quantified (null with fewer than ``min_shared``):
    ``Sample`` and one column per run."""
    samples = list(samples)
    values = _runs(long, samples).filter(pl.col('intensity').is_not_null() & (pl.col('intensity') > 0))
    wide = values.with_columns(pl.col('intensity').log(2)).pivot(
        on='sample', index='id', values='intensity', aggregate_function='first',
    )
    for sample in samples:
        if sample not in wide.columns:
            wide = wide.with_columns(pl.lit(None, pl.Float64).alias(sample))
    rows = {'Sample': samples}
    for other in samples:
        column = []
        for sample in samples:
            pair = wide.select(sample, other).drop_nulls() if sample != other else wide.select(sample).drop_nulls()
            if pair.height < min_shared:
                column.append(None)
            elif sample == other:
                column.append(1.0)
            else:
                r = pair.select(pl.corr(sample, other)).item()
                column.append(None if r is None or (isinstance(r, float) and math.isnan(r)) else float(r))
        rows[other] = column
    return pl.DataFrame(rows, schema={'Sample': pl.Utf8, **{s: pl.Float64 for s in samples}})


def coefficients_of_variation(wide, conditions):
    """CV (in percent) of every entity within each condition, from the raw
    intensities of ``wide`` (one column per sample): ``condition``, ``cv``.
    Only entities with values in at least two samples of a condition."""
    by_condition = {}
    for sample, condition in conditions.items():
        if sample in wide.columns:
            by_condition.setdefault(condition, []).append(sample)
    parts = []
    for condition, columns in sorted(by_condition.items()):
        if len(columns) < 2:
            continue
        values = wide.select([
            pl.when(pl.col(c).cast(pl.Float64) > 0).then(pl.col(c).cast(pl.Float64)).alias(c) for c in columns
        ])
        stats = values.select(
            pl.concat_list(columns).list.drop_nulls().alias('v'),
        ).with_columns(
            pl.col('v').list.len().alias('n'),
            pl.col('v').list.mean().alias('mean'),
            pl.col('v').list.std().alias('sd'),
        ).filter((pl.col('n') >= 2) & (pl.col('mean') > 0))
        parts.append(stats.select(
            pl.lit(condition).alias('condition'), (100.0 * pl.col('sd') / pl.col('mean')).alias('cv'),
        ))
    if not parts:
        return pl.DataFrame(schema={'condition': pl.Utf8, 'cv': pl.Float64})
    return pl.concat(parts)


def match_quality(long, samples):
    """The values matched between runs (match between runs) and the matches
    that were rejected: ``sample``, ``outcome`` (``matched, identified``,
    ``matched`` or ``not matched``), ``score`` and ``qvalue``; values a run
    detected itself have no match score."""
    samples = list(samples)
    rows = _runs(long, samples).filter((pl.col('status') != 'detected') & pl.col('score').is_not_null())
    return rows.select(
        'sample',
        pl.col('status').replace_strict({**STATUS_NAMES, 'missing': 'not matched'}, default='not matched')
        .alias('outcome'),
        pl.col('score').cast(pl.Float64), pl.col('qvalue').cast(pl.Float64),
    )
