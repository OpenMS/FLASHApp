"""The quantification of several runs, in the shared quantification format.

FLASHQuantConsensus (-out) gives one row per consensus feature with an
abundance, a status, a match score and a q-value per run; with -out_id it also
lists the PrSMs (FLASHTnT) that identify each consensus feature. This module
turns those two tables into the shared format the quantification views read:

- one long table per level (``feature``, ``proteoform``, ``protein``) with one
  row per entity and run: ``id``, ``sample``, ``intensity`` (linear; null when
  not quantified), ``status``, ``score``, ``qvalue`` and ``parent`` (the id at
  the next level up; null when there is none);
- one annotation table per level (``id``, ``label`` and level-specific
  columns such as mass, retention time, accession or ProForma);
- the experimental design (run, condition, replicate);
- a manifest (``quant.json``) with the contract version, the levels and their
  roll-ups.

Levels: a consensus feature is a FLASHDeconv feature measured in every run. A
proteoform is the ProForma (or, without one, the accession and mass) of the
PrSMs assigned to a feature; a feature whose PrSMs name two proteoforms is
ambiguous and has no parent. A proteoform's value in a run is the sum of its
features' values, a protein's the sum of its proteoforms'. Detected beats
transferred: a summed value is "detected" when any of its parts was detected
in that run.
"""

import json
import re

import polars as pl

CONTRACT_VERSION = 1
LEVELS = ('feature', 'proteoform', 'protein')

# FLASHQuantConsensus status -> contract status, best first
STATUS = {
    'detected': 'detected',
    'mbr_id': 'transferred_identified',
    'mbr': 'transferred',
    'nan': 'missing',
}
STATUS_ORDER = ['detected', 'transferred_identified', 'transferred', 'missing']

LONG_SCHEMA = {
    'id': pl.Utf8, 'sample': pl.Utf8, 'intensity': pl.Float64, 'status': pl.Utf8,
    'score': pl.Float64, 'qvalue': pl.Float64, 'parent': pl.Utf8,
}


def read_consensus(path):
    """FLASHQuantConsensus -out as (features, long): one row per consensus
    feature (ConsensusIndex, MonoisotopicMass, ApexRetentionTime, MinCharge,
    MaxCharge, DetectedRuns) and one row per feature and run (ConsensusIndex,
    sample, intensity, status, score, qvalue)."""
    table = pl.read_csv(
        path, separator='\t', infer_schema_length=0, null_values=['', 'nan', 'NaN'],
    )
    fixed = ['ConsensusIndex', 'MonoisotopicMass', 'ApexRetentionTime', 'MinCharge', 'MaxCharge', 'DetectedRuns']
    features = table.select(fixed).with_columns(
        pl.col('ConsensusIndex').cast(pl.Int64),
        pl.col('MonoisotopicMass', 'ApexRetentionTime').cast(pl.Float64),
        pl.col('MinCharge', 'MaxCharge', 'DetectedRuns').cast(pl.Int64),
    )
    runs = [c[len('Abundance_'):] for c in table.columns if c.startswith('Abundance_')]
    parts = []
    for run in runs:
        parts.append(table.select(
            pl.col('ConsensusIndex').cast(pl.Int64),
            pl.lit(run).alias('sample'),
            pl.col(f'Abundance_{run}').cast(pl.Float64).alias('intensity'),
            pl.col(f'Status_{run}').fill_null('nan').replace_strict(STATUS, default='missing').alias('status'),
            pl.col(f'MatchScore_{run}').cast(pl.Float64).alias('score'),
            pl.col(f'MatchQvalue_{run}').cast(pl.Float64).alias('qvalue'),
        ))
    long = pl.concat(parts) if parts else pl.DataFrame(
        schema={'ConsensusIndex': pl.Int64, 'sample': pl.Utf8, 'intensity': pl.Float64,
                'status': pl.Utf8, 'score': pl.Float64, 'qvalue': pl.Float64})
    # a value without a number is not quantified, whatever the status says
    long = long.with_columns(
        pl.when(pl.col('intensity').is_nan() | (pl.col('intensity') <= 0))
        .then(None).otherwise(pl.col('intensity')).alias('intensity'),
    )
    return features, long, runs


def read_identifications(path):
    """FLASHQuantConsensus -out_id: one row per PrSM assigned to a consensus
    feature, with the proteoform key used across runs."""
    table = pl.read_csv(path, separator='\t', infer_schema_length=0, null_values=[''])
    if table.is_empty():
        return pl.DataFrame(schema={
            'ConsensusIndex': pl.Int64, 'Run': pl.Utf8, 'ProteinAccession': pl.Utf8,
            'ProteinDescription': pl.Utf8, 'Proforma': pl.Utf8, 'ProteoformMass': pl.Float64,
            'ProteoformLevelQvalue': pl.Float64, 'Scan': pl.Utf8, 'proteoform': pl.Utf8,
        })
    for column in ('Proforma', 'ProteinDescription'):
        if column not in table.columns:
            table = table.with_columns(pl.lit(None, pl.Utf8).alias(column))
    table = table.with_columns(
        pl.col('ConsensusIndex').cast(pl.Int64),
        pl.col('ProteoformMass').cast(pl.Float64),
        pl.col('ProteoformLevelQvalue').cast(pl.Float64),
    )
    # ProteoformIndex is per run; the ProForma names a proteoform in every run
    return table.with_columns(
        pl.when(pl.col('Proforma').is_not_null())
        .then(pl.col('ProteinAccession') + '|' + pl.col('Proforma'))
        .otherwise(pl.col('ProteinAccession') + '|' + pl.col('ProteoformMass').round(1).cast(pl.Utf8) + ' Da')
        .alias('proteoform'),
    )



def read_sequences(path):
    """accession -> sequence of a FASTA file (empty without one)."""
    sequences, accession, parts = {}, None, []
    if not path:
        return sequences
    with open(path, errors='replace') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if accession is not None:
                    sequences[accession] = ''.join(parts)
                accession, parts = line[1:].split()[0] if len(line) > 1 else '', []
            elif accession is not None:
                parts.append(line)
    if accession is not None:
        sequences[accession] = ''.join(parts)
    return sequences


_TOKEN = re.compile(r'\(([A-Z]+)\)\[([^\]]+)\]|([A-Z])\[([^\]]+)\]|\[([^\]]+)\]-|([A-Z])')


def parse_proforma(proforma):
    """The bare sequence of a FLASHTnT ProForma and its mass shifts, each as
    (shift, first, last) with 1-based positions in that sequence: a localized
    shift has first == last, an unlocalized one (``(SEQ)[+mass]``) spans its
    residues."""
    sequence, shifts = [], []
    for block, block_shift, residue, residue_shift, nterm, plain in _TOKEN.findall(proforma or ''):
        if block:
            shifts.append((block_shift, len(sequence) + 1, len(sequence) + len(block)))
            sequence.append(block)
        elif residue:
            sequence.append(residue)
            shifts.append((residue_shift, len(sequence), len(sequence)))
        elif nterm:
            shifts.append((nterm, 1, 1))
        else:
            sequence.append(plain)
    return ''.join(sequence), shifts


def _shift_text(shift):
    """'+15.99' / '−83.18' for a mass, the name otherwise (Acetyl)."""
    try:
        value = float(shift)
    except ValueError:
        return shift
    return f'{value:+.2f}'.replace('-', '\u2212')


def describe_proteoform(accession, proforma, mass, sequences):
    """Short name and details of a proteoform: {label, entry, start, end,
    shifts}. The label is the protein's entry name, the residues the
    proteoform covers in it (when its sequence is in `sequences`) and its mass
    shifts, e.g. 'RS7_ECOLI 2\u2013157 (+446.14)'. Without the protein's
    sequence the start of the proteoform's stands for the residues; without
    a ProForma, the mass."""
    entry = (accession or '').split('|')[-1] or (accession or '')
    if not proforma:
        text = f'{mass:.1f} Da' if mass is not None and mass > 0 else 'mass not determined'
        return {'label': f'{entry} {text}', 'entry': entry, 'start': None, 'end': None, 'shifts': None}
    bare, shifts = parse_proforma(proforma)
    protein = sequences.get(accession, '')
    found = protein.find(bare) if protein and bare else -1
    start = found + 1 if found >= 0 else None
    end = found + len(bare) if found >= 0 else None
    offset = (start or 1) - 1
    details = []
    for shift, first, last in shifts:
        where = f'{first + offset}' if first == last else f'{first + offset}\u2013{last + offset}'
        details.append(f'{_shift_text(shift)} at {where}' if first == last else f'{_shift_text(shift)} in {where}')
    label = entry
    if start is not None:
        label += f' {start}\u2013{end}'
    else:
        # Not in the database (or none given): the start of its sequence
        label += ' ' + (bare[:10] + '\u2026' if len(bare) > 10 else bare)
    if shifts:
        label += ' (' + ', '.join(_shift_text(s) for s, _, _ in shifts) + ')'
    return {'label': label, 'entry': entry, 'start': start, 'end': end,
            'shifts': '; '.join(details) or None}


def describe_proteoforms(forms, sequences):
    """`forms` (accession, proforma, mass) with the columns of
    describe_proteoform: label, entry, start, end and shifts."""
    rows = [
        describe_proteoform(a, p, m, sequences)
        for a, p, m in zip(forms['accession'], forms['proforma'], forms['mass'])
    ]
    described = pl.DataFrame(rows, schema={
        'label': pl.Utf8, 'entry': pl.Utf8, 'start': pl.Int64, 'end': pl.Int64, 'shifts': pl.Utf8,
    })
    return forms.with_columns(described.get_columns())


def _rollup(long, mapping):
    """Sum `long` (id, sample, intensity, status) into the parent ids of
    `mapping` (id, parent_id). Status: the best status of the parts that have
    a value (missing when none has)."""
    rank = {s: i for i, s in enumerate(STATUS_ORDER)}
    worst = len(STATUS_ORDER) - 1
    return (
        long.join(mapping, on='id', how='inner')
        .with_columns(
            pl.when(pl.col('intensity').is_not_null())
            .then(pl.col('status').replace_strict(rank, default=worst))
            .otherwise(worst).alias('_rank'),
        )
        .group_by('parent_id', 'sample')
        .agg(
            pl.col('intensity').sum().alias('intensity'),
            pl.col('intensity').is_not_null().any().alias('_any'),
            pl.col('_rank').min().alias('_rank'),
        )
        .select(
            pl.col('parent_id').alias('id'), 'sample',
            pl.when(pl.col('_any')).then(pl.col('intensity')).otherwise(None).alias('intensity'),
            pl.col('_rank').replace_strict(dict(enumerate(STATUS_ORDER)), return_dtype=pl.Utf8).alias('status'),
            pl.lit(None, pl.Float64).alias('score'),
            pl.lit(None, pl.Float64).alias('qvalue'),
        )
        .sort('id', 'sample')
    )


def build(consensus_tsv, id_tsv=None, database=None):
    """The levels of one comparison: {level: (long, annotations)} and the run
    names. Without identifications (or with an empty -out_id) only the
    feature level exists. With the `database` (FASTA) searched, a
    proteoform's label names the residues it covers."""
    features, long, runs = read_consensus(consensus_tsv)
    feature_annotations = features.with_columns(
        pl.col('ConsensusIndex').cast(pl.Utf8).alias('id'),
        (pl.col('MonoisotopicMass').round(2).cast(pl.Utf8) + ' Da @ '
         + (pl.col('ApexRetentionTime')).round(2).cast(pl.Utf8) + ' min').alias('label'),
    ).rename({'MonoisotopicMass': 'mass', 'ApexRetentionTime': 'rt'})
    long = long.with_columns(pl.col('ConsensusIndex').cast(pl.Utf8).alias('id')).drop('ConsensusIndex')

    levels = {}
    ids = read_identifications(id_tsv) if id_tsv else None
    parents = pl.DataFrame(schema={'id': pl.Utf8, 'parent': pl.Utf8, 'ambiguous': pl.Boolean})
    if ids is not None and not ids.is_empty():
        per_feature = ids.group_by('ConsensusIndex').agg(
            pl.col('proteoform').unique().alias('_forms'),
            pl.len().alias('PrSMs'),
            pl.col('Run').n_unique().alias('IdentifiedRuns'),
        )
        parents = per_feature.select(
            pl.col('ConsensusIndex').cast(pl.Utf8).alias('id'),
            pl.when(pl.col('_forms').list.len() == 1).then(pl.col('_forms').list.first()).otherwise(None).alias('parent'),
            (pl.col('_forms').list.len() > 1).alias('ambiguous'),
            pl.col('PrSMs'), pl.col('IdentifiedRuns'),
        )
    feature_annotations = feature_annotations.join(parents, on='id', how='left').with_columns(
        pl.col('ambiguous').fill_null(False),
    )
    long = long.join(parents.select('id', 'parent'), on='id', how='left')
    levels['feature'] = (long.select(list(LONG_SCHEMA)).cast(LONG_SCHEMA), feature_annotations)

    if ids is not None and not ids.is_empty():
        to_proteoform = parents.filter(pl.col('parent').is_not_null()).select('id', pl.col('parent').alias('parent_id'))
        forms = ids.group_by('proteoform').agg(
            pl.col('ProteinAccession').first().alias('accession'),
            pl.col('ProteinDescription').drop_nulls().first().alias('description'),
            pl.col('Proforma').drop_nulls().first().alias('proforma'),
            pl.col('ProteoformMass').median().alias('mass'),
            pl.col('ProteoformLevelQvalue').min().alias('qvalue'),
            pl.col('ConsensusIndex').n_unique().alias('features'),
            pl.len().alias('PrSMs'),
            pl.col('Run').n_unique().alias('IdentifiedRuns'),
        ).filter(pl.col('proteoform').is_in(to_proteoform['parent_id'].implode()))
        form_annotations = describe_proteoforms(
            forms.with_columns(pl.col('proteoform').alias('id')), read_sequences(database),
        )
        form_long = _rollup(long.select('id', 'sample', 'intensity', 'status'), to_proteoform)
        form_long = form_long.join(
            form_annotations.select('id', pl.col('accession').alias('parent')), on='id', how='left',
        )
        levels['proteoform'] = (form_long.select(list(LONG_SCHEMA)).cast(LONG_SCHEMA), form_annotations)

        to_protein = form_annotations.select('id', pl.col('accession').alias('parent_id'))
        protein_long = _rollup(form_long.select('id', 'sample', 'intensity', 'status'), to_protein).with_columns(
            pl.lit(None, pl.Utf8).alias('parent'),
        )
        protein_annotations = form_annotations.group_by('accession').agg(
            pl.col('description').drop_nulls().first(),
            pl.len().alias('proteoforms'),
            pl.col('PrSMs').sum(),
        ).with_columns(
            pl.col('accession').alias('id'),
            pl.col('accession').str.split('|').list.get(-1, null_on_oob=True).alias('label'),
        )
        levels['protein'] = (protein_long.select(list(LONG_SCHEMA)).cast(LONG_SCHEMA), protein_annotations)
    return levels, runs


def manifest(levels, sources):
    """quant.json: contract version, levels with their noun and roll-ups."""
    nouns = {'feature': 'feature', 'proteoform': 'proteoform', 'protein': 'protein'}
    return {
        'contract': CONTRACT_VERSION,
        'levels': [level for level in LEVELS if level in levels],
        'entity': {level: nouns[level] for level in levels},
        'rollup': {
            level: {'from': below, 'method': 'sum'}
            for below, level in zip(LEVELS, LEVELS[1:]) if level in levels
        },
        'scale': 'linear',
        'label_column': 'label',
        'sources': sources,
    }


def wide(long, annotations, samples):
    """Entity x sample matrix for the comparison chain: ConsensusIndex (an
    integer per entity, in label order), Proteoform (the label, unique),
    MonoisotopicMass, RT and one intensity column per sample."""
    matrix = long.filter(pl.col('sample').is_in(samples)).pivot(
        on='sample', index='id', values='intensity', aggregate_function='first',
    )
    for sample in samples:
        if sample not in matrix.columns:
            matrix = matrix.with_columns(pl.lit(None, pl.Float64).alias(sample))
    info = annotations.select(
        'id', 'label',
        (pl.col('mass') if 'mass' in annotations.columns else pl.lit(None, pl.Float64)).alias('MonoisotopicMass'),
        (pl.col('rt') if 'rt' in annotations.columns else pl.lit(None, pl.Float64)).alias('RT'),
    )
    matrix = info.join(matrix, on='id', how='inner').sort('label', 'id')
    # labels must be unique: the comparison chain indexes by them
    matrix = matrix.with_columns(
        pl.when(pl.col('label').is_duplicated())
        .then(pl.col('label') + ' (' + pl.col('id') + ')')
        .otherwise(pl.col('label')).alias('Proteoform'),
    )
    return matrix.with_row_index('ConsensusIndex').with_columns(
        pl.col('ConsensusIndex').cast(pl.Int64),
    ).select('ConsensusIndex', 'id', 'Proteoform', 'MonoisotopicMass', 'RT', *samples)


def write(directory, levels, design, sources):
    """Write the contract into `directory` (quant.json, design.tsv,
    <level>.parquet, <level>.features.parquet)."""
    directory.mkdir(parents=True, exist_ok=True)
    for level, (long, annotations) in levels.items():
        long.write_parquet(directory / f'{level}.parquet')
        annotations.write_parquet(directory / f'{level}.features.parquet')
    design.write_csv(directory / 'design.tsv', separator='\t')
    (directory / 'quant.json').write_text(json.dumps(manifest(levels, sources), indent=2))
