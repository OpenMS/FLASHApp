"""A comparison of runs: the quantification of one FLASHQuantConsensus run
(src/quant/contract.py) in the workspace cache, its experimental design and
analysis settings, and the comparison chain's results per level.

Everything lives under one dataset of the FLASHDeconv/FLASHTnT cache, so the
Viewer lists a comparison next to the runs it links. Stored names:

- ``quant_manifest``: quant.json (levels, roll-ups, sources) plus the runs
- ``quant_<level>`` / ``quant_<level>_annotations``: the contract tables
- ``quant_design``: sample, condition, replicate
- ``quant_settings``: the analysis settings and the level shown
- ``quant_<level>_<result>``: the comparison chain's results for that level,
  recomputed when the design or the settings change (``quant_<level>_stamp``)

Streamlit-free.
"""

import hashlib
import json

import polars as pl

from src.parse import quant_compare
from src.quant import contract

RESULTS = ('statistics', 'normalized', 'filtered', 'imputed', 'normalized_wide', 'zscores')
LEVEL_NAMES = {'feature': 'Features', 'proteoform': 'Proteoforms', 'protein': 'Proteins'}
STATUS_LABELS = {
    'detected': 'detected',
    'transferred_identified': 'transferred, identified',
    'transferred': 'transferred',
    'missing': 'missing',
}


def is_comparison(file_manager, dataset_id):
    """True for a dataset that holds a comparison of runs."""
    return file_manager.result_exists(dataset_id, 'quant_manifest')


def comparisons(file_manager, tool=None):
    """The comparison datasets of a cache, newest first; with `tool`, only
    those a run of that workflow ('flashdeconv', 'flashtnt') stored."""
    found = sorted(file_manager.get_results_list(['quant_manifest']), reverse=True)
    if tool is None:
        return found
    return [
        dataset_id for dataset_id in found
        if file_manager.get_results(dataset_id, ['quant_manifest'])['quant_manifest'].get('tool') == tool
    ]


def store(file_manager, dataset_id, levels, runs, design, sources, tool, extra=None):
    """Store a comparison of `tool` ('flashdeconv' or 'flashtnt'): the
    contract tables, the design (sample -> (condition, replicate)) and
    default settings. `runs` are the compared samples: the runs of a
    label-free comparison, the labelled samples of an isobaric one (`extra`
    then says kind 'isobaric', the plex, its channels and the fractions)."""
    for level, (long, annotations) in levels.items():
        file_manager.store_data(dataset_id, f'quant_{level}', long)
        file_manager.store_data(dataset_id, f'quant_{level}_annotations', annotations)
    file_manager.store_data(dataset_id, 'quant_design', design_frame(runs, design))
    file_manager.store_data(dataset_id, 'quant_settings', {
        'level': 'proteoform' if 'proteoform' in levels else 'feature',
    })
    file_manager.store_data(
        dataset_id, 'quant_manifest',
        {**contract.manifest(levels, sources), 'runs': runs, 'tool': tool, **(extra or {})},
    )


def design_frame(runs, design):
    """sample, condition, replicate for every run; `design` maps a run to
    (condition, replicate)."""
    return pl.DataFrame({
        'sample': runs,
        'condition': [str(design.get(run, ('', 0))[0] or '') for run in runs],
        'replicate': [int(design.get(run, ('', i + 1))[1] or i + 1) for i, run in enumerate(runs)],
    })


class QuantSession:
    """One comparison dataset: its tables, design and settings, and the
    comparison chain run on the level shown."""

    def __init__(self, file_manager, dataset_id):
        self.file_manager = file_manager
        self.dataset_id = dataset_id
        self.manifest = self._get('quant_manifest')

    def _get(self, name, polars=False):
        return self.file_manager.get_results(self.dataset_id, [name], use_polars=polars)[name]

    def exists(self, name):
        return self.file_manager.result_exists(self.dataset_id, name)

    # -- what is stored -------------------------------------------------
    @property
    def isobaric(self):
        """True for a comparison of labelled samples (TMT, iTRAQ)."""
        return self.manifest.get('kind') == 'isobaric'

    @property
    def sample_noun(self):
        return 'sample' if self.isobaric else 'run'

    @property
    def levels(self):
        return list(self.manifest['levels'])

    @property
    def runs(self):
        return list(self.manifest['runs'])

    @property
    def settings(self):
        from src.parse.differential import default_settings  # pylint: disable=import-outside-toplevel
        stored = self._get('quant_settings') if self.exists('quant_settings') else {}
        return {**default_settings(), 'level': self.levels[0], **stored}

    @property
    def level(self):
        level = self.settings.get('level')
        return level if level in self.levels else self.levels[0]

    def design(self):
        """sample, condition, replicate (a polars DataFrame)."""
        return self._get('quant_design', polars=True).collect()

    def conditions(self):
        """sample -> condition, for the runs with a condition."""
        design = self.design()
        return {s: c for s, c in zip(design['sample'], design['condition']) if c}

    @property
    def compared(self):
        """The runs (samples) in the comparison: those with a condition, in
        run order; every run while none has one."""
        conditions = self.conditions()
        return [r for r in self.runs if r in conditions] or self.runs

    def long(self, level=None):
        return self._get(f'quant_{level or self.level}', polars=True)

    def annotations(self, level=None):
        return self._get(f'quant_{level or self.level}_annotations', polars=True)

    def save_design(self, design):
        self.file_manager.store_data(self.dataset_id, 'quant_design', design)

    def save_settings(self, settings):
        self.file_manager.store_data(self.dataset_id, 'quant_settings', settings)

    # -- the comparison chain --------------------------------------------
    def problem(self):
        """Why the level cannot be compared (None when it can)."""
        conditions = self.conditions()
        groups = set(conditions.values())
        if len(groups) < 2:
            return ('Give the runs at least two conditions (**Edit comparison**) to compare them. '
                    'Runs of one condition are its replicates.')
        if len(conditions) < 3:
            return 'Comparing conditions needs at least three runs with a condition.'
        return None

    def _stamp(self, level):
        settings = {k: v for k, v in self.settings.items() if k != 'level'}
        parts = [self.conditions(), settings, level, contract.CONTRACT_VERSION]
        return hashlib.md5(json.dumps(parts, sort_keys=True, default=str).encode()).hexdigest()[:12]

    def result_name(self, name, level=None):
        return f'quant_{level or self.level}_{name}'

    def analyse(self, level=None):
        """Run the comparison chain on `level` unless its results are current.
        Returns the error message when it cannot run (None when the results
        are there)."""
        level = level or self.level
        problem = self.problem()
        if problem:
            return problem
        stamp = self._stamp(level)
        name = self.result_name('stamp', level)
        if self.exists(name) and self._get(name) == stamp and all(
            self.exists(self.result_name(r, level)) for r in RESULTS
        ):
            return None
        conditions = self.conditions()
        samples = [s for s in self.runs if s in conditions]
        matrix = contract.wide(self.long(level).collect(), self.annotations(level).collect(), samples)
        if matrix.height == 0:
            return 'Nothing to compare at this level.'
        try:
            statistics, normalized, steps = quant_compare.differential(matrix, conditions, self.settings)
        except ValueError as e:
            return str(e)
        ids = matrix.select('ConsensusIndex', 'id')
        statistics = statistics.join(ids, on='ConsensusIndex', how='left')
        status = self.long(level).collect().select(
            'id', pl.col('sample').alias('Run'),
            pl.col('status').replace_strict(STATUS_LABELS, default='missing').alias('Status'),
        )
        normalized = (
            normalized.join(ids, on='ConsensusIndex', how='left')
            .join(status, on=['id', 'Run'], how='left')
            .with_columns(pl.col('Status').fill_null('missing'))
        )
        results = {
            'statistics': statistics, 'normalized': normalized,
            'filtered': steps['filtered'], 'imputed': steps['imputed'],
            'normalized_wide': steps['normalized'],
            'zscores': quant_compare.zscores(statistics, samples),
        }
        for result, frame in results.items():
            self.file_manager.store_data(self.dataset_id, self.result_name(result, level), frame)
        self.file_manager.store_data(self.dataset_id, name, stamp)
        return None

    def result(self, name, level=None):
        return self._get(self.result_name(name, level), polars=True)

    def table(self, level=None):
        """The comparison table of a level: one row per entity with its label,
        annotations, test results and its status in every run (status
        columns named after the runs). `Detected` counts the compared runs
        (`compared`) it was detected in."""
        level = level or self.level
        annotations = self.annotations(level).collect()
        statistics = self.result('statistics', level).collect().select(
            'ConsensusIndex', 'id', 'Proteoform', 'log2FC', 'p-value', 'p-adj',
        )
        status = (
            self.long(level).collect()
            .select('id', 'sample', pl.col('status').replace_strict(STATUS_LABELS, default='missing'))
            .pivot(on='sample', index='id', values='status', aggregate_function='first')
        )
        keep = [c for c in ('accession', 'description', 'proforma', 'start', 'end', 'shifts', 'mass', 'rt', 'PrSMs',
                            'IdentifiedRuns', 'proteoforms', 'DetectedRuns', 'parent') if c in annotations.columns]
        table = (
            statistics.join(annotations.select('id', *keep), on='id', how='left')
            .join(status, on='id', how='left')
            .with_columns(
                pl.col('p-adj').is_not_null().alias('Tested'),
                (pl.concat_list([pl.col(r) == 'detected' for r in self.compared]).list.sum()).alias('Detected'),
            )
        )
        return table
