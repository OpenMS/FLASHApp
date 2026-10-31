"""Compare FLASHQuant results of several runs.

Two steps, both free of Streamlit so they can run anywhere:

1. link_runs: the feature groups (proteoforms) of every run are linked
   across runs by monoisotopic mass and apex retention time with OpenMS's
   KD-tree feature linker (FeatureGroupingAlgorithmKD, the algorithm of the
   TOPP tool FeatureLinkerUnlabeledKD, which also corrects retention time
   shifts between the runs). A feature group enters the linker as a feature
   whose "m/z" is its monoisotopic mass, so the m/z tolerance in ppm is a
   mass tolerance.
2. differential: abundances of the linked proteoforms are compared between
   two or more conditions with openms-insight's analysis functions (filter,
   impute, log2 transform, median normalization, moderated t-test or
   F-test, Benjamini-Hochberg).
"""

import numpy as np
import polars as pl
import pyopenms

# Defaults of the linker: FLASHQuant masses agree within a few ppm; apexes of
# one proteoform within a minute after the linker's retention time correction
MASS_TOLERANCE_PPM = 10.0
RT_TOLERANCE_S = 60.0
# A proteoform must be quantified in at least this share of the runs of one
# condition to be tested
MIN_PRESENT_RATIO = 0.5

ID_COLUMN = 'Proteoform'


def _feature_map(quant, run_index, uid_start):
    """FeatureMap of one run's feature groups; unique ids from uid_start."""
    feature_map = pyopenms.FeatureMap()
    for offset, row in enumerate(quant.iter_rows(named=True)):
        feature = pyopenms.Feature()
        feature.setMZ(float(row['MonoisotopicMass']))
        feature.setRT(float(row['HighestApexRetentionTime']))
        feature.setIntensity(float(row['FeatureGroupQuantity']))
        feature.setCharge(0)
        feature.setUniqueId(uid_start + offset)
        feature_map.push_back(feature)
    feature_map.setUniqueId(run_index + 1)
    return feature_map


def link_runs(runs, mass_tolerance_ppm=MASS_TOLERANCE_PPM, rt_tolerance=RT_TOLERANCE_S):
    """Link the feature groups of several runs.

    `runs` maps a run name to its FLASHQuant feature groups (columns
    FeatureGroupIndex, MonoisotopicMass, HighestApexRetentionTime,
    FeatureGroupQuantity). Returns (consensus, members):

    - consensus: one row per linked proteoform: ConsensusIndex, Proteoform
      (a label from its mass and retention time), MonoisotopicMass and RT
      (averages of the linked groups), Runs (in how many runs it was found)
      and one quantity column per run (null where it was not found).
    - members: every linked feature group: ConsensusIndex, Run,
      FeatureGroupIndex, MonoisotopicMass, RT (before retention time
      correction) and Quantity.
    """
    names = list(runs)
    if len(names) < 2:
        raise ValueError('Comparing runs needs at least two runs.')
    maps, group_indices, uid = [], {}, 1
    for run_index, name in enumerate(names):
        quant = runs[name].select(
            'FeatureGroupIndex', 'MonoisotopicMass', 'HighestApexRetentionTime',
            'FeatureGroupQuantity',
        )
        maps.append(_feature_map(quant, run_index, uid))
        for offset, row in enumerate(quant.iter_rows()):
            group_indices[uid + offset] = (name, int(row[0]), float(row[1]), float(row[2]))
        uid += quant.height

    linker = pyopenms.FeatureGroupingAlgorithmKD()
    params = linker.getDefaults()
    params.setValue('mz_unit', 'ppm')
    params.setValue('link:mz_tol', float(mass_tolerance_ppm))
    params.setValue('link:rt_tol', float(rt_tolerance))
    params.setValue('warp:mz_tol', float(mass_tolerance_ppm))
    linker.setParameters(params)
    consensus_map = pyopenms.ConsensusMap()
    linker.group(maps, consensus_map)

    consensus_rows, member_rows = [], []
    for consensus_index, feature in enumerate(consensus_map):
        quantities = dict.fromkeys(names)
        for handle in feature.getFeatureList():
            name, group_index, group_mass, group_rt = group_indices[handle.getUniqueId()]
            quantities[name] = float(handle.getIntensity())
            member_rows.append((
                consensus_index, name, group_index, group_mass, group_rt,
                float(handle.getIntensity()),
            ))
        mass, rt = feature.getMZ(), feature.getRT()
        consensus_rows.append({
            'ConsensusIndex': consensus_index,
            ID_COLUMN: f'{mass:.2f} Da @ {rt / 60:.1f} min',
            'MonoisotopicMass': mass,
            'RT': rt,
            'Runs': sum(v is not None for v in quantities.values()),
            **quantities,
        })
    schema = {
        'ConsensusIndex': pl.Int64, ID_COLUMN: pl.String, 'MonoisotopicMass': pl.Float64,
        'RT': pl.Float64, 'Runs': pl.Int64, **{name: pl.Float64 for name in names},
    }
    consensus = pl.DataFrame(consensus_rows, schema=schema)
    members = pl.DataFrame(
        member_rows, orient='row', schema={
            'ConsensusIndex': pl.Int64, 'Run': pl.String, 'FeatureGroupIndex': pl.Int64,
            'MonoisotopicMass': pl.Float64, 'RT': pl.Float64, 'Quantity': pl.Float64,
        },
    )
    return consensus, members


def _metadata(design):
    return pl.DataFrame(
        [{'sample_id': run, 'group': condition} for run, condition in design.items()],
        schema={'sample_id': pl.String, 'group': pl.String},
    )


def differential(consensus, design, min_present_ratio=MIN_PRESENT_RATIO):
    """Differential abundance of the linked proteoforms between conditions.

    `design` maps run names (quantity columns of `consensus`) to conditions;
    runs without a condition are left out. Returns (statistics, normalized):

    - statistics: one row per tested proteoform with ConsensusIndex,
      Proteoform, MonoisotopicMass, RT, the normalized log2 abundance per
      run, log2FC (second condition over the first, in sorted order, for two
      conditions), stat, p-value and p-adj.
    - normalized: the normalized log2 abundances in long form
      (ConsensusIndex, Run, Condition, Abundance) for abundance profiles.
    """
    # pylint: disable=import-outside-toplevel
    # Imported here: linking runs does not need openms-insight
    from openms_insight.analysis.filter import filter_low_repeatability
    from openms_insight.analysis.imputation import impute_smallest_value
    from openms_insight.analysis.normalization import normalize_samples, transform_data
    from openms_insight.analysis.statistics import adjust_fdr_lazy, calculate_statistical_tests

    design = {run: condition for run, condition in design.items() if condition}
    conditions = sorted(set(design.values()))
    if len(conditions) < 2:
        raise ValueError('Comparing conditions needs runs of at least two conditions.')
    metadata = _metadata(design)
    runs = list(design)
    data = consensus.select('ConsensusIndex', ID_COLUMN, 'MonoisotopicMass', 'RT', *runs).lazy()
    data = filter_low_repeatability(
        quantification_data=data, metadata=metadata, group_column='group',
        max_missing_ratio=1.0 - min_present_ratio,
    )
    data = impute_smallest_value(quantification_data=data, metadata=metadata, scope='row')
    data = transform_data(quantification_data=data, metadata=metadata, strategy='log2')
    data = normalize_samples(
        quantification_data=data, metadata=metadata, strategy='median', id_col=ID_COLUMN,
        reference_feature=None,
    )
    method = 'limma_like'
    data = calculate_statistical_tests(quantification_data=data, metadata=metadata, method=method)
    statistics = adjust_fdr_lazy(quantification_data=data, strategy='BH').collect()

    normalized = (
        statistics.select('ConsensusIndex', *runs)
        .unpivot(index='ConsensusIndex', variable_name='Run', value_name='Abundance')
        .join(metadata.rename({'sample_id': 'Run', 'group': 'Condition'}), on='Run')
        .with_columns(pl.col('Abundance').cast(pl.Float64))
        .sort(['ConsensusIndex', 'Condition', 'Run'])
    )
    return statistics, normalized


def zscores(statistics, runs, top=50):
    """Row z-scores of the `top` proteoforms with the smallest adjusted
    p-values, for the clustered heatmap (wide: Proteoform plus one column per
    run)."""
    ranked = statistics.sort('p-adj', nulls_last=True).head(top)
    values = ranked.select(runs).to_numpy().astype(float)
    mean = np.nanmean(values, axis=1, keepdims=True)
    std = np.nanstd(values, axis=1, keepdims=True)
    std[std == 0] = 1.0
    scaled = (values - mean) / std
    return pl.DataFrame(
        {ID_COLUMN: ranked[ID_COLUMN], **{run: scaled[:, i] for i, run in enumerate(runs)}}
    )
