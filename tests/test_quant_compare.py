"""Linking FLASHQuant runs and testing conditions (src/parse/quant_compare.py)."""

import numpy as np
import polars as pl
import pytest

pytest.importorskip('pyopenms')

from src.parse import quant_compare  # noqa: E402


def _runs(n=200, up=20, seed=0):
    """Six runs of the same proteoforms (conditions A and B, three runs each);
    the first `up` proteoforms are four times as abundant in B. Masses jitter
    by 2 ppm, retention times shift per run."""
    rng = np.random.default_rng(seed)
    mass = rng.uniform(3000, 30000, n)
    rt = rng.uniform(600, 6000, n)
    quantity = rng.lognormal(18, 1, n)
    runs = {}
    for condition in 'AB':
        for replicate in range(3):
            q = quantity * rng.lognormal(0, 0.1, n)
            if condition == 'B':
                q[:up] *= 4
            runs[f'{condition}{replicate}'] = pl.DataFrame({
                'FeatureGroupIndex': np.arange(n),
                'MonoisotopicMass': mass * (1 + rng.normal(0, 2e-6, n)),
                'HighestApexRetentionTime': rt + rng.normal(0, 3, n) + 10 * replicate,
                'FeatureGroupQuantity': q,
            })
    return runs


def test_link_runs_links_each_proteoform_across_runs():
    runs = _runs()
    consensus, members = quant_compare.link_runs(runs)
    assert consensus.height == 200
    assert (consensus['Runs'] == 6).all()
    assert members.height == 6 * 200
    # Every consensus feature holds the same feature group of every run
    per_consensus = members.group_by('ConsensusIndex').agg(pl.col('FeatureGroupIndex').n_unique())
    assert (per_consensus['FeatureGroupIndex'] == 1).all()
    assert set(members.columns) == {
        'ConsensusIndex', 'Run', 'FeatureGroupIndex', 'MonoisotopicMass', 'RT', 'Quantity',
    }


def test_link_runs_keeps_runs_apart_beyond_the_mass_tolerance():
    runs = _runs(n=50)
    runs['A0'] = runs['A0'].with_columns(pl.col('MonoisotopicMass') * (1 + 50e-6))
    consensus, _ = quant_compare.link_runs(runs, mass_tolerance_ppm=10)
    assert consensus.filter(pl.col('A0').is_not_null())['Runs'].max() == 1


def test_link_runs_needs_two_runs():
    with pytest.raises(ValueError):
        quant_compare.link_runs({'A0': _runs()['A0']})


def test_differential_finds_the_changed_proteoforms():
    pytest.importorskip('openms_insight')
    runs = _runs()
    consensus, _ = quant_compare.link_runs(runs)
    design = {name: name[0] for name in runs}
    statistics, normalized = quant_compare.differential(consensus, design)
    up = statistics.filter((pl.col('p-adj') < 0.05) & (pl.col('log2FC') > 1))
    # Exactly the proteoforms four times as abundant in B (log2FC near 2)
    assert up.height == 20
    assert up['log2FC'].is_between(1.5, 2.5).all()
    assert normalized.columns == ['ConsensusIndex', 'Run', 'Abundance', 'Condition']
    assert normalized.height == statistics.height * 6

    heat = quant_compare.zscores(statistics, list(design), top=10)
    assert heat.height == 10
    assert np.allclose(heat.select(list(design)).to_numpy().mean(axis=1), 0)


def test_differential_needs_two_conditions():
    pytest.importorskip('openms_insight')
    runs = _runs(n=20)
    consensus, _ = quant_compare.link_runs(runs)
    with pytest.raises(ValueError):
        quant_compare.differential(consensus, {name: 'A' for name in runs})
