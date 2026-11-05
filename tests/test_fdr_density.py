"""FDR plot densities of FLASHTnT q-values (``src/parse/tnt.py``).

A kernel density needs at least two distinct values. A run with every
identification at q = 0, or with a single one above, is a good run: its FDR
plot gets no curve on that side instead of the workflow failing.
"""
import os
import sys

import pandas as pd
import pytest

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(PROJECT_ROOT)

pytest.importorskip("scipy")
pytest.importorskip("pyopenms")

from src.parse.tnt import fdr_density_distribution  # noqa: E402


class _Logger:
    def __init__(self):
        self.lines = []

    def log(self, message, level=0):
        self.lines.append(message)


def _proteins(targets, decoys=()):
    return pd.DataFrame({
        'accession': [f'P{i}' for i in range(len(targets))] + [f'DECOY_P{i}' for i in range(len(decoys))],
        'ProteoformLevelQvalue': list(targets) + list(decoys),
    })


def test_spread_values_give_both_densities():
    target, decoy = fdr_density_distribution(_proteins([0.01, 0.02, 0.05], [0.3, 0.5, 0.9]))
    assert len(target) == 200 and len(decoy) == 200
    assert target['y'].gt(0).all() and decoy['y'].gt(0).all()


def test_single_target_above_zero_does_not_fail():
    logger = _Logger()
    target, decoy = fdr_density_distribution(_proteins([0.0, 0.0, 0.004]), logger=logger)
    assert target.empty and list(target.columns) == ['x', 'y']
    assert decoy.empty
    assert any('target' in line for line in logger.lines)


def test_identical_values_do_not_fail():
    target, decoy = fdr_density_distribution(_proteins([0.01, 0.01], [0.2, 0.2, 0.2]))
    assert target.empty and decoy.empty


def test_all_at_zero_gives_empty_densities():
    logger = _Logger()
    target, decoy = fdr_density_distribution(_proteins([0.0] * 26), logger=logger)
    assert target.empty and decoy.empty
    assert logger.lines == []


def test_one_side_empty_other_drawn():
    target, decoy = fdr_density_distribution(_proteins([0.001], [0.1, 0.4, 0.7]))
    assert target.empty
    assert len(decoy) == 200
