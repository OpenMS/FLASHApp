"""Explore results opens on what the last start made (src/analysis.py)."""

from src.analysis import newest_result


def test_newest_run_wins():
    ids = ['example_deconvolution', 'A_20261106-101500', 'B_20261106-191000', 'C_20261105-235959']
    assert newest_result(ids) == 'B_20261106-191000'


def test_comparison_made_with_its_runs_wins_a_tie():
    ids = ['A_20261106-191000', 'comparison_20261106-191000', 'B_20261106-190959']
    assert newest_result(ids) == 'comparison_20261106-191000'


def test_examples_only():
    assert newest_result(['example_native_identified', 'example_label_free']) is None
    assert newest_result([]) is None
