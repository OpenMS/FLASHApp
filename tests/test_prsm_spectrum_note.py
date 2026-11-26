"""The note shown when no PrSM of the picked identification has its spectrum
in the workspace names the runs whose results are missing."""

import polars as pl

from src.render import quant_cells


FRAME = pl.DataFrame({'RowIndex': [0, 0, 1], 'Run': ['R2', 'R3', 'R1']})


def test_names_the_runs_whose_results_are_missing():
    note = quant_cells._no_spectrum_note(FRAME, 0, {'R1': 'example'})
    assert 'R2, R3' in note
    assert 'not in this workspace' in note


def test_general_note_when_every_run_is_here():
    note = quant_cells._no_spectrum_note(FRAME, 1, {'R1': 'example'})
    assert 'R2' not in note
