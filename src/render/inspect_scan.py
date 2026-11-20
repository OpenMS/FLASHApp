"""Inspect a spectrum of a run on the Match fragments page.

The inspect button of a row of the Viewer's spectrum tables (``TABLES``)
takes that row's spectrum: its
deconvolved masses and precursor mass, and, when the scan's proteoform is
the one picked, the proteoform's residues, and hands them to the matcher
(``fragment_matcher.store_input``), where masses, sequence, modifications
and ion types can be changed. Streamlit-free: the page switches.
"""

from __future__ import annotations

import polars as pl

from src import fragment_matcher
from src.render import insight_data


# The tables whose rows are a spectrum, and how to read a row: (the
# spectrum's index, the proteoform of the row or None)
TABLES = {
    'scan_table': lambda row: (row.get('index'), None),
    'protein_table': lambda row: (row.get('scan_idx'), row.get('index')),
    'proteoform_table': lambda row: (row.get('scan_idx'), row.get('best_index')),
    'proteoform_prsm_table': lambda row: (row.get('scan_idx'), row.get('index')),
}


def from_row(table, row):
    """(scan_idx, proteoform_index) of a row of one of ``TABLES``."""
    return TABLES[table](row)


def proteoform_residues(file_manager, dataset_id, proteoform_index, scan_idx=None):
    """The residues of an identified proteoform (its stretch of the protein),
    or None when the run holds no such proteoform, or when ``scan_idx`` is
    given and the proteoform was found in another spectrum."""
    if proteoform_index is None or not file_manager.result_exists(dataset_id, 'sequence_data'):
        return None
    rows = (
        insight_data.proteoform_sequences(file_manager, dataset_id)
        .filter(pl.col('proteoform_index') == int(proteoform_index))
        .select('sequence', 'start_position', 'end_position', 'scan_idx')
        .collect()
    )
    if rows.height == 0:
        return None
    row = rows.row(0, named=True)
    if scan_idx is not None and row['scan_idx'] is not None and int(row['scan_idx']) != int(scan_idx):
        return None
    residues = list(row['sequence'] if row['sequence'] is not None else [])
    start = max(int(row['start_position'] or 0), 0)
    end = row['end_position']
    end = len(residues) - 1 if end is None or end < 0 else min(int(end), len(residues) - 1)
    return ''.join(residues[start:end + 1]) or None


def scan_input(file_manager, dataset_id, scan_idx, proteoform_index=None, run_name=None):
    """What the matcher shows for one spectrum of a run: {'masses',
    'precursor_mass', 'sequence', 'name'}, or None when the spectrum has
    no deconvolved masses."""
    if scan_idx is None:
        return None
    spectra = insight_data.sequence_view_spectra(file_manager, dataset_id)
    rows = spectra.filter(pl.col('index') == int(scan_idx)).collect()
    values = rows['MonoMass'][0] if rows.height else None
    masses = [] if values is None else [float(m) for m in values]
    if not masses:
        return None
    precursor = rows['PrecursorMass'][0]
    scan_number = _scan_number(file_manager, dataset_id, scan_idx)
    label = f'scan {scan_number}' if scan_number is not None else f'spectrum {int(scan_idx) + 1}'
    return {
        'masses': masses,
        'precursor_mass': float(precursor) if precursor else None,
        'sequence': proteoform_residues(file_manager, dataset_id, proteoform_index, scan_idx),
        'name': f'{run_name}, {label}' if run_name else label.capitalize(),
    }


def store(workspace, data):
    """Save ``scan_input``'s result as the matcher's input."""
    return fragment_matcher.store_input(
        workspace, data['masses'], None, data['precursor_mass'],
        fragment_matcher.MONOISOTOPIC, sequence=data['sequence'], name=data['name'],
        hand_over=True,
    )


def _scan_number(file_manager, dataset_id, scan_idx):
    if not file_manager.result_exists(dataset_id, 'scan_table'):
        return None
    rows = (
        insight_data.scan(file_manager, dataset_id, 'scan_table')
        .filter(pl.col('index') == int(scan_idx))
        .select('Scan')
        .collect()
    )
    return int(rows['Scan'][0]) if rows.height and rows['Scan'][0] is not None else None
