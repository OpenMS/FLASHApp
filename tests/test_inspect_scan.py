"""The Viewer's Inspect action: a picked spectrum as the input of Match fragments.

Runs on the example results of the demo workspace in example-data/.
"""

import os
import shutil
import sys
from pathlib import Path

import polars as pl
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src import fragment_matcher  # noqa: E402
from src.render import inspect_scan, insight_data  # noqa: E402
from src.workflow.FileManager import FileManager  # noqa: E402

WORKSPACE = Path(__file__).resolve().parent.parent / 'example-data' / 'workspaces' / 'default'
DECONVOLVED = 'example_deconvolution'
IDENTIFIED = 'example_lysate_identified'

pytestmark = pytest.mark.skipif(
    not (WORKSPACE / 'cache' / 'cache.db').exists(), reason='demo workspace not present'
)


@pytest.fixture
def fm(tmp_path):
    # Opening the cache database writes to it: work on a copy of it
    cache = tmp_path / 'cache'
    cache.mkdir()
    shutil.copy(WORKSPACE / 'cache' / 'cache.db', cache / 'cache.db')
    (cache / 'files').symlink_to(WORKSPACE / 'cache' / 'files')
    return FileManager(cache, cache)


def _identified(fm):
    """The proteoform index and spectrum of an identification, with the
    protein sequence and the proteoform's start and end."""
    row = (
        insight_data.proteoform_sequences(fm, IDENTIFIED)
        .select('proteoform_index', 'scan_idx', 'sequence', 'start_position', 'end_position')
        .collect()
        .row(0, named=True)
    )
    return row


def test_scan_input_holds_the_deconvolved_masses_of_the_spectrum(fm):
    spectrum = (
        insight_data.sequence_view_spectra(fm, DECONVOLVED)
        .filter(pl.col('index') == 5).collect()
    )
    data = inspect_scan.scan_input(fm, DECONVOLVED, 5, run_name='Run')
    assert data['masses'] == [float(m) for m in spectrum['MonoMass'][0]]
    assert data['sequence'] is None
    scan_number = (
        insight_data.scan(fm, DECONVOLVED, 'scan_table')
        .filter(pl.col('index') == 5).collect()['Scan'][0]
    )
    assert data['name'] == f'Run, scan {scan_number}'


def test_scan_input_takes_the_residues_of_the_picked_proteoform(fm):
    row = _identified(fm)
    data = inspect_scan.scan_input(fm, IDENTIFIED, row['scan_idx'], row['proteoform_index'])
    residues = ''.join(row['sequence'])
    assert data['sequence'] == residues[row['start_position']:row['end_position'] + 1]


def test_a_proteoform_of_another_spectrum_gives_no_sequence(fm):
    row = _identified(fm)
    other = (
        insight_data.sequence_view_spectra(fm, IDENTIFIED)
        .filter(pl.col('index') != row['scan_idx']).select('index').collect()['index'][0]
    )
    data = inspect_scan.scan_input(fm, IDENTIFIED, other, row['proteoform_index'])
    assert data is not None and data['sequence'] is None


def test_no_scan_or_an_unknown_scan_gives_nothing(fm):
    assert inspect_scan.scan_input(fm, DECONVOLVED, None) is None
    assert inspect_scan.scan_input(fm, DECONVOLVED, 10 ** 7) is None


def test_store_writes_the_matcher_input(fm, tmp_path):
    row = _identified(fm)
    data = inspect_scan.scan_input(fm, IDENTIFIED, row['scan_idx'], row['proteoform_index'], 'Lysate')
    inspect_scan.store(tmp_path, data)
    stored = fragment_matcher.load_input(tmp_path)
    assert stored['masses'] == sorted(data['masses'])
    assert stored['sequence'] == data['sequence']
    assert stored['mass_type'] == fragment_matcher.MONOISOTOPIC
    assert stored['name'].startswith('Lysate, scan ')


def test_the_viewer_opens_the_matcher_page_with_its_state():
    # The Viewer writes the sequence into the state the Match fragments page reads
    root = Path(__file__).resolve().parent.parent
    viewer = (root / 'content' / 'TopDown' / 'TopDownViewer.py').read_text()
    matcher = (root / 'content' / 'Tools' / 'FragmentMatcher.py').read_text()
    assert "MATCHER_STATE = 'fragment_matcher_state'" in viewer
    assert "STATE_KEY = 'fragment_matcher_state'" in matcher
    assert (root / 'content' / 'Tools' / 'FragmentMatcher.py').exists()
    assert "MATCHER_PAGE = 'content/Tools/FragmentMatcher.py'" in viewer
