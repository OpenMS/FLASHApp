"""
FLASHAnnotator output -> per-proteoform fragment annotations
(``src/parse/annotator.py``).

The viewer draws the full database sequence, while FLASHAnnotator reports
positions within the proteoform; these tests pin the coordinate conversion and
the parsing of terminal and internal envelopes.
"""
import io
import os
import sys

import pandas as pd

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(PROJECT_ROOT)

from src.parse.annotator import (
    build_fragments, group_annotation, proteoform_offset, read_annotation,
    read_summary,
)

ANNOTATION_TSV = (
    "PrSMIndex\tScan\tIonType\tLength\tBond\tCharge\tMonoMass\tMostAbundantMz\tMatchedPeaks"
    "\tPredictedPeaks\tCosine\tPearson\tIntensity\tPPMError\tIntensityRank\tLocalDensity\tDecoyShift\tInternalFirst\n"
    "0\t5\tb\t3\t3\t1\t300.10000\t301.10000\t3\t4\t0.9500\t0.9\t1000.0\t1.00\t0.9\t1.0\t0.000\t0\n"
    "0\t5\ty\t4\t6\t2\t500.20000\t251.10000\t4\t5\t0.9000\t0.8\t2000.0\t-2.00\t0.9\t1.0\t0.000\t0\n"
    "0\t5\tint(b/y)\t5\t0\t1\t600.30000\t601.30000\t5\t5\t0.8000\t0.7\t500.0\t0.50\t0.9\t1.0\t0.000\t3\n"
    "0\t5\tb\t2\t2\t1\t200.00000\t201.00000\t2\t2\t0.9000\t0.9\t10.0\t0.00\t0.9\t1.0\t1.500\t0\n"
    "1\t6\tzp1\t2\t8\t1\t250.00000\t251.00000\t2\t3\t0.9900\t0.9\t50.0\t0.00\t0.9\t1.0\t0.000\t0\n"
)


def _annotation():
    return read_annotation(io.StringIO(ANNOTATION_TSV))


def test_decoy_envelopes_are_dropped():
    df = _annotation()
    assert len(df) == 4
    assert (df['DecoyShift'] == 0).all()


def test_proteoform_offset():
    db = 'MAAAKPEPTIDEKR'
    assert proteoform_offset(db, 'AAAKPEP', 2, 8) == 1
    # open N-terminus: ends at the known end
    assert proteoform_offset(db, 'PEPTIDE', -1, 12) == 5
    # both open: where the proteoform occurs
    assert proteoform_offset(db, 'TIDE', -1, -1) == 8


def test_build_fragments_coordinates():
    groups = group_annotation(_annotation())
    # proteoform of 10 residues starting at database residue 2 (offset 2)
    fragments = build_fragments(groups[0], offset=2, proteoform_length=10)
    by_name = {f['name']: f for f in fragments}
    b3 = by_name['b3']
    assert (b3['kind'], b3['first'], b3['last'], b3['start'], b3['end']) == ('N', 0, 2, 2, 4)
    y4 = by_name['y4']
    assert (y4['kind'], y4['first'], y4['last'], y4['start'], y4['end']) == ('C', 6, 9, 8, 11)
    internal = by_name['b/y[3-7]']
    assert internal['kind'] == 'I'
    assert internal['ion_type'] == 'b/y'
    assert (internal['first'], internal['last'], internal['start'], internal['end']) == (2, 6, 4, 8)
    assert abs(y4['observed_mass'] - 500.2 * (1 - 2e-6)) < 1e-9
    # terminal fragments first, internal last
    assert [f['kind'] for f in fragments] == ['N', 'C', 'I']


def test_z_dot_ion_label():
    groups = group_annotation(_annotation())
    fragments = build_fragments(groups[1], offset=0, proteoform_length=10)
    assert fragments[0]['name'] == 'z+1(2)'
    assert fragments[0]['ion_type'] == 'z+1'
    assert fragments[0]['kind'] == 'C'


def test_missing_files_are_empty():
    assert read_annotation(None).empty
    assert read_summary(None) == {}
    assert group_annotation(read_annotation(None)) == {}


def test_read_summary_marks_uncomputed_values():
    tsv = io.StringIO(
        "ProteoformIndex\tScan\tFAScore\tRawCoverage(%)\tRawCalibrationPPM\tRawCoverageWithInternal(%)\n"
        "0\t5\t12.5\t40.0\tnan\t55.0\n"
        "1\t6\t-1\t-1\tnan\t-1\n"
    )
    summary = read_summary(tsv)
    assert summary[0]['fa_score'] == 12.5
    assert summary[0]['coverage_with_internal'] == 55.0
    assert summary[0]['calibration_ppm'] is None
    assert summary[1]['fa_score'] is None
