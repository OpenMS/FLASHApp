"""Tests for the frames the openms-insight viewer components are built from.

They run on the demo workspace shipped in example-data/, so they also check
that results written by earlier versions of the workflows still load.
"""

import os
import shutil
import sys
from pathlib import Path

import polars as pl
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.render import insight_data  # noqa: E402
from src.workflow.FileManager import FileManager  # noqa: E402

WORKSPACE = Path(__file__).resolve().parent.parent / 'example-data' / 'workspaces' / 'default'
FD_DATASET = 'example_fd'
TNT_DATASET = 'example_spectrum_aqpz'


def _file_manager(tmp_path, tool):
    # Opening the cache database writes to it: work on a copy of it, next to
    # the stored results
    cache = tmp_path / tool
    cache.mkdir()
    shutil.copy(WORKSPACE / tool / 'cache' / 'cache.db', cache / 'cache.db')
    (cache / 'files').symlink_to(WORKSPACE / tool / 'cache' / 'files')
    return FileManager(cache, cache)


@pytest.fixture
def fd(tmp_path):
    return _file_manager(tmp_path, 'flashdeconv')


@pytest.fixture
def tnt(tmp_path):
    return _file_manager(tmp_path, 'flashtnt')


@pytest.fixture
def quant(tmp_path):
    return _file_manager(tmp_path, 'flashquant')


def test_mass_table_has_one_row_per_mass(fd):
    masses = insight_data.mass_table(fd, FD_DATASET).collect()
    per_scan = insight_data.scan(fd, FD_DATASET, 'mass_table').collect()
    assert masses.height == per_scan['MonoMass'].list.len().sum()
    first = masses.filter(pl.col('scan_idx') == 0).sort('mass_idx')
    assert first['mass_idx'].to_list() == list(range(first.height))
    assert first['MonoMass'].to_list() == per_scan['MonoMass'][0].to_list()
    assert masses.schema['MonoMass'] == pl.Float64


def test_signal_peaks_per_mass_and_precursor(fd):
    peaks = insight_data.signal_peaks(fd, FD_DATASET).collect()
    raw = insight_data.scan(fd, FD_DATASET, 'threedim_SN_plot').collect()

    # Signal peaks of mass 0 of spectrum 0: mass = m/z * charge
    row = peaks.filter(
        (pl.col('scan_idx') == 0) & (pl.col('mass_idx') == 0) & (pl.col('kind') == 'Signal')
    )
    signal = raw['SignalPeaks'][0][0].to_list()
    assert row['mass'][0].to_list() == pytest.approx([p[1] * p[3] for p in signal])
    assert row['charge'][0].to_list() == [p[3] for p in signal]
    assert row['intensity'][0].to_list() == [p[2] for p in signal]

    # Fragmentation spectra get the peaks of their precursor under mass_idx -1
    precursor_rows = peaks.filter(pl.col('mass_idx') == -1)
    assert precursor_rows.height > 0
    ms1 = insight_data.scan_table(fd, FD_DATASET).filter(pl.col('MSLevel') == 1).collect()
    assert precursor_rows.filter(pl.col('scan_idx').is_in(ms1['index'].implode())).height == 0


def test_fdr_densities_with_empty_decoys(fd):
    densities = insight_data.fdr_densities(fd, FD_DATASET, 'density_target', 'density_decoy').collect()
    assert set(densities['kind']) == {'target'}
    assert densities.schema['x'] == pl.Float64


def test_proteins_best_per_spectrum(tnt):
    proteins = insight_data.proteins(tnt, TNT_DATASET).collect()
    raw = insight_data.scan(tnt, TNT_DATASET, 'protein_dfs').collect()
    assert proteins['index'].to_list() == raw['index'].to_list()
    best = proteins.filter('is_best')
    assert best.height == raw['Scan'].n_unique()
    for scan_number, group in raw.group_by('Scan'):
        top = group.sort('Score', descending=True, maintain_order=True)['index'][0]
        assert best.filter(pl.col('Scan') == scan_number[0])['index'].to_list() == [top]
    assert proteins['scan_idx'].null_count() == 0


def test_tags_get_scan_index(tnt):
    tags = insight_data.tags(tnt, TNT_DATASET).collect()
    assert tags.height == insight_data.scan(tnt, TNT_DATASET, 'tag_dfs').collect().height
    assert tags['scan_idx'].null_count() == 0


def test_proteoform_sequences_from_pickled_results(tnt):
    # The demo results predate the parquet sequence store
    assert insight_data.result_path(tnt, TNT_DATASET, 'sequence_data').suffix == '.gz'
    rows = insight_data.proteoform_sequences(tnt, TNT_DATASET).collect()
    assert rows.height == insight_data.scan(tnt, TNT_DATASET, 'protein_dfs').collect().height
    row = rows.row(0, named=True)
    assert isinstance(row['sequence'], list) and len(row['sequence']) > 0
    assert len(row['MonoMass']) > 0
    # Named as Insight's ProteinSequenceView reads them
    assert {'start_position', 'end_position', 'observed_mass'} <= set(rows.columns)
    assert not {'proteoform_start', 'proteoform_end', 'computed_mass'} & set(rows.columns)


def test_quant_traces(quant):
    traces = insight_data.quant_traces(quant, 'example').collect()
    raw = insight_data.scan(quant, 'example', 'quant_dfs').collect()
    assert traces.height == raw['Charges'].list.len().sum()
    first = traces.filter((pl.col('FeatureGroupIndex') == raw['FeatureGroupIndex'][0]) & (pl.col('trace') == 0))
    assert first['charge'][0] == raw['Charges'][0][0]
    assert first['mz'][0].to_list() == pytest.approx([float(v) for v in raw['MZs'][0][0].split(',')])
    table = insight_data.quant_table(quant, 'example').collect()
    assert 'Charges' not in table.columns and table.height == raw.height


def test_annotation_summary_fields_renamed():
    lf = pl.LazyFrame({
        'proteoform_end': [3],
        'annotation_summary': [{'fa_score': 1.0, 'coverage': 2.0}],
    })
    out = insight_data._generic_sequence_schema(lf).collect()
    assert out.columns == ['end_position', 'annotation_summary']
    assert out['annotation_summary'][0] == {'score': 1.0, 'coverage': 2.0}
