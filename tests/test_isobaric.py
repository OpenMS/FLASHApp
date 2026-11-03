"""Isobaric labelling: FLASHDeconv's reporter ion table, the channel table and
the comparison levels built from them (src/quant/isobaric.py), and what the
analysis choice on Describe experiment means (src/analysis.py)."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pathlib import Path  # noqa: E402

import polars as pl  # noqa: E402

from src import analysis  # noqa: E402
from src.quant import contract, isobaric  # noqa: E402

HEADER = ('ScanNum\tPrecursorScanNum\tPrecursorFeatureIndex\tPrecursorMZ\tPrecursorCharge\t'
          'PrecursorMonoisotopicMass\tPrecursorQscore\tRetentionTime\tQuantifiedChannelCount\t'
          + '\t'.join(f'QuantForCh{i}' for i in range(1, 7)))


def _quant_file(tmp_path, rows):
    path = Path(tmp_path, 'quant.tsv')
    path.write_text(HEADER + '\n' + ''.join('\t'.join(map(str, r)) + '\n' for r in rows))
    return path


def _row(scan, mass, channels, rt=120.0):
    mz = 1001.0 if isinstance(mass, str) else mass / 10 + 1
    return [scan, scan - 1, 0, mz, 10, mass, 0.9, rt, len(channels), *channels]


def test_read_quant_names_the_channels_and_gives_minutes(tmp_path):
    path = _quant_file(tmp_path, [_row(2, 10000.0, [1, 2, 3, 4, 5, 6])])
    frame = isobaric.read_quant(path, 'tmt6plex')
    assert frame.columns[:9] == ['Scan', 'PrecursorScan', 'Feature', 'PrecursorMZ', 'Charge',
                                 'PrecursorMass', 'Qscore', 'RT', 'QuantifiedChannels']
    assert frame.columns[9:] == isobaric.channel_names('tmt6plex')
    assert frame['RT'][0] == 2.0 and frame['131'][0] == 6.0


def test_read_quant_of_an_empty_or_missing_file(tmp_path):
    empty = Path(tmp_path, 'empty.tsv')
    empty.write_text('')
    for path in (empty, Path(tmp_path, 'missing.tsv'), None):
        frame = isobaric.read_quant(path, 'tmt10plex')
        assert frame.is_empty() and '131' in frame.columns


def test_channel_table_keeps_known_channels_and_unique_samples():
    table = pl.DataFrame({
        'channel': ['126', '127', 'unknown'], 'sample': ['A', 'A', 'B'], 'condition': ['ctrl', ' x ', 'y'],
    })
    aligned = isobaric.align_channels(table, 'tmt6plex')
    assert aligned['channel'].to_list() == isobaric.channel_names('tmt6plex')
    assert aligned['sample'].to_list()[:3] == ['A (126)', 'A (127)', '128']
    assert aligned['condition'].to_list()[:3] == ['ctrl', 'x', '']


def test_channel_table_round_trip(tmp_path):
    table = isobaric.default_channels('itraq4plex').with_columns(
        condition=pl.Series(['a', 'a', 'b', 'b']))
    isobaric.save_channels(tmp_path, table, 'itraq4plex')
    loaded = isobaric.load_channels(tmp_path, 'itraq4plex')
    assert loaded['condition'].to_list() == ['a', 'a', 'b', 'b']
    # Another plex falls back to its own channels
    assert isobaric.load_channels(tmp_path, 'tmt6plex')['channel'].to_list() == isobaric.channel_names('tmt6plex')
    assert isobaric.design(loaded) == {'114': ('a', 1), '115': ('a', 2), '116': ('b', 1), '117': ('b', 2)}


def _prsms(tmp_path, rows):
    path = Path(tmp_path, 'prsms.tsv')
    header = 'Scan\tProteinAccession\tProteinDescription\tProforma\tProteoformMass\tProteoformLevelQvalue\tScore\n'
    path.write_text(header + ''.join('\t'.join(map(str, r)) + '\n' for r in rows))
    return isobaric.read_prsms(path)


def test_read_prsms_keeps_the_best_target_at_one_percent(tmp_path):
    prsms = _prsms(tmp_path, [
        (2, 'sp|P1|A', 'A', 'PEPTIDE', 10000.0, 0.0, 10),
        (2, 'sp|P2|B', 'B', 'PEPTIDES', 10087.0, 0.0, 20),
        (3, 'DECOY_sp|P1|A', 'A', 'X', 10000.0, 0.0, 30),
        (4, 'sp|P3|C', 'C', 'Y', 11000.0, 0.5, 30),
    ])
    assert prsms['Scan'].to_list() == [2]
    assert prsms['accession'][0] == 'sp|P2|B' and prsms['proteoform'][0] == 'sp|P2|B|PEPTIDES'


def test_build_links_spectra_across_runs_and_sums_the_channels(tmp_path):
    channels = isobaric.default_channels('tmt6plex').with_columns(
        condition=pl.Series(['c'] * 3 + ['t'] * 3))
    run1 = isobaric.read_quant(_quant_file(tmp_path, [
        _row(2, 10000.0, [1, 1, 1, 2, 2, 2]),
        _row(3, 10000.05, [1, 1, 1, 2, 2, 2]),   # 5 ppm: the same feature
        _row(4, 12000.0, [0, 0, 0, 0, 0, 0]),    # no reporter ions: left out
    ]), 'tmt6plex')
    run2 = isobaric.read_quant(_quant_file(tmp_path, [_row(2, 10000.02, [1, 1, 1, 2, 2, 2])]), 'tmt6plex')
    prsms = _prsms(tmp_path, [(2, 'sp|P1|A', 'A', 'PEPTIDE', 10000.0, 0.0, 10)])
    levels = isobaric.build([('r1', run1, prsms), ('r2', run2, None)], channels)
    assert set(levels) == {'feature', 'proteoform', 'protein'}
    long, info = levels['feature']
    assert long.columns == list(contract.LONG_SCHEMA)
    assert info.height == 1 and info['Spectra'][0] == 3 and info['Runs'][0] == 2
    assert dict(zip(long['sample'], long['intensity'])) == {'126': 3, '127': 3, '128': 3,
                                                            '129': 6, '130': 6, '131': 6}
    form_long, forms = levels['proteoform']
    assert forms['PrSMs'][0] == 1 and form_long.filter(pl.col('sample') == '131')['intensity'][0] == 2


def test_build_uses_the_identified_mass_without_a_precursor_mass(tmp_path):
    channels = isobaric.default_channels('tmt6plex')
    quant = isobaric.read_quant(_quant_file(tmp_path, [_row(2, 'NaN', [1, 1, 1, 1, 1, 1])]), 'tmt6plex')
    assert isobaric.build([('r', quant, None)], channels) == {}
    prsms = _prsms(tmp_path, [(2, 'sp|P1|A', 'A', 'PEPTIDE', 10000.0, 0.0, 10)])
    levels = isobaric.build([('r', quant, prsms)], channels)
    assert levels['feature'][1]['mass'][0] == 10000.0


def test_quantification_choice_and_its_former_switch():
    assert analysis.quantification({}) == analysis.NO_QUANT
    assert analysis.quantification({analysis.COMPARE_KEY: True}) == analysis.LABEL_FREE
    assert analysis.quantification({analysis.QUANT_KEY: analysis.ISOBARIC, analysis.COMPARE_KEY: True}) == analysis.ISOBARIC
    assert analysis.plex({analysis.PLEX_KEY: 'nonsense'}) == isobaric.DEFAULT_PLEX


def test_analysis_names_every_mode():
    assert analysis.name(False, analysis.NO_QUANT) == 'Deconvolution'
    assert analysis.name(True, analysis.LABEL_FREE) == 'Deconvolution + identification + label-free quantification'
    assert analysis.name(False, analysis.ISOBARIC, 'tmt16plex') == 'Deconvolution + TMTpro 16-plex quantification'


def test_summary_says_what_a_start_does():
    assert analysis.summary(0, True, 'Native', None, analysis.NO_QUANT) == 'Select the MS data to analyse.'
    text = analysis.summary(1, True, 'Lysate', 'x/db.fasta', analysis.LABEL_FREE)
    assert '**db.fasta**' in text and 'lysate settings' in text and 'at least two runs' in text
    text = analysis.summary(3, False, 'Native', None, analysis.ISOBARIC, 'tmt10plex')
    assert 'FLASHTnT' not in text and 'TMT 10-plex' in text and 'fractions' in text


def test_labelled_runs_get_the_reporter_ions_in_their_default():
    from src.render import topdown_layout
    from src.render.layout_editor.catalog import (
        DECONV_DEFAULT_LAYOUT, TNT_DEFAULT_LAYOUTS, TOPDOWN_COMPONENTS, TOPDOWN_SEQUENCE_COMPONENTS,
        validate_layout,
    )
    catalog = TOPDOWN_COMPONENTS + TOPDOWN_SEQUENCE_COMPONENTS
    for layout, after in [(DECONV_DEFAULT_LAYOUT, 'scan_table'), *[(l, None) for l in TNT_DEFAULT_LAYOUTS.values()]]:
        rows = topdown_layout.with_reporters(layout, topdown_layout.REPORTER_AFTER)
        assert list(topdown_layout.REPORTER_ROW) in rows
        assert validate_layout([rows], catalog) == ''
        if after:
            assert rows[rows.index(list(topdown_layout.REPORTER_ROW)) - 1].count(after) == 1
        # Never twice
        assert topdown_layout.with_reporters(rows, topdown_layout.REPORTER_AFTER) == rows


def test_a_proteoform_without_any_mass_is_its_own_feature(tmp_path):
    channels = isobaric.default_channels('tmt6plex')
    quant = isobaric.read_quant(_quant_file(tmp_path, [_row(1, 'NaN', [1, 1, 1, 2, 2, 2])]), 'tmt6plex')
    prsms = _prsms(tmp_path, [(1, 'sp|P1|A', 'A', 'PEPTIDE', -1, 0.0, 10)])
    assert prsms['mass'][0] is None
    levels = isobaric.build([('r1', quant, prsms), ('r2', quant, prsms)], channels)
    info = levels['feature'][1]
    assert info.height == 1 and info['Runs'][0] == 2 and info['label'][0].startswith('PEPTIDE')
    assert levels['proteoform'][1]['label'][0] == 'A PEPTIDE'
