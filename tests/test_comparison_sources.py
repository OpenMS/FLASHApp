"""What a comparison was made from (src/quant/sources.py): the stored runs,
the PrSMs of every run, the reporter ions and the consensus features."""

import math

import polars as pl

from src.quant import sources


class _Results:
    """The part of FileManager the sources read."""

    def __init__(self, settings):
        self.settings = settings

    def result_exists(self, dataset_id, name):
        return name == sources.SETTINGS_KEY and dataset_id in self.settings

    def get_results(self, dataset_id, names):
        return {sources.SETTINGS_KEY: self.settings[dataset_id]}

    def get_results_list(self, names):
        return list(self.settings)


def test_runs_are_found_by_id_and_by_their_mzml():
    fm = _Results({
        'A_20260101-101010': {'mzml': '/x/A.mzML'},
        'example_b': {'mzml': '/y/B.mzML'},
        'older_b': {'mzml': '/y/B.mzML'},
    })
    manifest = {'runs': ['A', 'B', 'C'], 'sources': {'runs': ['A_20260101-101010', 'B_20260101-111111']}}
    assert sources.run_datasets(fm, manifest) == {'A': 'A_20260101-101010', 'B': 'older_b'}


def test_the_fractions_of_an_isobaric_comparison_are_its_runs():
    fm = _Results({'example_tmt': {'mzml': 'Set1.mzML'}})
    manifest = {'kind': 'isobaric', 'runs': ['s1', 's2'], 'fractions': ['Set1_20260101-101010'],
                'sources': {'runs': ['Set1_20260101-101010']}}
    assert sources.run_datasets(fm, manifest) == {'Set1': 'example_tmt'}


def _identifications():
    return pl.DataFrame({
        'ConsensusIndex': [3, 3, 5], 'Run': ['R1', 'R2', 'R1'], 'Scan': [10, 12, 20],
        'RetentionTime': [600.0, 660.0, 1200.0], 'ProteinAccession': ['P1', 'P1', 'P2'],
        'Proforma': ['PEP', 'PEP', 'TIDE'], 'ProteoformMass': [100.0, 100.0, -1.0],
        'ProteoformLevelQvalue': [0.0, 0.001, 0.002], 'proteoform': ['P1|PEP', 'P1|PEP', 'P2|TIDE'],
    })


def test_prsms_are_keyed_by_the_level_shown():
    for level, keys in (('proteoform', ['P1|PEP', 'P1|PEP', 'P2|TIDE']), ('protein', ['P1', 'P1', 'P2']),
                        ('feature', ['3', '3', '5'])):
        frame = sources.prsms(_identifications(), level)
        assert frame['key'].to_list() == keys
        assert frame['prsm_row'].to_list() == [0, 1, 2]
    frame = sources.prsms(_identifications(), 'proteoform')
    assert frame['RT'].to_list() == [10.0, 11.0, 20.0]
    assert frame['ProteoformMass'].to_list() == [100.0, 100.0, None]
    assert frame['PrecursorMass'].null_count() == 3
    assert sources.prsms(None, 'proteoform').height == 0


def test_a_prsm_is_found_by_its_scan_in_the_run():
    hits = pl.DataFrame({'index': [0, 1, 2], 'Scan': [10, 10, 20], 'scan_idx': [4, 4, 9],
                         'is_best': [False, True, True]})
    assert sources.prsm_in_run(hits, 10) == (1, 4)
    assert sources.prsm_in_run(hits, 20) == (2, 9)
    assert sources.prsm_in_run(hits, 30) is None


def _ions():
    frames = {'F1': pl.DataFrame({'Scan': [1, 2, 3], '126': [100.0, 0.0, 50.0], '127': [200.0, 0.0, 0.0],
                                  'Charge': [5, 6, 7]})}
    return sources.reporter_ions(frames, ['126', '127'], ['a', 'b'])


def test_reporter_ions_are_one_row_per_spectrum_and_label():
    ions = _ions()
    assert ions.height == 6
    assert ions.filter(pl.col('channel') == '127')['sample'].unique().to_list() == ['b']
    assert ions['intensity'].null_count() == 3
    labels = sources.labels_found(ions, ['126', '127'])
    assert labels['spectra'].to_list() == [1, 1, 1]


def test_channel_presence_counts_the_labelled_spectra():
    presence = sources.channel_presence(_ions(), ['126', '127'])
    assert presence['channel'].to_list() == ['126', '127']
    assert presence['present'].to_list() == [100.0, 50.0]
    assert abs(presence['median'][1] - math.log2(200.0)) < 1e-9


def test_reporter_ratios_need_both_ions():
    ratios = sources.reporter_ratios(_ions(), '126')
    assert ratios['log2_ratio'].to_list() == [1.0]


def test_feature_overview_counts_how_each_run_found_a_feature():
    long = pl.DataFrame({
        'id': ['1', '1', '1', '2', '2', '2'], 'sample': ['A', 'B', 'C'] * 2,
        'intensity': [8.0, 4.0, None, 2.0, None, None],
        'status': ['detected', 'transferred', 'missing', 'detected', 'missing', 'missing'],
    })
    annotations = pl.DataFrame({'id': ['1', '2'], 'label': ['f1', 'f2'], 'mass': [1000.0, 900.0],
                                'rt': [10.0, 20.0]})
    overview = sources.feature_overview(long, annotations, ['A', 'B', 'C'])
    assert overview.select('id', 'Detected', 'Matched', 'Missing').rows() == [('1', 1, 1, 1), ('2', 1, 0, 2)]
    assert overview['RowIndex'].to_list() == [0, 1]
    assert overview['proteoform'].null_count() == 2
    values = sources.feature_values(long, overview, ['A', 'B', 'C'])
    first = values.filter(pl.col('RowIndex') == 0)
    assert first['Found'].to_list() == ['detected', 'matched', 'missing']
    assert first['log2_intensity'].to_list() == [3.0, 2.0, None]
