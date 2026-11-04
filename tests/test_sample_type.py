"""Native or lysate recommendation (src/sample_type.py)."""

from pathlib import Path

import pytest

from src import sample_type
from src.sample_type import LYSATE, NATIVE, RunSpectra, classify, count_spectra, recommend

TNT_EXAMPLES = Path('example-data/workspaces/default/flashtnt/input-files/mzML-files')


def _mzml(path, levels):
    spectra = ''.join(
        f'<spectrum index="{i}" id="scan={i + 1}" defaultArrayLength="0">'
        f'<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="{level}"/>'
        '</spectrum>'
        for i, level in enumerate(levels)
    )
    path.write_text(f'<mzML><run><spectrumList count="{len(levels)}">{spectra}</spectrumList></run></mzML>')
    return path


def test_classify():
    assert classify(0, 500) == NATIVE  # no MS1 survey
    assert classify(3, 1) == NATIVE  # one targeted spectrum
    assert classify(200, 2000) == LYSATE


def test_count_spectra_and_early_stop(tmp_path, monkeypatch):
    run = _mzml(tmp_path / 'run.mzML', [1, 2, 2, 1, 2])
    assert (count_spectra(run).ms1, count_spectra(run).ms2) == (2, 3)
    # Tags split across read chunks are still counted
    monkeypatch.setattr(sample_type, '_CHUNK', 7)
    assert (count_spectra(run).ms1, count_spectra(run).ms2) == (2, 3)
    lysate = _mzml(tmp_path / 'lysate.mzML', ([1] + [2] * 5) * 400)
    counted = count_spectra(lysate)
    assert counted.sample_type == LYSATE and not counted.complete


def test_recommend(tmp_path):
    native = _mzml(tmp_path / 'a.mzML', [2])
    lysate = _mzml(tmp_path / 'b.mzML', ([1] + [2] * 5) * 10)
    assert recommend([native]).sample_type == NATIVE
    assert 'no MS1 survey scans' in recommend([native]).reason()
    assert recommend([lysate]).sample_type == LYSATE
    mixed = recommend([native, lysate])
    assert mixed.sample_type is None and 'separately' in mixed.reason()
    # A .raw file is counted from its converted mzML, if there is one
    raw = tmp_path / 'c.raw'
    raw.write_bytes(b'')
    assert recommend([raw]).unreadable == [str(raw)]
    converted = tmp_path / 'converted'
    converted.mkdir()
    _mzml(converted / 'c.mzML', [2])
    assert recommend([raw], converted).sample_type == NATIVE


@pytest.mark.skipif(not TNT_EXAMPLES.exists(), reason='example data not checked out')
def test_the_examples_are_named_after_their_sample_type():
    native = TNT_EXAMPLES / 'Native-identified_GapA_HCD.mzML'
    lc_ms = TNT_EXAMPLES / 'Deconvolution-only_ProteinMix_LC-MS.mzML'
    assert recommend([native]).sample_type == NATIVE
    assert recommend([lc_ms]).sample_type == LYSATE


class _Results:
    def __init__(self, settings, path):
        self.settings, self.path = settings, path

    def result_exists(self, dataset_id, name):
        return self.settings is not None

    def get_results(self, dataset_id, names):
        return {'FTnT_parameters': self.settings}

    def get_result_path(self, dataset_id, name):
        return self.path


def test_of_result(tmp_path):
    import polars as pl
    path = tmp_path / 'scan_table.pq'
    pl.DataFrame({'MSLevel': [1, 2, 2] * 20}).write_parquet(path)
    assert sample_type.of_result(_Results({}, path), 'd') == LYSATE
    # Searched with -native: native, whatever the spectra
    assert sample_type.of_result(_Results({'native': True}, path), 'd') == NATIVE
    pl.DataFrame({'MSLevel': [2]}).write_parquet(path)
    # Searched without -native: run as a lysate, whatever the spectra
    assert sample_type.of_result(_Results({}, path), 'd') == LYSATE
    # No parameters (a manual upload): from the spectra
    assert sample_type.of_result(_Results(None, path), 'd') == NATIVE
    assert sample_type.of_result(_Results(None, None), 'd') is None


def test_a_flag_stored_as_text_is_read_as_its_value():
    from src.sample_type import is_on
    assert is_on(True) and is_on('true') and is_on('True')
    assert not is_on('false') and not is_on(False) and not is_on(None) and not is_on('')
