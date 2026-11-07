"""Which presets fit the chosen analysis (src/analysis.py)."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import analysis  # noqa: E402

PRESETS = json.loads(Path(__file__).resolve().parents[1].joinpath('presets.json').read_text())['flashtnt']


def test_sample_presets_fit_only_with_identification():
    for name in ('Native', 'Lysate'):
        assert analysis.preset_fit(PRESETS[name], True, analysis.NO_QUANT) == []
        assert analysis.preset_fit(PRESETS[name], False, analysis.NO_QUANT) == [analysis.IDENTIFICATION]
        assert analysis.preset_fit(PRESETS[name], False, analysis.LABEL_FREE) == [analysis.IDENTIFICATION]


def test_quantification_presets_need_their_quantification():
    label_free = {'FLASHQuantConsensus': {'rt_tol': 60.0}}
    labels = {'FLASHDeconv': {'iq:isotope_correction': 'false'}, '_general': {analysis.PLEX_KEY: 'tmt10plex'}}
    assert analysis.preset_fit(label_free, True, analysis.NO_QUANT) == [analysis.LABEL_FREE_STEP]
    assert analysis.preset_fit(label_free, False, analysis.LABEL_FREE) == []
    assert analysis.preset_fit(labels, True, analysis.LABEL_FREE) == [analysis.ISOBARIC_STEP]
    assert analysis.preset_fit(labels, True, analysis.ISOBARIC) == []


def test_deconvolution_settings_always_fit():
    preset = {'_description': 'x', 'FLASHDeconv': {'SD:min_charge': 2}, '_general': {'max_threads': 2}}
    assert analysis.preset_fit(preset, False, analysis.NO_QUANT) == []
    assert analysis.preset_warning('X', preset, False, analysis.NO_QUANT) is None


def test_warning_names_the_switch():
    text = analysis.preset_warning('Native', PRESETS['Native'], False, analysis.NO_QUANT)
    assert 'identification' in text and 'Identify proteoforms' in text
    text = analysis.preset_warning('Q', {'FLASHQuantConsensus': {'a': 1}}, True, analysis.ISOBARIC)
    assert 'Label-free' in text


def test_sample_presets_support_analyses_without_identification():
    # Deconvolution only (DDA runs) and label-free without identification are
    # supported: the sample type is a note there, never a warning
    for key in ('Native', 'Lysate'):
        assert analysis.preset_note(key, PRESETS[key], True, analysis.LABEL_FREE) is None
        for quant in (analysis.NO_QUANT, analysis.LABEL_FREE, analysis.ISOBARIC):
            level, text = analysis.preset_note(key, PRESETS[key], False, quant)
            assert level == 'info'
            assert 'Identify proteoforms' in text


def test_other_presets_still_warn():
    level, _ = analysis.preset_note('Q', {'FLASHQuantConsensus': {'a': 1}}, False, analysis.NO_QUANT)
    assert level == 'warning'
    level, _ = analysis.preset_note('Native', {'FLASHQuantConsensus': {'a': 1}, **PRESETS['Native']},
                                    False, analysis.NO_QUANT)
    assert level == 'warning'


def test_internal_fragments_follow_the_sample_type():
    # A lysate is annotated from the termini only; native top-down with internal fragments too
    assert PRESETS['Lysate']['FLASHAnnotator']['algorithm:internal_ions'] == 'false'
    assert PRESETS['Native']['FLASHAnnotator']['algorithm:internal_ions'] == 'true'
