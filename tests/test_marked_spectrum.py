"""The raw heatmaps mark the peaks of the selected spectrum."""

import pytest

pytest.importorskip("openms_insight")

from src.render import render  # noqa: E402


def test_ms2_map_marks_the_selected_spectrum():
    assert render._marked_spectrum('ms2_raw_heatmap', None, 'd', 'topdown', object()) == 'scan'
    assert render._marked_spectrum('ms2_raw_heatmap', None, 'd', 'flashtnt', object()) == 'pscan'


def test_other_maps_and_no_state_mark_nothing():
    assert render._marked_spectrum('ms1_deconv_heat_map', None, 'd', 'topdown', object()) is None
    assert render._marked_spectrum('ms1_raw_heatmap', None, 'd', 'topdown', None) is None


def test_ms1_map_marks_the_precursor_spectrum(monkeypatch):
    followed = []
    monkeypatch.setattr(render.insight_data, 'result_path', lambda fm, ds, name: 'x')
    monkeypatch.setattr(render, '_follow_selected_spectrum',
                        lambda fm, ds, sm, identifier: followed.append(identifier))
    assert render._marked_spectrum('ms1_raw_heatmap', None, 'd', 'topdown', object()) == render.SURVEY
    assert followed == ['scan']


def test_ms1_map_without_ms1_spectra_marks_nothing(monkeypatch):
    def no_ms1(*args):
        raise render.MissingResult('no MS1 spectra')
    monkeypatch.setattr(render.insight_data, 'result_path', lambda fm, ds, name: 'x')
    monkeypatch.setattr(render, '_follow_selected_spectrum', no_ms1)
    assert render._marked_spectrum('ms1_raw_heatmap', None, 'd', 'topdown', object()) is None
    monkeypatch.setattr(render.insight_data, 'result_path', lambda fm, ds, name: None)
    assert render._marked_spectrum('ms1_raw_heatmap', None, 'd', 'topdown', object()) is None
