"""The demo workspace ships one example per default layout of the Viewer
(example-data/build_examples.py), each under a name that says what it shows."""

from pathlib import Path

import polars as pl

from src import sample_type
from src.quant.session import QuantSession, comparisons
from src.render.quant_view import comparison_names
from src.workflow.FileManager import FileManager

DEMO = Path('example-data', 'workspaces', 'default')


def _file_manager(tmp_path):
    return FileManager(tmp_path, DEMO / 'cache')


def test_every_example_is_named_after_what_it_shows(tmp_path):
    file_manager = _file_manager(tmp_path)
    runs = file_manager.get_results_list(['threedim_SN_plot'])
    names = {dataset_id: file_manager.get_display_name(dataset_id) for dataset_id in runs}
    assert names == {
        'example_deconvolution': 'Example: deconvolution only (protein mix, LC-MS)',
        'example_native_identified': 'Example: native, identified (GapA)',
        'example_lysate_identified': 'Example: lysate, identified (E. coli + spike-in)',
        'example_isobaric': 'Example: TMT 6-plex, deconvolution only (E. coli, ratio mix)',
    }


def test_identified_examples_open_their_sample_type_layout(tmp_path):
    file_manager = _file_manager(tmp_path)
    assert not file_manager.result_exists('example_deconvolution', 'FTnT_parameters')
    assert sample_type.of_result(file_manager, 'example_native_identified') == sample_type.NATIVE
    assert sample_type.of_result(file_manager, 'example_lysate_identified') == sample_type.LYSATE


def test_the_label_free_example_compares_two_conditions(tmp_path):
    file_manager = _file_manager(tmp_path)
    assert sorted(comparisons(file_manager)) == ['example_isobaric_comparison', 'example_label_free']
    assert comparison_names(file_manager) == {
        'Example: label-free comparison (spike-in 1x vs 0.1x)': 'example_label_free',
        'Example: TMT comparison (8x vs 4x vs 1x)': 'example_isobaric_comparison',
    }
    session = QuantSession(file_manager, 'example_label_free')
    assert not session.isobaric
    assert len(session.runs) == 4
    design = session._get('quant_design')
    assert sorted(set(design['condition'])) == ['Spike-in 0.1x', 'Spike-in 1x']


def test_the_tmt_example_compares_its_reporter_channels(tmp_path):
    file_manager = _file_manager(tmp_path)
    assert file_manager.get_results('example_isobaric', ['isobaric_plex'])['isobaric_plex']['plex'] == 'tmt6plex'
    session = QuantSession(file_manager, 'example_isobaric_comparison')
    assert session.isobaric
    assert session.conditions() == {
        '8x_a': '8x', '4x_a': '4x', '1x_a': '1x', '1x_b': '1x', '4x_b': '4x', '8x_b': '8x',
    }
    # The channels hold 8:4:1 of the same lysate: normalizing them would
    # remove the difference the example shows
    assert session.settings['postproc-normalization'] == 'None'


def test_the_reporter_plot_shows_room_beyond_the_outermost_labels(tmp_path):
    from src.quant import isobaric
    from src.render import insight_data
    file_manager = _file_manager(tmp_path)
    channels = isobaric.channel_names('tmt6plex')
    mzs = dict(isobaric.PLEXES['tmt6plex'][1])
    ions = insight_data.reporter_ions(
        file_manager, 'example_isobaric', channels, [mzs[c] for c in channels],
    ).collect()
    first = ions.filter(pl.col('scan_idx') == ions['scan_idx'][0])
    labelled = first.filter(pl.col('label') != '')
    assert labelled.height == 6
    assert first['mz'].min() <= labelled['mz'].min() - insight_data.REPORTER_MARGIN + 1e-9
    assert first['mz'].max() >= labelled['mz'].max() + insight_data.REPORTER_MARGIN - 1e-9
