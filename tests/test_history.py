"""The History of a result: its timeline, kept processing report and changed settings."""

from pathlib import Path

import polars as pl

from src import analysis, history
from src.workflow.FileManager import FileManager
from src.workflow.processing_report import command_io
from src.workflow.progress import RunProgress

INI = """<?xml version="1.0" encoding="ISO-8859-1"?>
<PARAMETERS version="1.7.0">
  <NODE name="FLASHDeconv">
    <ITEM name="version" value="3.6.0" type="string"/>
    <NODE name="1">
      <ITEM name="in" value="" type="input-file"/>
      <ITEM name="threads" value="1" type="int"/>
      <NODE name="SD">
        <ITEM name="tol" value="10.0" type="double"/>
        <ITEM name="min_charge" value="1" type="int"/>
        <ITEMLIST name="isolation" type="double">
          <LISTITEM value="1.0"/>
          <LISTITEM value="2.0"/>
        </ITEMLIST>
      </NODE>
    </NODE>
  </NODE>
</PARAMETERS>
"""


def _file_manager(tmp_path):
    cache = Path(tmp_path, 'cache')
    return FileManager(cache, cache)


def _settings(fm, dataset_id, mzml, identify=False, quant='none'):
    fm.store_data(dataset_id, history.SETTINGS_KEY,
                  {'identify': identify, 'quantification': quant, 'plex': None, 'mzml': mzml})


def _start(progress, tool, inputs, outputs, params):
    step = progress.start_step(tool, [Path(i).name for i in inputs], provenance={
        'kind': 'topp', 'program': tool, 'version': '3.6.0', 'params': params,
        'io': command_io({'in': inputs, 'out': outputs}, len(inputs))})
    for i in range(len(inputs)):
        progress.file_started(step, i)
        progress.file_finished(step, i, True, returncode=0)
    progress.end_step(step)


def test_the_timeline_follows_reused_runs_added_files_and_the_comparison(tmp_path):
    fm = _file_manager(tmp_path)
    a, b, d, c = ('A_20261101-100000', 'B_20261102-100000', 'D_20261103-100000', 'comparison_20261104-100000')
    _settings(fm, a, '/in/a.mzML')
    history.record(fm, a, history.start_record(1000.0, ['/in/a.mzML'], False, 'none', None))
    _settings(fm, b, '/in/a.mzML', identify=True)
    history.record(fm, b, history.start_record(
        2000.0, ['/in/a.mzML'], True, 'none', None,
        follows={'name': 'A', 'note': 'to add identification', 'reuses': [a]}, source=a))
    _settings(fm, d, '/in/d.mzML', identify=True)
    history.record(fm, d, history.start_record(
        3000.0, ['/in/d.mzML'], True, analysis.LABEL_FREE, None, follows={'name': 'B', 'adds_to': [b]}))
    fm.store_data(c, 'quant_manifest', {'kind': 'label-free', 'sources': {'runs': [b, d]}})
    fm.store_data(c, 'quant_design', pl.DataFrame({'condition': ['ctrl', 'drug']}))
    history.record(fm, c, history.start_record(3000.0, [], True, analysis.LABEL_FREE, None))

    events = history.timeline(fm, c)
    assert [e['what'] for e in events] == [history.STARTED, history.IDENTIFICATION,
                                           history.FILES_ADDED, history.COMPARED]
    assert events[0]['files'] == ['a.mzML']
    assert 'reused' in events[1]['detail']
    assert '*ctrl* / *drug*' in events[3]['detail'] and '2 runs' in events[3]['detail']
    # A run's own history stops at its own start
    assert [e['what'] for e in history.timeline(fm, a)] == [history.STARTED]
    assert history.comparisons_of(fm, d) == [c]
    assert '`a.mzML`' in history.event_markdown(events[0])


def test_the_runs_of_one_start_are_one_event(tmp_path):
    fm = _file_manager(tmp_path)
    runs = ['X_20261101-100000', 'Y_20261101-100001']
    for run in runs:
        _settings(fm, run, f'/in/{run[0]}.mzML', quant=analysis.LABEL_FREE)
        history.record(fm, run, history.start_record(500.0, [f'/in/{run[0]}.mzML'], False, analysis.LABEL_FREE, None))
    c = 'comparison_20261101-100500'
    fm.store_data(c, 'quant_manifest', {'kind': 'label-free', 'sources': {'runs': runs}})
    events = history.timeline(fm, c)
    assert [e['what'] for e in events] == [history.STARTED, history.COMPARED]
    assert events[0]['files'] == ['X.mzML', 'Y.mzML']


def test_a_run_made_before_the_history_shows_what_its_settings_tell(tmp_path):
    fm = _file_manager(tmp_path)
    run = 'old_20250102-030405'
    _settings(fm, run, '/in/old.mzML', identify=True)
    fm.store_data(run, 'FD_parameters', {'SD:tol': 5.0, 'SD:min_charge': 1})
    wf = Path(tmp_path, 'flashtnt')
    Path(wf, 'ini').mkdir(parents=True)
    Path(wf, 'ini', 'FLASHDeconv.ini').write_text(INI)

    [event] = history.timeline(fm, run)
    assert event['what'] == history.PROCESSED and event['files'] == ['old.mzML']
    assert event['when'] == history.stamp_time(run)
    assert 'identification' in event['detail']
    rows = history.changed_settings(fm, run, wf)
    assert [(r[0], r[2], r[3]) for r in rows] == [('FLASHDeconv', '10.0', '5.0')]
    assert not history.has_report(fm, run)


def test_ini_defaults_read_the_tool_parameters(tmp_path):
    ini = Path(tmp_path, 'FLASHDeconv.ini')
    ini.write_text(INI)
    defaults = history.ini_defaults(ini)
    assert defaults['SD:tol'] == '10.0'
    assert defaults['SD:isolation'] == '1.0, 2.0'
    assert history.ini_defaults(Path(tmp_path, 'missing.ini')) == {}


def test_each_result_keeps_the_report_of_the_start_that_made_it(tmp_path):
    fm = _file_manager(tmp_path)
    wf = Path(tmp_path, 'flashtnt')
    Path(wf, 'ini').mkdir(parents=True)
    Path(wf, 'ini', 'FLASHDeconv.ini').write_text(INI)
    _start(RunProgress(wf), 'FLASHDeconv', ['/in/a.mzML'], ['/out/a_deconv.mzML'],
           {'SD:tol': '5.0', 'SD:min_charge': '1', 'threads': '4'})
    run = 'A_20261101-100000'
    _settings(fm, run, '/in/a.mzML')
    history.keep(wf, fm, [run])
    assert history.has_report(fm, run)
    assert fm.result_exists(run, history.REPORT_TAG)
    assert 'FLASHDeconv' in Path(fm.get_result_path(run, history.REPORT_TAG)).read_text(encoding='utf-8')

    # The next start writes over the logs; the kept report stays
    Path(wf, 'logs', 'progress.json').unlink()
    _start(RunProgress(wf), 'FLASHDeconv', ['/in/b.mzML'], ['/out/b_deconv.mzML'], {'SD:tol': '20.0'})
    rows = history.changed_settings(fm, run, wf)
    assert [(r[0], r[2], r[3]) for r in rows] == [('FLASHDeconv', '10.0', '5.0')]


def test_nothing_is_kept_without_a_recorded_start(tmp_path):
    fm = _file_manager(tmp_path)
    history.keep(Path(tmp_path, 'flashtnt'), fm, ['A_20261101-100000'])
    assert not history.has_report(fm, 'A_20261101-100000')
