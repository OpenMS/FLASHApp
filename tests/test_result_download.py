"""The ZIP archives of the Download page and the Viewer's Download button."""

import os
import zipfile
from io import BytesIO
from pathlib import Path

from src.common import result_download as rd
from src.workflow.FileManager import FileManager


def _manager(tmp_path):
    fm = FileManager(tmp_path, Path(tmp_path, 'cache'))
    for run in ('run_a', 'run_b'):
        for tag, name in (('out_tsv', 'out.tsv'), ('spec1_tsv', 'spec1.tsv')):
            fm.store_file(run, tag, BytesIO(f'{run} {tag}'.encode()), file_name=name)
    return fm


def test_runs_and_files(tmp_path):
    fm = _manager(tmp_path)
    assert sorted(rd.downloadable_runs(fm)) == ['run_a', 'run_b']
    assert sorted(f.name for f in rd.result_files(fm, 'run_a')) == ['out.tsv', 'spec1.tsv']


def test_archive_is_reused_and_rebuilt_when_a_file_is_newer(tmp_path):
    fm = _manager(tmp_path)
    files = rd.result_files(fm, 'run_a')
    first = rd.build_archive(files, rd.archive_path(fm, 'run_a'))
    with zipfile.ZipFile(first) as zf:
        assert sorted(zf.namelist()) == ['out.tsv', 'spec1.tsv']
    built = first.stat().st_mtime_ns
    assert rd.build_archive(files, first).stat().st_mtime_ns == built
    # The archive is not a result file of its own
    assert rd.ARCHIVE_NAME not in [f.name for f in rd.result_files(fm, 'run_a')]

    newer = next(f for f in rd.result_files(fm, 'run_a') if f.name == 'out.tsv')
    newer.write_text('changed')
    os.utime(newer, (first.stat().st_mtime + 10,) * 2)
    with zipfile.ZipFile(rd.build_archive(files, first)) as zf:
        assert zf.read('out.tsv') == b'changed'


def test_all_runs_archive_has_a_folder_per_run(tmp_path):
    fm = _manager(tmp_path)
    fm.rename_dataset('run_b', 'My run/2')
    with zipfile.ZipFile(BytesIO(rd.zip_folders(rd.run_folders(fm, ['run_a', 'run_b'])))) as zf:
        assert sorted(zf.namelist()) == [
            'My run_2/out.tsv', 'My run_2/spec1.tsv', 'run_a/out.tsv', 'run_a/spec1.tsv',
        ]


def test_format_size():
    assert rd.format_size(512) == '512 B'
    assert rd.format_size(2048) == '2.0 KB'
    assert rd.format_size(3 * 1024 ** 3) == '3.0 GB'


def test_runs_of_the_same_name_get_their_own_folder(tmp_path):
    fm = _manager(tmp_path)
    fm.rename_dataset('run_a', 'Same')
    fm.rename_dataset('run_b', 'Same')
    assert list(rd.run_folders(fm, ['run_a', 'run_b'])) == ['Same', 'Same (2)']


def test_safe_name():
    assert rd.safe_name('Example: lysate, identified (E. coli / spike-in)') == \
        'Example - lysate, identified (E. coli _ spike-in)'
    assert rd.safe_name('..') == 'results'


def test_compressed_files_are_stored_and_text_is_deflated(tmp_path):
    text = Path(tmp_path, 'out.tsv')
    text.write_text('mass\tintensity\n' * 1000)
    cache = Path(tmp_path, 'mass_table.pq')
    cache.write_bytes(os.urandom(4096))
    target = rd.build_archive([text, cache], Path(tmp_path, 'results.zip'))
    with zipfile.ZipFile(target) as zf:
        assert zf.getinfo('out.tsv').compress_type == zipfile.ZIP_DEFLATED
        assert zf.getinfo('mass_table.pq').compress_type == zipfile.ZIP_STORED
        assert zf.read('out.tsv') == text.read_bytes()
        assert zf.read('mass_table.pq') == cache.read_bytes()
