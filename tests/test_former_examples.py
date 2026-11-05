"""Example runs of earlier releases are taken out of an existing workspace,
so the upload can show today's examples."""

import hashlib

from src import former_examples


def _former(monkeypatch, *contents):
    monkeypatch.setattr(
        former_examples, 'FORMER_EXAMPLES',
        {(len(c), hashlib.sha256(c).hexdigest()) for c in contents})


def test_retire_removes_earlier_examples_only(tmp_path, monkeypatch):
    _former(monkeypatch, b'old example run')
    shipped = tmp_path / 'shipped'
    shipped.mkdir()
    (shipped / 'New_example.mzML').write_bytes(b'new example run')
    folder = tmp_path / 'mzML-files'
    folder.mkdir()
    (folder / 'Old_example.mzML').write_bytes(b'old example run')
    (folder / 'Renamed.mzML').write_bytes(b'new example run')  # renamed since
    (folder / 'New_example.mzML').write_bytes(b'new example run')
    (folder / 'mine.mzML').write_bytes(b'a run of my own')
    removed = former_examples.retire(folder, [str(shipped / 'New_example.mzML')], tmp_path / 'marker')
    assert sorted(removed) == ['Old_example.mzML', 'Renamed.mzML']
    assert sorted(p.name for p in folder.iterdir() ) == [
        'New_example.mzML', 'mine.mzML']


def test_retire_runs_once(tmp_path, monkeypatch):
    _former(monkeypatch, b'old example run')
    folder = tmp_path / 'mzML-files'
    folder.mkdir()
    former_examples.retire(folder, [], tmp_path / 'marker')
    (folder / 'Old_example.mzML').write_bytes(b'old example run')
    assert former_examples.retire(folder, [], tmp_path / 'marker') == []
    assert (folder / 'Old_example.mzML').exists()


def test_retire_drops_referenced_examples(tmp_path, monkeypatch):
    _former(monkeypatch, b'old example run')
    elsewhere = tmp_path / 'flashdeconv'
    elsewhere.mkdir()
    (elsewhere / 'example_fd.mzML').write_bytes(b'old example run')
    (elsewhere / 'mine.mzML').write_bytes(b'a run of my own')
    folder = tmp_path / 'mzML-files'
    folder.mkdir()
    (folder / 'external_files.txt').write_text(
        f"{elsewhere / 'example_fd.mzML'}\n{elsewhere / 'mine.mzML'}\n")
    assert former_examples.retire(folder, [], tmp_path / 'marker') == ['example_fd.mzML']
    assert (folder / 'external_files.txt').read_text() == f"{elsewhere / 'mine.mzML'}\n"
    assert (elsewhere / 'example_fd.mzML').exists()  # only the reference goes


def test_same_size_other_content_is_kept(tmp_path, monkeypatch):
    _former(monkeypatch, b'old example run')
    folder = tmp_path / 'mzML-files'
    folder.mkdir()
    (folder / 'mine.mzML').write_bytes(b'old example XXX')
    assert former_examples.retire(folder, [], tmp_path / 'marker') == []


def _cache(path, rows):
    import sqlite3
    (path / 'files').mkdir(parents=True)
    con = sqlite3.connect(path / 'cache.db')
    con.execute('CREATE TABLE stored_data (id TEXT PRIMARY KEY, display_name TEXT, scan_table TEXT)')
    con.execute('CREATE TABLE stored_files (id TEXT PRIMARY KEY, display_name TEXT)')
    for run in rows:
        con.execute('INSERT INTO stored_data VALUES (?, ?, ?)', (run, run.title(), f'files/{run}/scan_table.pq'))
        con.execute('INSERT INTO stored_files VALUES (?, ?)', (run, run.title()))
        (path / 'files' / run).mkdir()
        (path / 'files' / run / 'scan_table.pq').write_text(run)
    con.commit()
    con.close()


def test_example_results_added_once(tmp_path):
    import sqlite3
    example = tmp_path / 'example'
    _cache(example, ['example_a', 'example_b'])
    cache = tmp_path / 'workspace' / 'cache'
    cache.mkdir(parents=True)
    con = sqlite3.connect(cache / 'cache.db')
    con.execute('CREATE TABLE stored_data (id TEXT PRIMARY KEY, display_name TEXT)')
    con.execute("INSERT INTO stored_data VALUES ('example_a', 'mine')")
    con.commit()
    con.close()

    assert former_examples.add_example_results(cache, example) == ['example_b']
    con = sqlite3.connect(cache / 'cache.db')
    assert dict(con.execute('SELECT id, display_name FROM stored_data')) == {
        'example_a': 'mine', 'example_b': 'Example_B'}
    assert con.execute("SELECT scan_table FROM stored_data WHERE id='example_b'").fetchone() == (
        'files/example_b/scan_table.pq',)
    assert [r[0] for r in con.execute('SELECT id FROM stored_files')] == ['example_b']
    con.close()
    assert (cache / 'files' / 'example_b' / 'scan_table.pq').read_text() == 'example_b'
    assert not (cache / 'files' / 'example_a').exists()

    # Deleted by the user, it stays deleted
    import shutil
    shutil.rmtree(cache / 'files' / 'example_b')
    assert former_examples.add_example_results(cache, example) == []


def test_example_results_shipped():
    """The demo workspace's cache holds the example runs the Viewer lists."""
    import sqlite3
    con = sqlite3.connect(former_examples.EXAMPLE_CACHE / 'cache.db')
    runs = {r[0] for r in con.execute('SELECT id FROM stored_data')}
    con.close()
    assert {'example_deconvolution', 'example_native_identified'} <= runs
    for run in runs:
        assert (former_examples.EXAMPLE_CACHE / 'files' / run).is_dir()
