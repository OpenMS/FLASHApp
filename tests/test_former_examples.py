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
