"""
A workspace started from a demo links the demo's large result files
(src/common/common.py, copy_demo_workspace). Storing a result again must
replace the link with a file of the workspace's own, never write through it
into the shared demo.
"""
from io import BytesIO
from pathlib import Path

import polars as pl

from src.workflow.FileManager import FileManager


def _linked(tmp_path, name):
    demo = tmp_path / "demo" / name
    demo.parent.mkdir(parents=True, exist_ok=True)
    demo.write_bytes(b"demo")
    fm = FileManager(tmp_path / "wf")
    target = Path(fm.cache_path, "files", "run", name)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.symlink_to(demo)
    return fm, demo, target


def test_store_data_replaces_a_link(tmp_path):
    fm, demo, target = _linked(tmp_path, "masses.pq")
    fm.store_data("run", "masses", pl.DataFrame({"a": [1]}))
    assert demo.read_bytes() == b"demo"
    assert not target.is_symlink()
    assert fm.get_results("run", ["masses"], use_polars=True)["masses"].collect()["a"][0] == 1


def test_store_file_replaces_a_link(tmp_path):
    fm, demo, target = _linked(tmp_path, "spectra.mzML")
    fm.store_file("run", "spectra", BytesIO(b"new"), file_name="spectra.mzML")
    assert demo.read_bytes() == b"demo"
    assert not target.is_symlink() and target.read_bytes() == b"new"


def _linked_folder(tmp_path):
    """A workspace whose result folder 'run' is linked whole from the demo
    (``link_folders`` in settings.json)."""
    demo_run = tmp_path / "demo" / "run"
    demo_run.mkdir(parents=True)
    (demo_run / "masses.pq").write_bytes(b"demo")
    fm = FileManager(tmp_path / "wf")
    folder = Path(fm.cache_path, "files", "run")
    folder.symlink_to(demo_run, target_is_directory=True)
    return fm, demo_run, folder


def test_storing_into_a_linked_folder_leaves_the_demo_alone(tmp_path):
    fm, demo_run, folder = _linked_folder(tmp_path)
    fm.store_data("run", "edits", pl.DataFrame({"a": [1]}))
    fm.store_file("run", "report", BytesIO(b"new"), file_name="report.md")
    assert sorted(p.name for p in demo_run.iterdir()) == ["masses.pq"]
    assert not folder.is_symlink()
    assert (folder / "masses.pq").read_bytes() == b"demo"
    assert fm.get_results("run", ["edits"], use_polars=True)["edits"].collect()["a"][0] == 1


def test_parquet_sink_into_a_linked_folder_leaves_the_demo_alone(tmp_path):
    fm, demo_run, folder = _linked_folder(tmp_path)
    with fm.parquet_sink("run", "sunk") as tmp:
        pl.DataFrame({"a": [2]}).write_parquet(tmp)
    assert not (demo_run / "sunk.pq").exists()
    assert (folder / "sunk.pq").is_file()


def test_removing_a_linked_result_keeps_the_demo(tmp_path):
    fm, demo_run, folder = _linked_folder(tmp_path)
    fm.remove_results("run")
    assert not folder.exists()
    assert (demo_run / "masses.pq").read_bytes() == b"demo"


def test_archive_of_a_linked_result_is_built_in_the_workspace(tmp_path):
    from src.common.result_download import archive_path, build_archive

    fm, demo_run, folder = _linked_folder(tmp_path)
    target = build_archive([folder / "masses.pq"], archive_path(fm, "run"))
    assert target.is_file()
    assert sorted(p.name for p in demo_run.iterdir()) == ["masses.pq"]
