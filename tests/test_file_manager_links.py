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
