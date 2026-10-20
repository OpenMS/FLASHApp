"""
Tests for FileManager copy-on-write protection of demo "ground truth".

Bug being fixed: when a demo workspace is loaded in online mode on Linux, the
workspace is materialized by symlinking every committed demo file back to the
read-only ground truth under ``example-data/`` (see ``_symlink_tree`` /
``copy_demo_workspace`` in ``src/common/common.py``). The committed demo ships
writable cache artifacts (``cache.db`` and ``cache/files/<dataset>/*``). Because
writing to a symlink follows it to its target, reprocessing data used to write
through those symlinks and overwrite the committed ground truth:

  * ``cache.db`` (guaranteed): ``sqlite3`` writes the database in place on every
    ``store_*`` call.
  * result files: clobbered whenever the reprocessed ``dataset_id`` collided
    with a demo dataset id.

FileManager now performs copy-on-write: it materializes a real file in the
workspace before any in-place write, so the ground truth is never modified.

These tests pin that behaviour by simulating a symlinked cache and asserting the
ground-truth bytes survive a reprocess while the workspace diverges.
"""

import os
import sys

import pytest

# FileManager imports these at module import time; skip the suite if a stripped
# environment lacks them (matches the convention of the other test modules).
pd = pytest.importorskip("pandas")
pl = pytest.importorskip("polars")
pytest.importorskip("pyarrow")

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pathlib import Path

from src.workflow.FileManager import FileManager


# Result artifacts the demo ships, one per FileManager write path:
PICKLE_TAG = "orig_pickle"   # store_data -> .pkl.gz   (gzip.open write branch)
FRAME_TAG = "orig_frame"     # store_data -> .pq        (pandas to_parquet branch)
POLARS_TAG = "orig_polars"   # store_data -> .pq        (polars write_parquet branch)
FILE_TAG = "orig_file"       # store_file -> .bin       (shutil.copy branch)
DEMO_DS = "demods"


def _seed_ground_truth(gt_cache: Path) -> None:
    """Populate a cache dir the way the committed demo does: a real cache.db
    plus real result files under files/<DEMO_DS>/."""
    fm = FileManager(gt_cache.parent / "gt-wf", cache_path=gt_cache)
    fm.store_data(DEMO_DS, PICKLE_TAG, {"demo": True})
    fm.store_data(DEMO_DS, FRAME_TAG, pd.DataFrame({"a": [1, 2, 3]}))
    fm.store_data(DEMO_DS, POLARS_TAG, pl.DataFrame({"b": [4, 5, 6]}))
    src = gt_cache.parent / "seed_src.bin"
    src.write_bytes(b"GROUND TRUTH BIN")
    fm.store_file(DEMO_DS, FILE_TAG, src, file_name="orig_file.bin")
    fm.cache_connection.close()


def _symlink_workspace(gt_cache: Path, ws_cache: Path) -> None:
    """Mimic _symlink_tree: real directories, but the cache files are symlinks
    pointing at the ground-truth source (absolute, like item.resolve())."""
    (ws_cache / "files" / DEMO_DS).mkdir(parents=True)
    (ws_cache / "cache.db").symlink_to((gt_cache / "cache.db").resolve())
    for name in (
        f"{PICKLE_TAG}.pkl.gz",
        f"{FRAME_TAG}.pq",
        f"{POLARS_TAG}.pq",
        "orig_file.bin",
    ):
        link = ws_cache / "files" / DEMO_DS / name
        link.symlink_to((gt_cache / "files" / DEMO_DS / name).resolve())


def _gt_files(gt_cache: Path) -> list[Path]:
    base = gt_cache / "files" / DEMO_DS
    return [
        gt_cache / "cache.db",
        base / f"{PICKLE_TAG}.pkl.gz",
        base / f"{FRAME_TAG}.pq",
        base / f"{POLARS_TAG}.pq",
        base / "orig_file.bin",
    ]


def test_cache_db_symlink_is_materialized_on_connect(tmp_path):
    """Opening a workspace whose cache.db is a symlink must replace it with an
    independent real copy (preserving the demo index) and leave ground truth
    untouched."""
    gt_cache = tmp_path / "ground_truth" / "cache"
    gt_cache.mkdir(parents=True)
    _seed_ground_truth(gt_cache)

    gt_db = gt_cache / "cache.db"
    gt_db_bytes = gt_db.read_bytes()

    ws_cache = tmp_path / "workspace" / "cache"
    _symlink_workspace(gt_cache, ws_cache)
    assert (ws_cache / "cache.db").is_symlink()  # precondition

    fm = FileManager(tmp_path / "workspace" / "wf", cache_path=ws_cache)
    try:
        # Workspace cache.db is now a real, independent file...
        assert not (ws_cache / "cache.db").is_symlink()
        assert (ws_cache / "cache.db").is_file()
        # ...the demo's index rows were preserved in the workspace copy...
        assert fm.result_exists(DEMO_DS, PICKLE_TAG)
        # ...and the ground-truth database is byte-for-byte unchanged.
        assert gt_db.read_bytes() == gt_db_bytes
    finally:
        fm.cache_connection.close()


def test_reprocess_does_not_write_through_symlinks(tmp_path):
    """Reprocessing a dataset whose id collides with the demo must not modify
    any ground-truth file; the workspace copies diverge instead."""
    gt_cache = tmp_path / "ground_truth" / "cache"
    gt_cache.mkdir(parents=True)
    _seed_ground_truth(gt_cache)

    snapshot = {p: p.read_bytes() for p in _gt_files(gt_cache)}

    ws_cache = tmp_path / "workspace" / "cache"
    _symlink_workspace(gt_cache, ws_cache)

    fm = FileManager(tmp_path / "workspace" / "wf", cache_path=ws_cache)
    try:
        # Reprocess: overwrite every artifact through what are currently symlinks.
        fm.store_data(DEMO_DS, PICKLE_TAG, {"reprocessed": 999})
        fm.store_data(DEMO_DS, FRAME_TAG, pd.DataFrame({"a": [10, 20]}))
        fm.store_data(DEMO_DS, POLARS_TAG, pl.DataFrame({"b": [70, 80]}))
        new_src = tmp_path / "new_src.bin"
        new_src.write_bytes(b"WORKSPACE NEW BIN")
        fm.store_file(DEMO_DS, FILE_TAG, new_src, file_name="orig_file.bin")

        # Ground truth is completely unchanged.
        for p, original in snapshot.items():
            assert p.read_bytes() == original, f"ground truth modified: {p}"
        assert (gt_cache / "files" / DEMO_DS / "orig_file.bin").read_bytes() == (
            b"GROUND TRUTH BIN"
        )

        # Workspace files are now real files (not symlinks) with new content.
        ws_files = ws_cache / "files" / DEMO_DS
        for name in (
            f"{PICKLE_TAG}.pkl.gz",
            f"{FRAME_TAG}.pq",
            f"{POLARS_TAG}.pq",
            "orig_file.bin",
        ):
            wp = ws_files / name
            assert not wp.is_symlink(), f"{name} is still a symlink into ground truth"
            assert wp.is_file()
        assert (ws_files / "orig_file.bin").read_bytes() == b"WORKSPACE NEW BIN"

        # Reading back through the workspace FileManager returns new content.
        assert fm.get_results(DEMO_DS, [PICKLE_TAG])[PICKLE_TAG] == {"reprocessed": 999}
        assert fm.get_results(DEMO_DS, [FRAME_TAG])[FRAME_TAG]["a"].tolist() == [10, 20]
    finally:
        fm.cache_connection.close()
