"""
Tests for the pure IDA-log -> mzML mapping helpers.

These live in src/workflow/_ida_log.py so the log discovery and filename
auto-match logic is unit-testable without booting Streamlit or pulling in
pyopenms (mirrors tests/test_log_status.py).
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.workflow._ida_log import (
    IDA_LOG_KEY,
    IDA_NONE,
    available_ida_logs,
    auto_match_log,
)


def _make_log_dir(tmp_path):
    d = tmp_path / "input-files" / IDA_LOG_KEY
    d.mkdir(parents=True)
    return d


def test_available_logs_empty_when_no_dir(tmp_path):
    assert available_ida_logs(tmp_path) == []


def test_available_logs_lists_only_log_files(tmp_path):
    d = _make_log_dir(tmp_path)
    (d / "a.log").write_text("x")
    (d / "b.log").write_text("y")
    (d / "notes.txt").write_text("z")
    found = available_ida_logs(tmp_path)
    assert sorted(os.path.basename(p) for p in found) == ["a.log", "b.log"]


def test_available_logs_skips_external_files_marker(tmp_path):
    d = _make_log_dir(tmp_path)
    (d / "a.log").write_text("x")
    (d / "external_files.txt").write_text("")
    found = available_ida_logs(tmp_path)
    assert [os.path.basename(p) for p in found] == ["a.log"]


def test_available_logs_includes_existing_external_paths_only(tmp_path):
    d = _make_log_dir(tmp_path)
    external = tmp_path / "ext.log"
    external.write_text("x")
    missing = tmp_path / "gone.log"  # referenced but never created
    (d / "external_files.txt").write_text(f"{external}\n{missing}\n")
    found = available_ida_logs(tmp_path)
    assert str(external) in found
    assert str(missing) not in found


def test_auto_match_by_stem():
    logs = ["/data/sample1.log", "/data/other.log"]
    assert auto_match_log("sample1.mzML", logs) == "/data/sample1.log"


def test_auto_match_returns_none_sentinel_when_no_match():
    assert auto_match_log("sample1.mzML", ["/data/other.log"]) == IDA_NONE


def test_auto_match_handles_full_path_mzml():
    assert auto_match_log("/inputs/run3.mzML", ["/logs/run3.log"]) == "/logs/run3.log"


def test_auto_match_empty_logs_returns_none_sentinel():
    assert auto_match_log("anything.mzML", []) == IDA_NONE
