"""FLASHIda: shown only in the local app in developer mode; raw file names,
the log folder of a method and the run logs the live viewer reads."""
import hashlib
import json
import os
from types import SimpleNamespace

import pytest

from src import flashida
from src.common import developer
from src.parse import ida_log


class _State(dict):
    __getattr__ = dict.__getitem__


@pytest.fixture
def state(monkeypatch, tmp_path):
    s = _State(settings={
        "online_deployment": False,
        "developer_mode": {"password_sha256": hashlib.sha256(b"open sesame").hexdigest()},
    })
    monkeypatch.setattr(developer, "st", SimpleNamespace(session_state=s))
    monkeypatch.setattr("streamlit.session_state", s)
    monkeypatch.setenv("LOCAL_WORKSPACES_DIR", str(tmp_path))
    return s


def test_only_in_the_local_app_with_developer_mode(state):
    assert not flashida.available()
    state["developer-mode-password"] = "open sesame"
    developer._unlock()
    assert flashida.available()
    # stays on after a restart
    for key in [k for k in state if k != "settings"]:
        del state[key]
    assert flashida.available()


def test_never_online(state):
    state["settings"]["online_deployment"] = True
    state["developer-mode-password"] = "open sesame"
    developer._unlock()
    assert developer.is_developer()
    assert not flashida.available()


@pytest.mark.parametrize("name, expected", [
    ("sample1_FLASHIda_TopN3.raw", "TopN3"),
    ("sample1_FLASHIda_TopN3_2.raw", "TopN3"),
    ("sample1_FLASHIda_TopN3.RAW", "TopN3"),
    ("sample1.raw", None),
    ("sample1_FLASHIda_TopN3.mzML", None),
])
def test_raw_file_names(name, expected):
    assert flashida.parse_raw_name(name) == expected


def test_find_raws(tmp_path):
    (tmp_path / "a_FLASHIda_M1.raw").write_text("")
    (tmp_path / "b_FLASHIda_M2_rep.raw").write_text("")
    (tmp_path / "other.raw").write_text("")
    (tmp_path / "x_FLASHIda_M3.raw").mkdir()
    assert flashida.find_raws(tmp_path) == {"a_FLASHIda_M1.raw": "M1", "b_FLASHIda_M2_rep.raw": "M2"}


def test_log_folder_of_a_method(tmp_path):
    exe_dir = tmp_path / "bin"
    exe_dir.mkdir()
    method = tmp_path / "TopN3.json"
    method.write_text(json.dumps({"runtime": {"log_dir": "logs"}}))
    assert flashida.method_log_dir(method, exe_dir) == (exe_dir / "logs").resolve()
    method.write_text(json.dumps({"runtime": {"log_dir": str(tmp_path / "abs")}}))
    assert flashida.method_log_dir(method, exe_dir) == (tmp_path / "abs").resolve()
    method.write_text("{}")
    assert flashida.method_log_dir(method, exe_dir) == exe_dir.resolve()

    flashida.remember_log_folder(tmp_path, "C:/a")
    flashida.remember_log_folder(tmp_path, "C:/b")
    flashida.remember_log_folder(tmp_path, "C:/a")
    assert flashida.log_folders(tmp_path) == ["C:/a", "C:/b"]


COMMANDS = [
    "tracking_id\tscan_type\tms_level\tparent_tracking_id\tprecursor_id\tpriority\tmono_mass\tcharge"
    "\tprecursor_mz\tisolation_width\tqscore\tactivation\tcollision_energy\tion_type\tion_index",
    "T0\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\tnone\t0\t\t0",
    "T1\texploration\t2\tT0\t1\t2\t12351.39\t15\t824.967\t1.86\t0.92\tHCD\t20\t\t0",
    "T2\texploration\t2\tT0\t1\t2\t12351.39\t15\t824.967\t1.86\t0.92\tHCD\t30\t\t0",
    "T3\trecording\t3\tT2\t1\t1\t12351.39;1498.85\t15;2\t824.967;750.934\t1.86;2\t0.92;0.86"
    "\tHCD;CID\t30;35\tb\t14",
    "T4\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\tnone\t0\t\t0",
    "T5\trecording\t2\tT4\t\"7\"\t2\t5000.5\t8\t626.07\t1.2\t0.5\tHCD\t30\t\t0",
]
RESULTS = [
    "tracking_id\tms_level\trt\tmass_count\ttag_count\texploration_score\tvariant_index"
    "\ttotal_variants\twinner_tracking_id\tdeconv_masses\tdeconv_intensities\tgroup_verdict",
    "T0\t1\t0.5\t2\t-1\t-1\t-1\t0\t\t12351.39;5000.5\t10,20;5\t",
    "T1\t2\t0.6\t100\t40\t100\t0\t2\t\t\t\t",
    "T2\t2\t0.6\t120\t50\t120\t1\t2\tT2\t\t\t",
    "T4\t1\t0.9\t0\t-1\t-1\t-1\t0\t\t\t\t",
]
IDENTIFICATION = [
    "tracking_id\tms_level\tprecursor_id\tproteoform\tflash_extender_score",
    "T1\t2\t1\tPEPTIDE\t200",
    "T2\t2\t1\tPEPTIDEK\t260",
]
POOLED = [
    "trigger_scan_id\ttrigger\tprecursor_id\tupdate_index\tproteoform\tflash_extender_score\tcoverage_pct",
    "T1\tMS2\t1\t1\tPEPTIDE\t200\t0.3",
    "T2\tMS2\t1\t2\tPEPTIDEK\t260\t0.45",
]


def _run(folder, partial=""):
    folder.mkdir(parents=True)
    for name, lines in (("scan_commands.tsv", COMMANDS), ("scan_results.tsv", RESULTS),
                        ("identification.tsv", IDENTIFICATION),
                        ("pooled_identification.tsv", POOLED)):
        (folder / name).write_text("\n".join(lines) + "\n" + (partial if name == "scan_results.tsv" else ""))
    return folder


def test_find_runs_newest_first(tmp_path):
    old = _run(tmp_path / "a_FLASHIda_M1_2026-11-03-10-00-00")
    new = _run(tmp_path / "b_FLASHIda_M1_2026-11-03-11-00-00")
    (tmp_path / "not-a-run").mkdir()
    os.utime(old / "scan_commands.tsv", (1, 1))
    assert ida_log.find_runs(tmp_path) == [new, old]
    assert ida_log.find_runs(new) == [new]
    assert ida_log.find_runs(tmp_path / "missing") == []


def test_a_line_being_written_is_not_read(tmp_path):
    run = _run(tmp_path / "run", partial="T5\t2\t1.0\t3")
    results = ida_log.read_table(run / "scan_results.tsv")
    assert results["tracking_id"].tolist() == ["T0", "T1", "T2", "T4"]
    before = ida_log.stamp(run)
    with open(run / "scan_results.tsv", "a") as f:
        f.write("\t-1\t-1\t-1\t0\t\t\t\t\n")
    assert ida_log.read_table(run / "scan_results.tsv")["tracking_id"].tolist()[-1] == "T5"
    assert ida_log.stamp(run) != before


def test_decisions(tmp_path):
    result = ida_log.decisions(ida_log.read_run(_run(tmp_path / "run")))

    assert result["summary"] == {
        "surveys": 2, "surveys_acquired": 2, "ms2": 3, "ms3": 1, "waiting": 2,
        "precursors": 2, "identified": 1, "rt": 0.9,
    }

    scans = result["scans"].set_index("Scan")
    assert scans.loc["T1", "precursor"] == 1 and scans.loc["T5", "precursor"] == 2
    assert scans.loc["T2", "Winner"] and not scans.loc["T1", "Winner"]
    assert scans.loc["T1", "Variant"] == "1/2"
    assert scans.loc["T3", "Fragment"] == "b14"
    assert scans.loc["T3", "Mass"] == pytest.approx(1498.85)
    assert scans.loc["T3", "Activation"] == "CID"
    assert scans.loc["T3", "Status"] == "waiting"
    assert scans.loc["T2", "Proteoform"] == "PEPTIDEK"

    precursors = result["precursors"].set_index("precursor")
    assert precursors.loc[1, "Mass"] == pytest.approx(12351.39)
    assert precursors.loc[1, "RT"] == pytest.approx(0.5)
    assert (precursors.loc[1, "MS2"], precursors.loc[1, "MS3"]) == (2, 1)
    assert precursors.loc[1, "Best energy"] == 30
    assert precursors.loc[1, "Proteoform"] == "PEPTIDEK"
    assert precursors.loc[1, "Coverage"] == pytest.approx(45)
    assert precursors.loc[2, "RT"] == pytest.approx(0.9)
    assert not precursors.loc[2, "Identified"]

    points = result["points"]
    survey = points[points["Decision"] == "Deconvolved"]
    assert survey["mass"].tolist() == [12351.39, 5000.5]
    assert survey["intensity"].tolist() == [30, 5]
    assert set(points["Decision"]) == {"Deconvolved", "Identified", "Chosen"}


def test_an_empty_run(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    (run / "scan_commands.tsv").write_text(COMMANDS[0] + "\n")
    result = ida_log.decisions(ida_log.read_run(run))
    assert result["summary"]["precursors"] == 0
    assert result["points"].empty
