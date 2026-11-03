"""FLASHIda: shown only in the local app in developer mode; raw file names,
the method file and the MS2 table of the viewer."""
import hashlib
import xml.etree.ElementTree as ET
from types import SimpleNamespace

import pandas as pd
import pytest

from src import flashida
from src.common import developer
from src.parse import ida


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
    ("sample1_FLASHIda_TopN3.raw", ("TopN3", False)),
    ("sample1_FLASHIda_TopN3_2.raw", ("TopN3", True)),
    ("sample1_FLASHIda_TopN3_1_rep.raw", ("TopN3", False)),
    ("sample1_FLASHIda_TopN3.RAW", ("TopN3", False)),
    ("sample1.raw", None),
    ("sample1_FLASHIda_TopN3.mzML", None),
])
def test_raw_file_names(name, expected):
    assert flashida.parse_raw_name(name) == expected


def test_find_raws(tmp_path):
    (tmp_path / "a_FLASHIda_M1.raw").write_text("")
    (tmp_path / "b_FLASHIda_M2_2.raw").write_text("")
    (tmp_path / "other.raw").write_text("")
    (tmp_path / "x_FLASHIda_M3.raw").mkdir()
    assert flashida.find_raws(tmp_path) == {
        "a_FLASHIda_M1.raw": ("M1", False), "b_FLASHIda_M2_2.raw": ("M2", True),
    }


def test_method_xml(tmp_path):
    path = tmp_path / "method.xml"
    flashida.write_method_xml({
        flashida.method_key("MS1", "TopN"): 5,
        flashida.method_key("IDA", "Tolerances"): "5.0\n7.5\n",
        flashida.method_key("IDA", "UseFAIMS"): True,
        flashida.method_key("IDA", "TargetLogs"): r"C:\logs\a.log",
    }, path)
    root = ET.parse(path).getroot()
    assert root.tag == "MethodParameters"
    assert [c.tag for c in root] == ["MS1", "MS2", "IDA"]
    assert root.find("MS1/TopN").text == "5"
    assert root.find("MS2/Activation").text == "HCD"
    assert [d.text for d in root.find("IDA/Tolerances")] == ["5.0", "7.5"]
    assert [d.tag for d in root.find("IDA/CVValues")] == ["double"] * 5
    assert root.find("IDA/UseFAIMS").text == "true"
    assert root.find("IDA/UseCVQScore").text == "false"
    assert [d.text for d in root.find("IDA/TargetLogs")] == [r"C:\logs\a.log"]


def test_ms2_table_links_scans_and_precursors():
    scan_table = pd.DataFrame({"index": [0, 1, 2, 3], "Scan": [10, 11, 12, 13]})
    spec2 = pd.DataFrame({
        "ScanNum": [11, 11, 12, 13], "RetentionTime": [1.0, 1.0, 2.0, 3.0],
        "PrecursorScanNum": [10, 10, 10, 99], "PrecursorMz": [500.0, 500.0, 600.0, 700.0],
        "MassCountInSpec": [2, 2, 1, 1],
    })
    table = ida.ms2_table(spec2, scan_table, identified={12})
    assert table["Scan"].tolist() == [11, 12, 13]
    assert table["index"].tolist() == [1, 2, 3]
    assert table["precursor_index"].tolist() == [0, 0, -1]
    assert table["Identified"].tolist() == [False, True, False]


def test_identified_scans_skip_the_toppic_header(tmp_path):
    path = tmp_path / "run_ms2_toppic_prsm_single.tsv"
    header = "********************** Parameters **********************\n" + "Key:\tvalue\n" * 7
    path.write_text(
        header + "Data file name\tPrSM ID\tScan(s)\tE-value\n"
        "run.msalign\t0\t3099\t1e-10\nrun.msalign\t1\t3105 3106\t1e-5\n"
    )
    assert ida.identified_scans(path) == {3099, 3105, 3106}


def test_simulation_heatmap():
    sim = pd.DataFrame({"monoMasses": [1000.0, 2000.0], "rt": [1.0, 2.0],
                        "precursorIntensity": [5.0, 6.0], "other": [1, 2]})
    heat = ida.simulation_heatmap(sim)
    assert list(heat.columns) == ["row", "mass", "rt", "intensity"]
    assert heat["row"].tolist() == [0, 1]
