"""The AI assistant reads FLASHApp's results from the workspace cache."""

import json
import shutil
from pathlib import Path

import pytest

from src.agent_access.config import agent_settings
from src.agent_access.tools import ToolError, Tools
from src.assistant_results import CacheRuns

EXAMPLE = Path("example-data/workspaces/default")


@pytest.fixture
def workspace(tmp_path):
    ws = tmp_path / "ws"
    shutil.copytree(EXAMPLE, ws)
    return ws


def test_cache_datasets_are_runs(workspace):
    runs = CacheRuns(workspace / "flashtnt")
    listed = {r["id"]: r for r in runs.list_runs()}
    assert listed["example_lysate_identified"]["kind"] == "identified"
    assert listed["example_deconvolution"]["kind"] == "deconvolved"
    assert listed["example_label_free"]["kind"] == "comparison"
    assert runs.get(runs.latest())["kind"] != "comparison"
    assert (runs.results_dir("example_lysate_identified") / "protein.tsv").is_file()


def test_no_cache_no_runs(tmp_path):
    assert CacheRuns(tmp_path / "flashtnt").list_runs() == []


def test_assistant_reads_and_shows_a_cached_run(workspace):
    tools = Tools(Path.cwd(), workspace, agent_settings(json.load(open("settings.json"))))
    ids = [r["id"] for r in tools.call("list_runs", {"workflow": "flashtnt"})["runs"]]
    assert "example_native_identified" in ids
    table = tools.call("read_result", {"workflow": "flashtnt", "run_id": "example_lysate_identified",
                                       "path": "protein.tsv", "limit": 2})
    assert "ProteinAccession" in table["columns"]
    view = tools.call("save_view", {"workflow": "flashtnt", "label": "Proteins",
                                    "components": [["protein_table"]],
                                    "runs": ["example_lysate_identified"],
                                    "question": "Which proteins?", "explanation": ["Listed."]})
    saved = json.loads((workspace / "layouts" / "viewer.json").read_text())
    assert saved["current"] == view["view"]
    with pytest.raises(ToolError):
        tools.call("save_view", {"workflow": "flashtnt", "label": "x", "components": [["protein_table"]],
                                 "runs": ["no_such_run"], "question": "q", "explanation": ["a"]})
