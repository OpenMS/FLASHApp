"""The Match fragments page opens in a workspace without masses."""

import json

from streamlit.testing.v1 import AppTest


def test_fragment_matcher_without_masses(tmp_path):
    # With an openms-insight that can compare proteoforms, the page used to
    # read the stored masses before any were entered (stored is None).
    at = AppTest.from_file("content/Tools/FragmentMatcher.py")
    with open("settings.json", "r") as f:
        at.session_state.settings = json.load(f)
    at.session_state.settings["test"] = True
    at.session_state.settings["workspaces_dir"] = str(tmp_path)
    at.secrets["workspace"] = "test-fragment-matcher-empty"
    at.run(timeout=60)
    assert not at.exception
