"""OpenMS-Insight image downloads follow the sidebar's "image export format".

page_setup() hands the selectbox value to openms_insight.set_image_format(),
which keeps it in the session's state for every Insight component the page
draws. The call is skipped without error when openms-insight is missing or is
a release without set_image_format; those cases replace the package with a
stand-in in ``sys.modules`` that monkeypatch restores after each test.
"""

import json
import sys
import types

import pytest
from streamlit.testing.v1 import AppTest

PAGE = "content/user_guide.py"


@pytest.fixture
def app(tmp_path):
    at = AppTest.from_file(PAGE)
    with open("settings.json", "r", encoding="utf-8") as f:
        at.session_state.settings = json.load(f)
    at.session_state.settings["test"] = True
    # Keep the workspace this test creates out of the real workspaces folder.
    at.session_state.settings["workspaces_dir"] = str(tmp_path)
    at.secrets["workspace"] = "test"
    return at


def test_sidebar_format_reaches_insight(app):
    image_format = pytest.importorskip("openms_insight.core.image_format")
    app.run(timeout=60)
    assert not app.exception
    assert app.session_state[image_format._SESSION_KEY] == "svg"

    app.selectbox(key="image-format").set_value("png").run(timeout=60)
    assert not app.exception
    assert app.session_state[image_format._SESSION_KEY] == "png"


def test_insight_release_without_setting(app, monkeypatch):
    # An older release: the package imports, the function does not exist
    monkeypatch.setitem(sys.modules, "openms_insight", types.ModuleType("openms_insight"))
    app.run(timeout=60)
    assert not app.exception


def test_insight_not_installed(app, monkeypatch):
    # None in sys.modules makes any import of the name raise ImportError
    monkeypatch.setitem(sys.modules, "openms_insight", None)
    app.run(timeout=60)
    assert not app.exception


def test_options_follow_insight_list(app, monkeypatch):
    # The sidebar offers exactly the formats Insight lists, in its order
    fake = types.ModuleType("openms_insight")
    fake.IMAGE_FORMATS = ("png", "svg")
    fake.set_image_format = lambda fmt: None
    monkeypatch.setitem(sys.modules, "openms_insight", fake)
    app.run(timeout=60)
    assert not app.exception
    box = app.selectbox(key="image-format")
    assert box.options == ["png", "svg"]
    assert box.value == "svg"


def test_saved_format_insight_dropped(app, monkeypatch):
    # A saved format Insight no longer lists falls back to its first one
    fake = types.ModuleType("openms_insight")
    fake.IMAGE_FORMATS = ("png",)
    fake.set_image_format = lambda fmt: None
    monkeypatch.setitem(sys.modules, "openms_insight", fake)
    app.run(timeout=60)
    assert not app.exception
    assert app.selectbox(key="image-format").value == "png"
