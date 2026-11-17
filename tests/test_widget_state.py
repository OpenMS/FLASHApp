"""A widget pointed at a value keeps what the browser shows unless the value changes.

A value assigned in session state goes to the browser with the widget (its
``set_value``) and replaces the user's pick there. The Viewer's experiment
selector was assigned on every run, so a run that had started before the user
picked another experiment put the old one back, and the switch was undone.
"""

from streamlit.testing.v1 import AppTest


def _page():
    import streamlit as st
    from src.render.widget_state import follow

    follow("picked", st.session_state.get("shown", "a"))
    st.selectbox("experiment", ["a", "b", "c"], key="picked")


def test_unchanged_value_is_not_sent_again():
    at = AppTest.from_function(_page).run()
    assert at.selectbox(key="picked").value == "a"
    at.run()
    assert not at.selectbox(key="picked").proto.set_value


def test_changed_value_is_sent():
    at = AppTest.from_function(_page).run()
    at.session_state["shown"] = "c"
    at.run()
    assert at.selectbox(key="picked").value == "c"
    assert at.selectbox(key="picked").proto.set_value


def test_pick_by_user_stays():
    at = AppTest.from_function(_page).run()
    at.session_state["shown"] = "b"
    at.selectbox(key="picked").set_value("b").run()
    assert at.selectbox(key="picked").value == "b"
    assert not at.selectbox(key="picked").proto.set_value
