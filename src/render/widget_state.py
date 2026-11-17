"""Point a widget at a value without overwriting what the user just picked."""

import streamlit as st


def follow(key, value):
    """Make the widget `key` show `value` on this run.

    Only assigned when the widget shows something else. A value assigned in
    session state is sent to the browser with the widget and replaces what the
    browser shows. Assigned on every run, a run that started before the user
    picked another value (one of the reruns the Insight components ask for after
    a switch) would put the old value back in the browser, and the next request
    from the page would carry it to the server, so the switch was undone.
    """
    if st.session_state.get(key) != value:
        st.session_state[key] = value
