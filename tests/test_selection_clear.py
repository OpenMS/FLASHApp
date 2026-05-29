"""
Tests for the selection-clearing round-trip used by the FLASHViewer grid.

Each view (Sequence View, Tag Table, Protein Table, ...) is a separate Streamlit
component instance with its own frontend store; they share selection state only by
round-tripping through Python's StateTracker. Clearing a selection (e.g. deselecting
an amino acid, or switching proteoform) must therefore propagate back to every view.

The frontend sends a cleared field as `null`/`None` (App.vue maps `undefined -> null`
so the clear survives JSON serialization). These tests pin the two invariants the fix
relies on:

  1. A cleared field is echoed back as `None` so every component can clear it.
  2. render_component strips `None`-valued keys for the data computation, preserving
     update.py's "key not in selection_store" convention.

They also document the original bug: when the cleared key was *dropped* from the
payload entirely, the merge-only StateTracker kept echoing the stale value.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.render.StateTracker import StateTracker


def _echo_with(tracker, **overrides):
    """Mimic a component returning the echoed state with `overrides` applied."""
    state = tracker.getState()  # includes counter + id, like getState() -> frontend
    state.update(overrides)
    return state


def _active_state(state):
    """The view render.py passes to update/filter: None == "not selected" == absent."""
    return {k: v for k, v in state.items() if v is not None}


def test_selecting_a_value_round_trips():
    tracker = StateTracker()
    tracker.updateState(_echo_with(tracker, AApos=5))
    assert tracker.getState()["AApos"] == 5
    assert _active_state(tracker.getState())["AApos"] == 5


def test_clearing_a_selection_round_trips_as_none():
    tracker = StateTracker()
    tracker.updateState(_echo_with(tracker, AApos=5))
    assert tracker.getState()["AApos"] == 5

    # Deselect: the frontend sends AApos=None (App.vue maps undefined -> null).
    tracker.updateState(_echo_with(tracker, AApos=None))
    echoed = tracker.getState()

    # (1) Echoed back as None so every component clears the field locally.
    assert echoed["AApos"] is None
    # (2) The data-computation view treats None as absent (not selected).
    assert "AApos" not in _active_state(echoed)


def test_dropped_key_keeps_stale_value_regression():
    """Pre-fix behavior: `undefined` was dropped from the payload, so the merge-only
    StateTracker never learned about the clear and kept echoing the stale value.
    This is exactly the bug the null-bridge (send None instead of dropping) fixes."""
    tracker = StateTracker()
    tracker.updateState(_echo_with(tracker, AApos=5))

    payload = tracker.getState()
    payload.pop("AApos")  # simulate the JSON-dropped undefined key
    tracker.updateState(payload)

    assert tracker.getState()["AApos"] == 5  # stale value survives -> the original bug
