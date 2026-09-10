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
from src.render.util import payload_hash


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


def test_all_none_first_message_does_not_claim_keys():
    """The frontend sends null for every unset field, so the first cell to report
    (e.g. a spectrum) carries scanIndex=None, massIndex=None, ... Those must not be
    adopted as state: otherwise the Scan Table's default row selection, sent with the
    same counter, is treated as a stale conflict and the dependent cells stay empty
    until the user clicks a row (OpenMS/FLASHApp#100)."""
    tracker = StateTracker()
    all_none = _echo_with(tracker, scanIndex=None, massIndex=None, AApos=None)
    assert tracker.updateState(all_none) is False
    assert tracker.currentStateCounter == 0
    assert "scanIndex" not in tracker.getState()

    default_row = _echo_with(tracker, scanIndex=0, massIndex=0)  # still counter 0
    assert tracker.updateState(default_row) is True
    assert tracker.getState()["scanIndex"] == 0
    assert tracker.getState()["massIndex"] == 0


def test_none_for_a_never_set_key_is_a_no_op():
    """Clearing an existing key keeps its semantics (see above); a None for a key that
    was never set changes nothing and must not raise."""
    tracker = StateTracker()
    tracker.updateState(_echo_with(tracker, AApos=5))
    assert tracker.updateState(_echo_with(tracker, AApos=None, tagIndex=None)) is True
    echoed = tracker.getState()
    assert echoed["AApos"] is None
    assert "tagIndex" not in echoed


def test_payload_hash_changes_with_the_tracker():
    """Two runs of the same input give byte-identical tables. The payload hash must
    still differ per experiment (= per StateTracker), or the frontend skips the
    render, never rebuilds its tables and never seeds the default row selection."""
    data = {"per_scan_data": [{"index": 0, "Scan": 3098}]}
    first, second = StateTracker(), StateTracker()
    assert payload_hash(data, first.id) != payload_hash(data, second.id)
    assert payload_hash(data, first.id) == payload_hash(dict(data), first.id)
