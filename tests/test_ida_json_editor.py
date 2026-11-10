"""The JSON tab of the FLASHIda Method dialog: a JSON tree (streamlit-extras'
json_editor) whose changes are kept at once.

`edited_method` decides what the tree's value means for the run's method.
"""

import pytest

from src.render.ida_runs import edited_method

METHOD = {"global": {"method_name": "Native", "duration": 60}, "deconvolution": {"min_charge": 5}}


def test_unchanged_tree_keeps_the_method():
    assert edited_method({"global": {"method_name": "Native", "duration": 60},
                          "deconvolution": {"min_charge": 5}}, METHOD) is None


def test_a_changed_value_is_the_new_method():
    edited = {"global": {"method_name": "Native", "duration": 90}, "deconvolution": {"min_charge": 5}}
    new = edited_method(edited, METHOD)
    assert new == edited
    # A copy: the run keeps it even when the tree's value changes again
    edited["global"]["duration"] = 120
    assert new["global"]["duration"] == 90


def test_added_and_removed_keys_are_kept():
    edited = {"global": {"method_name": "Native", "duration": 60, "note": "new key"}}
    assert edited_method(edited, METHOD) == edited


@pytest.mark.parametrize("edited", [[1, 2], "text", 3, None])
def test_a_method_must_stay_an_object(edited):
    with pytest.raises(ValueError, match="JSON object"):
        edited_method(edited, METHOD)


def test_the_editor_is_streamlit_extras_json_editor():
    from streamlit_extras.json_editor import json_editor

    import src.render.ida_runs as ida_runs

    assert ida_runs.json_editor is json_editor
