"""The JSON tab of the FLASHIda Method dialog: a JSON tree (streamlit-extras'
json_editor) whose changes are kept at once.

The tree returns the method as it was when it was drawn, and the other tabs
may add keys to the method afterwards. `apply_tree_edit` therefore applies only
what changed in the tree, onto the method as it is now.
"""

import pytest

from src.render.ida_runs import apply_tree_edit

DRAWN = {"global": {"method_name": "Native", "duration": 60}, "deconvolution": {"min_charge": 5}}
# The method now: the MS1 tab added a key after the tree was drawn
METHOD = {**DRAWN, "ms_settings": {"ms1": {"analyzer": ""}}}


def test_a_changed_value_is_set_and_later_keys_stay():
    after = {"global": {"method_name": "Native", "duration": 90}, "deconvolution": {"min_charge": 5}}
    new = apply_tree_edit(METHOD, DRAWN, after)
    assert new["global"]["duration"] == 90
    assert new["ms_settings"] == {"ms1": {"analyzer": ""}}


def test_added_and_removed_keys():
    after = {"global": {"method_name": "Native", "note": "new"}, "deconvolution": {"min_charge": 5}}
    new = apply_tree_edit(METHOD, DRAWN, after)
    assert new["global"] == {"method_name": "Native", "note": "new"}
    assert "ms_settings" in new


def test_a_removed_section():
    after = {"global": DRAWN["global"]}
    new = apply_tree_edit(METHOD, DRAWN, after)
    assert "deconvolution" not in new and "ms_settings" in new


def test_the_method_is_not_changed_in_place():
    after = {"global": {"method_name": "X", "duration": 60}, "deconvolution": {"min_charge": 5}}
    apply_tree_edit(METHOD, DRAWN, after)
    assert METHOD["global"]["method_name"] == "Native"


@pytest.mark.parametrize("after", [[1, 2], "text", 3, None])
def test_a_method_must_stay_an_object(after):
    with pytest.raises(ValueError, match="JSON object"):
        apply_tree_edit(METHOD, DRAWN, after)


def test_the_editor_is_streamlit_extras_json_editor():
    from streamlit_extras.json_editor import json_editor

    import src.render.ida_runs as ida_runs

    assert ida_runs.json_editor is json_editor
