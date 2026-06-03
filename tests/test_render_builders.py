"""Construct-smoke for ``src.render.render.make_builders`` + the frozen grid.

For each tool: build synthetic FileManager caches, run ``build_insight_caches``,
then ``make_builders``; call every builder to actually construct the OpenMS-Insight
component (which triggers subprocess preprocessing over ``data_path=`` and a disk
cache), and assert ``_prepare_vue_data`` / ``_get_component_args`` run over that
cached data. Then drive the frozen ``render_linked_grid`` with a patched render
bridge so the grid wiring (shared StateManager + per-cell keys) is exercised
without touching the Vue layer.

This is intentionally NOT a Streamlit ``AppTest`` (Insight's spawn-multiprocessing
preprocessing is incompatible with AppTest's runtime).
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from openms_insight import StateManager

from src.workflow.FileManager import FileManager
from src.render.render import make_builders
from src.render.schema import build_insight_caches
from src.view.grid import render_linked_grid
from tests.conftest import (
    make_deconv_caches,
    make_tnt_caches,
    make_quant_caches,
    make_sequence_cache,
)


def _fm(workspace):
    return FileManager(workspace, Path(workspace, "cache"))


# Layout per tool -> the comp_names the smoke must construct + render.
DECONV_COMPS = [
    "scan_table", "mass_table", "deconv_spectrum", "anno_spectrum",
    "combined_spectrum", "3D_SN_plot", "ms1_deconv_heat_map", "ms2_deconv_heat_map",
    "ms1_raw_heatmap", "ms2_raw_heatmap", "fdr_plot", "sequence_view",
]
TNT_COMPS = [
    "protein_table", "tag_table", "sequence_view", "combined_spectrum",
    "id_fdr_plot", "scan_table", "mass_table",
]
QUANT_COMPS = ["quant_visualization", "quant_traces_3d"]


def _exercise_builder(builder, sm):
    """Construct one component and run its two data-shaping hooks over its cache.

    Components are duck-typed: every Insight visualization is callable and exposes
    ``_prepare_vue_data`` / ``_get_component_args`` (``SequenceView`` is the one
    component that is not a ``BaseComponent`` subclass but honors the same surface).
    """
    comp = builder()
    assert callable(comp)
    assert hasattr(comp, "_prepare_vue_data") and hasattr(comp, "_get_component_args")
    state = sm.get_state_for_vue()
    vue_data = comp._prepare_vue_data(state)
    assert isinstance(vue_data, dict) and len(vue_data) > 0
    args = comp._get_component_args()
    assert "componentType" in args
    return comp


# --------------------------------------------------------------------------- #
# make_builders signature + per-component construction
# --------------------------------------------------------------------------- #
def test_make_builders_returns_zero_arg_factories(mock_streamlit, temp_workspace):
    fm = _fm(temp_workspace)
    ds = make_deconv_caches(fm)
    make_sequence_cache(fm)
    build_insight_caches(fm, ds, "flashdeconv")

    builders = make_builders(fm, ds, "flashdeconv")
    assert isinstance(builders, dict)
    # every value is a zero-arg callable factory
    for name, factory in builders.items():
        assert callable(factory), name


def test_builders_construct_and_prepare_flashdeconv(mock_streamlit, temp_workspace):
    fm = _fm(temp_workspace)
    ds = make_deconv_caches(fm)
    make_sequence_cache(fm)
    build_insight_caches(fm, ds, "flashdeconv")

    sm = StateManager(session_key=f"flashdeconv__{ds}")
    builders = make_builders(fm, ds, "flashdeconv")
    for name in DECONV_COMPS:
        assert name in builders, name
        comp = _exercise_builder(builders[name], sm)
        # cache_id carries the dataset -> per-dataset reset guarantee
        assert comp._cache_id == f"flashdeconv__{ds}__{name}"


def test_builders_construct_and_prepare_flashtnt(mock_streamlit, temp_workspace):
    fm = _fm(temp_workspace)
    ds = make_tnt_caches(fm)
    settings = fm.get_results(ds, ["settings"])["settings"]
    build_insight_caches(fm, ds, "flashtnt")

    sm = StateManager(session_key=f"flashtnt__{ds}")
    builders = make_builders(fm, ds, "flashtnt", settings=settings)
    for name in TNT_COMPS:
        assert name in builders, name
        _exercise_builder(builders[name], sm)


def test_builders_construct_and_prepare_flashquant(mock_streamlit, temp_workspace):
    fm = _fm(temp_workspace)
    ds = make_quant_caches(fm)
    build_insight_caches(fm, ds, "flashquant")

    sm = StateManager(session_key=f"flashquant__{ds}")
    builders = make_builders(fm, ds, "flashquant")
    for name in QUANT_COMPS:
        assert name in builders, name
        _exercise_builder(builders[name], sm)


# --------------------------------------------------------------------------- #
# value-based cross-link selection (index -> value migration)
# --------------------------------------------------------------------------- #
def test_filters_interactivity_value_based(mock_streamlit, temp_workspace):
    """scan/mass/protein selection is value-based via filters/interactivity."""
    fm = _fm(temp_workspace)
    ds = make_tnt_caches(fm)
    build_insight_caches(fm, ds, "flashtnt")
    builders = make_builders(fm, ds, "flashtnt")

    scan_table = builders["scan_table"]()
    assert scan_table.get_interactivity_mapping() == {"scan": "scan_id"}

    mass_table = builders["mass_table"]()
    assert mass_table.get_filters_mapping() == {"scan": "scan_id"}
    assert mass_table.get_interactivity_mapping() == {"mass": "mass_id"}

    plot3d = builders["3D_SN_plot"]()
    # massIndex -> value filter on mass_in_scan; scanIndex -> scan
    assert plot3d.get_filters_mapping() == {"scan": "scan_id", "mass": "mass_in_scan"}

    tag_table = builders["tag_table"]()
    # proteinIndex + proteoform_scan_map collapse to a precomputed protein_id filter
    assert tag_table.get_filters_mapping() == {"protein": "protein_id"}
    assert tag_table.get_interactivity_mapping() == {"tag": "tag_id"}


def test_scan_to_mass_filter_applies(mock_streamlit, temp_workspace):
    """Selecting a scan filters the mass table to that scan's masses (value-based)."""
    fm = _fm(temp_workspace)
    ds = make_deconv_caches(fm)
    build_insight_caches(fm, ds, "flashdeconv")
    builders = make_builders(fm, ds, "flashdeconv")

    mass_table = builders["mass_table"]()
    # scan_id 0 has 2 masses, scan_id 1 has 1 mass
    d0 = mass_table._prepare_vue_data({"scan": 0})["tableData"]
    d1 = mass_table._prepare_vue_data({"scan": 1})["tableData"]
    assert len(d0) == 2
    assert len(d1) == 1


# --------------------------------------------------------------------------- #
# the frozen grid renders the builders against a shared StateManager
# --------------------------------------------------------------------------- #
def test_render_linked_grid_exercises_components(mock_streamlit, temp_workspace):
    """render_linked_grid builds each cell's component + runs its data hooks.

    The Vue render bridge is patched out; the patch calls each component's
    ``_prepare_vue_data`` / ``_get_component_args`` so the grid's
    build->prepare->render path is exercised end-to-end without spawning the
    front-end. Asserts a single shared StateManager and per-cell keys.
    """
    fm = _fm(temp_workspace)
    ds = make_deconv_caches(fm)
    make_sequence_cache(fm)
    build_insight_caches(fm, ds, "flashdeconv")
    builders = make_builders(fm, ds, "flashdeconv")

    rendered = []  # (cache_id, key, state_manager_session_key)

    def fake_render(component, state_manager, key=None, height=None):
        state = state_manager.get_state_for_vue()
        component._prepare_vue_data(state)
        component._get_component_args()
        rendered.append((component._cache_id, key, state_manager._session_key))
        return None

    layout = [
        ["scan_table", "mass_table"],
        ["anno_spectrum", "deconv_spectrum"],
        ["3D_SN_plot"],
    ]
    with patch("openms_insight.rendering.bridge.render_component", fake_render):
        sm = render_linked_grid(layout, builders, state_key=f"flashdeconv__{ds}")

    assert isinstance(sm, StateManager)
    # every cell rendered (5 panels)
    assert len(rendered) == 5
    # all panels shared ONE StateManager session_key (cross-linking)
    assert {r[2] for r in rendered} == {f"flashdeconv__{ds}"}
    # per-cell keys follow the f"{grid_key}_{r}_{c}" pattern
    keys = {r[1] for r in rendered}
    assert "linked_grid_0_0" in keys and "linked_grid_2_0" in keys


def test_render_linked_grid_warns_on_unknown_component(mock_streamlit, temp_workspace):
    """An unknown comp_name is skipped (on_missing='warn') without raising."""
    fm = _fm(temp_workspace)
    ds = make_deconv_caches(fm)
    build_insight_caches(fm, ds, "flashdeconv")
    builders = make_builders(fm, ds, "flashdeconv")

    def fake_render(component, state_manager, key=None, height=None):
        component._prepare_vue_data(state_manager.get_state_for_vue())
        return None

    with patch("openms_insight.rendering.bridge.render_component", fake_render):
        sm = render_linked_grid(
            [["scan_table", "does_not_exist"]], builders,
            state_key=f"flashdeconv__{ds}",
        )
    assert isinstance(sm, StateManager)
