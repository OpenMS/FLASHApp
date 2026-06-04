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

import pandas as pd
import polars as pl
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
    # massIndex == the per-scan ordinal the 3D S/N plot consumes (SignalPeaks[i]);
    # the oracle mass-table click selected the row's within-scan index, NOT a
    # global id, so the "mass" slot must carry mass_in_scan.
    assert mass_table.get_interactivity_mapping() == {"mass": "mass_in_scan"}

    plot3d = builders["3D_SN_plot"]()
    # massIndex -> value filter on mass_in_scan; scanIndex -> scan
    assert plot3d.get_filters_mapping() == {"scan": "scan_id", "mass": "mass_in_scan"}
    # 3D x-axis is the oracle "Mass" = mz*charge (not raw m/z); y=charge, z=intensity
    p3d_args = plot3d._get_component_args()
    assert (p3d_args["xColumn"], p3d_args["yColumn"], p3d_args["zColumn"]) == (
        "mass", "charge", "intensity",
    )
    assert p3d_args["xLabel"] == "Mass"  # Plot3D default matches oracle axis title

    tag_table = builders["tag_table"]()
    # tags are scan (spectrum) data: the oracle filtered by Scan and showed ALL of
    # a scan's tags for ANY proteoform on that scan, so the tag table follows the
    # protein->scan selection via scan_id (not a collapsed per-scan protein_id).
    assert tag_table.get_filters_mapping() == {"scan": "scan_id"}
    assert tag_table.get_interactivity_mapping() == {"tag": "tag_id"}

    # the protein-row click resolves to its scan (value-based proteoform_scan_map):
    # it sets BOTH protein and scan so all scan-keyed panels follow the proteoform.
    protein_table = builders["protein_table"]()
    assert protein_table.get_interactivity_mapping() == {
        "protein": "protein_id", "scan": "scan_id",
    }
    # round-10 finding 3-cascade-001: a protein click also cascade-clears the
    # dependent aa (residue) + tag selections (oracle updateSelectedProtein).
    assert protein_table._get_component_args()["clearsSelections"] == ["aa", "tag"]


def test_tnt_tagger_resolves_tag_payload(mock_streamlit, temp_workspace):
    """The augmented (tagger) spectrum resolves a scalar tag_id (from the tag-table
    click) to the tag's masses/sequence/selectedAA via the tags frame -- the
    value-based replacement for the oracle's opaque TagData payload.
    """
    fm = _fm(temp_workspace)
    ds = make_tnt_caches(fm)
    build_insight_caches(fm, ds, "flashtnt")
    builders = make_builders(fm, ds, "flashtnt")

    tagger = builders["combined_spectrum"]()
    # tag_id 0: Scan 10, TagSequence "PEP", mzs "1,2,3", StartPos 0.
    payload = tagger._resolve_tag_payload(0, {"aa": 2})
    assert payload is not None
    assert payload["sequence"] == "PEP"
    assert payload["masses"] == [1.0, 2.0, 3.0]
    # selectedAA = residue position (aa) - tag StartPos = 2 - 0.
    assert payload["selectedAA"] == 2
    # tag_id 1: StartPos 3 -> selectedAA = 5 - 3 = 2.
    assert tagger._resolve_tag_payload(1, {"aa": 5})["selectedAA"] == 2
    # cleared / unknown selection -> no payload (no crash).
    assert tagger._resolve_tag_payload(None, {}) is None
    assert tagger._resolve_tag_payload(999, {}) is None
    # tag + residue selections drive a re-render.
    deps = tagger.get_state_dependencies()
    assert "tag" in deps and "aa" in deps

    # The SequenceView publishes residue clicks as the "aa" selection the tagger
    # consumes (closing the residue -> selectedAA cross-link).
    # round-12 finding 3-seqview-001: oracle two-path residue click -- PATH 1 aa is
    # coverage-gated + toggling (coverage_column set), PATH 2 publishes the matched
    # fragment's mass_in_scan to "mass" (fragment_mass_identifier="mass").
    sv = builders["sequence_view"]()
    assert sv._residue_identifier == "aa"
    assert sv._coverage_column == "coverage"
    assert sv._fragment_mass_identifier == "mass"
    # round-13 findings 3-seqview-003/004: mass-info header (observed proteoform mass)
    # + inbound mass->fragment-table-row highlight.
    assert sv._observed_mass_column == "observed_mass"
    assert sv._mass_selection_identifier == "mass"
    # round-14 finding 3-seqview-005: oracle proteoform-branch header labels.
    assert sv._theoretical_mass_label == "Theoretical protein mass"
    assert sv._observed_mass_label == "Observed proteoform mass"

    # round-13 finding 3-seqview-002: the FLASHDeconv sequence view (global sequence,
    # no tags/coverage -> PATH 2 only) must ALSO publish the fragment's mass on a
    # residue click (oracle aminoAcidSelected -> updateSelectedMass runs on every tool).
    dfm = _fm(temp_workspace)
    make_deconv_caches(dfm, ds="deconv_seqmass")
    make_sequence_cache(dfm)  # global deconv sequence ("sequence" dataset)
    build_insight_caches(dfm, "deconv_seqmass", "flashdeconv")
    dsv = make_builders(dfm, "deconv_seqmass", "flashdeconv")["sequence_view"]()
    assert dsv._fragment_mass_identifier == "mass"
    assert dsv._coverage_column is None  # no tags on the global deconv sequence
    assert dsv._mass_selection_identifier == "mass"  # inbound mass->fragment highlight
    # round-15 finding 3-seqview-006: deconv shows the PRECURSOR mass-info header
    # (per-scan PrecursorMass -> observed_mass), with the generic "Precursor" title.
    assert dsv._observed_mass_column == "observed_mass"
    assert dsv._mass_header_title == "Precursor"
    # seq_deconv carries per-scan observed_mass (PrecursorMass, NULL for MS1).
    sdf = pl.read_parquet(dfm.result_path("deconv_seqmass", "seq_deconv"))
    assert "observed_mass" in sdf.columns

    # In FLASHDeconv (no tags frame) the tagger has no tag resolution wired.
    dds = make_deconv_caches(_fm(temp_workspace), ds="deconv1")
    fm2 = _fm(temp_workspace)
    build_insight_caches(fm2, "deconv1", "flashdeconv")
    deconv_tagger = make_builders(fm2, "deconv1", "flashdeconv")["combined_spectrum"]()
    assert deconv_tagger._tag_data is None


def test_tnt_residue_narrows_tag_table(mock_streamlit, temp_workspace):
    """Clicking a sequence residue ('aa') narrows the tag table to tags spanning it
    (StartPos <= aa <= EndPos), on top of the scan filter; shows all when unset.
    """
    fm = _fm(temp_workspace)
    ds = make_tnt_caches(fm)
    build_insight_caches(fm, ds, "flashtnt")
    tag_table = make_builders(fm, ds, "flashtnt")["tag_table"]()
    # fixture: scan 0 has tag 0 (StartPos 0,EndPos 2) and tag 1 (StartPos 3,EndPos 5).
    assert "aa" in tag_table.get_state_dependencies()
    both = tag_table._prepare_vue_data({"scan": 0})["tableData"]
    assert sorted(both["tag_id"].tolist()) == [0, 1]
    only0 = tag_table._prepare_vue_data({"scan": 0, "aa": 1})["tableData"]
    assert only0["tag_id"].tolist() == [0]
    only1 = tag_table._prepare_vue_data({"scan": 0, "aa": 4})["tableData"]
    assert only1["tag_id"].tolist() == [1]


def test_quant_3d_axes_match_oracle(mock_streamlit, temp_workspace):
    """Quant feature-trace 3D uses oracle axes: x=m/z, y=RT, z=intensity (labeled),
    drawn as connected per-charge elution lines (stem off), not per-point spikes."""
    fm = _fm(temp_workspace)
    ds = make_quant_caches(fm)
    build_insight_caches(fm, ds, "flashquant")
    p3d = make_builders(fm, ds, "flashquant")["quant_traces_3d"]()
    args = p3d._get_component_args()
    assert (args["xColumn"], args["yColumn"], args["zColumn"]) == ("mz", "rt", "intensity")
    assert args["xLabel"] == "m/z"
    assert args["yLabel"] == "retention time"
    assert args["stem"] is False  # connected elution lines per charge, not spikes
    # oracle FLASHQuantView draws one trace per CHARGE but breaks the polyline
    # between isotopes within a charge (it pushes a -1000 z sentinel before/after
    # each isotope's points) and labels each trace `Charge: ${charge}`.
    assert args["categoryColumn"] == "charge"
    # round-8 finding 3-quant-005: the polyline breaks per ACTUAL trace, not per
    # (charge, isotope) -- two traces of one feature can share (charge, isotope), so
    # keying on "isotope" would merge them. series_column="trace_in_feature" (stable
    # per-feature running trace id) breaks each real trace into its own line.
    assert args["seriesColumn"] == "trace_in_feature"
    assert args["categoryNameTemplate"] == "Charge: {}"  # legend "Charge: 2"


def test_axis_titles_match_oracle(mock_streamlit, temp_workspace):
    """Spectra + heatmaps carry the oracle's human-readable axis titles (not raw
    column names)."""
    fm = _fm(temp_workspace)
    ds = make_deconv_caches(fm)
    build_insight_caches(fm, ds, "flashdeconv")
    b = make_builders(fm, ds, "flashdeconv")

    dec = b["deconv_spectrum"]()._get_component_args()
    assert dec["xLabel"] == "Monoisotopic Mass" and dec["yLabel"] == "Intensity"
    ann = b["anno_spectrum"]()._get_component_args()
    assert ann["xLabel"] == "m/z" and ann["yLabel"] == "Intensity"
    # round-16 finding 3-heatmap-001: deconv heatmaps -> "Monoisotopic Mass";
    # RAW heatmaps -> "m/z" (raw m/z data), matching oracle PlotlyHeatmap yAxisLabel.
    for h in ("ms1_deconv_heat_map", "ms2_deconv_heat_map",
              "ms1_raw_heatmap", "ms2_raw_heatmap"):
        a = b[h]()._get_component_args()
        assert a["xLabel"] == "Retention Time", h
        expected_y = "m/z" if h.endswith("raw_heatmap") else "Monoisotopic Mass"
        assert a["yLabel"] == expected_y, h


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


# --------------------------------------------------------------------------- #
# oracle Tabulator column chrome (titles + formatters + sorters + initialSort)
# --------------------------------------------------------------------------- #
# Ported from TabulatorScanTable / TabulatorMassTable / TabulatorProteinTable /
# TabulatorTagTable.vue + FLASHQuantView.vue. These lock that the migrated Insight
# Tables present the SAME curated columns (titles + number formatters + per-table
# initial sort) and HIDE the internal carrier columns, while keeping the existing
# value-based cross-link wiring (covered by the tests above) intact.
def _col_defs(comp):
    """Displayed column-definition list as it reaches Vue."""
    return comp._get_component_args()["columnDefinitions"]


def _by_title(defs):
    return {c["title"]: c for c in defs}


def _by_field(defs):
    return {c["field"]: c for c in defs}


def test_scan_table_column_chrome(mock_streamlit, temp_workspace):
    """Scan Table: oracle titles/fields, guarded-toFixed on RT/PrecursorMass, no
    initialSort; the per-scan ordinal carrier (mass_in_scan) is not displayed."""
    fm = _fm(temp_workspace)
    ds = make_deconv_caches(fm)
    build_insight_caches(fm, ds, "flashdeconv")
    defs = _col_defs(make_builders(fm, ds, "flashdeconv")["scan_table"]())

    bt = _by_title(defs)
    # title -> field parity (oracle "Index" maps to the schema id column scan_id)
    assert bt["Index"]["field"] == "scan_id"
    assert bt["Scan Number"]["field"] == "Scan"
    assert bt["MS Level"]["field"] == "MSLevel"
    assert bt["Retention time"]["field"] == "RT"
    assert bt["Precursor Mass"]["field"] == "PrecursorMass"
    assert bt["#Masses"]["field"] == "#Masses"
    # toFixedFormatter() -> the guarded "fixed" named formatter
    assert bt["Retention time"]["formatter"] == "fixed"
    assert bt["Retention time"]["formatterParams"] == {"precision": 4, "minLength": 4}
    assert bt["Precursor Mass"]["formatter"] == "fixed"
    # exactly the oracle's 6 columns, in order; no carriers (mass_in_scan) shown
    shown = [c["field"] for c in defs]
    assert shown == ["scan_id", "Scan", "MSLevel", "RT", "PrecursorMass", "#Masses"]
    assert "mass_in_scan" not in shown


def test_mass_table_column_chrome(mock_streamlit, temp_workspace):
    """Mass Table: oracle titles, fixed formatter on the 5 score/mass columns; the
    interactivity carrier (mass_in_scan) stays in the data but is not displayed."""
    fm = _fm(temp_workspace)
    ds = make_deconv_caches(fm)
    build_insight_caches(fm, ds, "flashdeconv")
    defs = _col_defs(make_builders(fm, ds, "flashdeconv")["mass_table"]())

    bt = _by_title(defs)
    assert bt["Index"]["field"] == "mass_id"
    assert bt["Monoisotopic mass"]["field"] == "MonoMass"
    assert bt["Sum intensity"]["field"] == "SumIntensity"
    assert bt["Min charge"]["field"] == "MinCharges"
    assert bt["Max charge"]["field"] == "MaxCharges"
    assert bt["Min isotope"]["field"] == "MinIsotopes"
    assert bt["Max isotope"]["field"] == "MaxIsotopes"
    # the five toFixed'd columns carry the "fixed" formatter
    for title in ("Monoisotopic mass", "Sum intensity", "Cosine score", "SNR", "QScore"):
        assert bt[title]["formatter"] == "fixed", title
        assert bt[title]["formatterParams"] == {"precision": 4, "minLength": 4}
    # charge/isotope columns are plain (no formatter), matching the oracle
    assert "formatter" not in bt["Min charge"]
    # carrier hidden
    assert "mass_in_scan" not in {c["field"] for c in defs}


def test_protein_table_column_chrome(mock_streamlit, temp_workspace):
    """Protein Table: oracle titles, -1->'-' placeholder on Mass/Q-Value, initialSort
    Score desc; Coverage(%) (commented out in the oracle) is omitted; the protein_id
    / scan_id carriers (cross-link) are not displayed (no 'Index' column)."""
    fm = _fm(temp_workspace)
    ds = make_tnt_caches(fm)
    build_insight_caches(fm, ds, "flashtnt")
    comp = make_builders(fm, ds, "flashtnt")["protein_table"]()
    defs = _col_defs(comp)

    bt = _by_title(defs)
    assert bt["Scan No."]["field"] == "Scan"
    assert bt["Accession"]["field"] == "accession"
    assert bt["Description"]["field"] == "description"
    assert bt["Length"]["field"] == "length"
    assert bt["Mass"]["field"] == "ProteoformMass"
    assert bt["No. of Matched Fragments"]["field"] == "MatchingFragments"
    assert bt["No. of Modifications"]["field"] == "ModCount"
    assert bt["No. of Tags"]["field"] == "TagCount"
    assert bt["Score"]["field"] == "Score"
    assert bt["Q-Value (Proteoform Level)"]["field"] == "ProteoformLevelQvalue"
    # inline -1 -> '-' becomes the "placeholder" named formatter
    assert bt["Mass"]["formatter"] == "placeholder"
    assert bt["Mass"]["formatterParams"] == {
        "sentinels": [-1], "text": "-", "loose": True,
    }
    assert bt["Q-Value (Proteoform Level)"]["formatter"] == "placeholder"
    # initialSort ported verbatim (Score desc)
    assert comp._get_component_args()["initialSort"] == [{"column": "Score", "dir": "desc"}]
    # Coverage(%) is commented out in the oracle -> not displayed; carriers hidden
    shown = {c["field"] for c in defs}
    assert "Coverage(%)" not in shown
    assert "protein_id" not in shown and "scan_id" not in shown
    # no synthetic "Index" column on the protein table (oracle leads with Scan No.)
    assert "Index" not in {c["title"] for c in defs}


def test_tag_table_column_chrome(mock_streamlit, temp_workspace):
    """Tag Table: oracle titles, -1->'-' placeholder on N mass / C mass, initialSort
    Score desc; StartPos/EndPos ARE displayed (and drive the residue interval filter)
    while tag_id / mzs / ProteinIndex carriers are not displayed."""
    fm = _fm(temp_workspace)
    ds = make_tnt_caches(fm)
    build_insight_caches(fm, ds, "flashtnt")
    comp = make_builders(fm, ds, "flashtnt")["tag_table"]()
    defs = _col_defs(comp)

    bt = _by_title(defs)
    assert bt["Scan Number"]["field"] == "Scan"
    assert bt["Start Position"]["field"] == "StartPos"
    assert bt["End Position"]["field"] == "EndPos"
    assert bt["Sequence"]["field"] == "TagSequence"
    assert bt["Length"]["field"] == "Length"
    assert bt["Tag Score"]["field"] == "Score"
    assert bt["N mass"]["field"] == "Nmass"
    assert bt["C mass"]["field"] == "Cmass"
    # the unicode Delta title is preserved verbatim
    assert "Δ mass" in bt and bt["Δ mass"]["field"] == "DeltaMass"
    # N mass / C mass use the -1 -> '-' placeholder; Delta mass is plain
    assert bt["N mass"]["formatter"] == "placeholder"
    assert bt["C mass"]["formatter"] == "placeholder"
    assert "formatter" not in bt["Δ mass"]
    assert comp._get_component_args()["initialSort"] == [{"column": "Score", "dir": "desc"}]
    shown = {c["field"] for c in defs}
    # StartPos/EndPos shown (also the interval-filter bounds); carriers hidden
    assert {"StartPos", "EndPos"} <= shown
    assert not ({"tag_id", "mzs", "ProteinIndex"} & shown)


def test_tag_table_placeholder_renders_dash_data(mock_streamlit, temp_workspace):
    """The N mass / C mass placeholder columns carry the -1 sentinel data the
    formatter renders as '-' (fixture has Nmass=-1 on tag 0, Cmass=-1 on tag 1)."""
    fm = _fm(temp_workspace)
    ds = make_tnt_caches(fm)
    build_insight_caches(fm, ds, "flashtnt")
    comp = make_builders(fm, ds, "flashtnt")["tag_table"]()
    rows = comp._prepare_vue_data({"scan": 0})["tableData"]
    # both Nmass and Cmass are projected (displayed) and carry the -1 sentinel
    assert "Nmass" in rows.columns and "Cmass" in rows.columns
    assert -1.0 in rows["Nmass"].tolist()
    assert -1.0 in rows["Cmass"].tolist()


def test_quant_feature_table_column_chrome(mock_streamlit, temp_workspace):
    """FLASHQuant feature table: oracle titles, FWHM RT fields mapped to the schema's
    StartRT/EndRT, no formatters, no initialSort; the duplicate 'Feature Group
    Quantity' from the oracle is de-duplicated to one column."""
    fm = _fm(temp_workspace)
    ds = make_quant_caches(fm)
    build_insight_caches(fm, ds, "flashquant")
    comp = make_builders(fm, ds, "flashquant")["quant_visualization"]()
    args = comp._get_component_args()
    defs = args["columnDefinitions"]

    bt = _by_title(defs)
    assert bt["Index"]["field"] == "feature_id"
    assert bt["Monoisotopic Mass"]["field"] == "MonoisotopicMass"
    assert bt["Average Mass"]["field"] == "AverageMass"
    # oracle StartRetentionTime(FWHM)/EndRetentionTime(FWHM) -> schema StartRT/EndRT
    assert bt["Start Retention Time (FWHM)"]["field"] == "StartRT"
    assert bt["End Retention Time (FWHM)"]["field"] == "EndRT"
    assert bt["Most Abundant Charge"]["field"] == "MostAbundantFeatureCharge"
    assert bt["Isotope Cosine Score"]["field"] == "IsotopeCosineScore"
    # no number formatters in the oracle quant table
    assert all("formatter" not in c for c in defs)
    # no initialSort for the quant table
    assert "initialSort" not in args
    # the oracle's duplicate "Feature Group Quantity" collapses to a single column
    assert [c["title"] for c in defs].count("Feature Group Quantity") == 1


# --------------------------------------------------------------------------- #
# round-8 wiring findings (selective highlight / dynamic title / per-trace 3D /
# go-to fields / FDR chrome / feature-group title / best-per-spectrum toggle)
# --------------------------------------------------------------------------- #
def test_anno_spectrum_selective_highlight_wiring(mock_streamlit, temp_workspace):
    """finding 3-anno-001: the annotated spectrum drops the static is_signal
    highlight and instead highlights the SELECTED mass's signal peaks via the
    highlight LINK frame (z=N labels + deconv-peaks toggle). It MUST expose peak_id
    as its first interactivity column so the link key-set maps onto drawn peaks."""
    fm = _fm(temp_workspace)
    ds = make_tnt_caches(fm)
    build_insight_caches(fm, ds, "flashtnt")
    comp = make_builders(fm, ds, "flashtnt")["anno_spectrum"]()

    # static is_signal highlight is REMOVED.
    assert comp._highlight_column is None
    args = comp._get_component_args()
    assert args["highlightColumn"] is None
    # selection-driven LINK highlight params.
    assert comp._highlight_selection == "mass"
    assert comp._highlight_link_path == fm.result_path(ds, "anno_highlight_link")
    assert comp._highlight_link_key_column == "peak_id"
    assert comp._highlight_link_match_column == "mass_in_scan"
    assert comp._highlight_charge_column == "charge"
    assert comp._highlight_annotation_template == "z={}"
    assert comp._deconv_peaks_toggle is True
    # peak_id is the FIRST (only) interactivity column == the highlight id_column
    # (lineplot keys the highlight key-set off list(interactivity.values())[0]); the
    # private "anno_peak" slot is NOT the shared "mass" slot (parity-bug fix kept).
    assert args["interactivity"] == {"anno_peak": "peak_id"}
    assert list(comp.get_interactivity_mapping().values())[0] == "peak_id"
    assert "mass" not in comp.get_interactivity_mapping()
    # the selective-highlight modebar wiring is enabled with the deconv-peaks toggle.
    assert args["selectiveHighlightEnabled"] is True
    assert args["deconvPeaksToggle"] is True
    # the highlight is a state dependency on "mass" (selection change -> recompute).
    assert "mass" in comp.get_state_dependencies()


def test_anno_spectrum_highlight_maps_onto_peaks(mock_streamlit, temp_workspace):
    """Selecting a mass highlights that mass's signal peaks on the annotated
    spectrum (the link key-set maps onto drawn peaks via peak_id) and emits the
    client-side toggle payload keyed on peak_id."""
    fm = _fm(temp_workspace)
    ds = make_tnt_caches(fm)
    build_insight_caches(fm, ds, "flashtnt")
    comp = make_builders(fm, ds, "flashtnt")["anno_spectrum"]()

    link = pl.read_parquet(fm.result_path(ds, "anno_highlight_link"))
    assert link.height > 0
    row = link.row(0, named=True)
    vd = comp._prepare_vue_data({"scan": row["scan_id"], "mass": row["mass_in_scan"]})
    hl_col = vd["_plotConfig"]["highlightColumn"]
    pdf = vd["plotData"]
    # at least one annotated peak is highlighted for the selected mass.
    assert hl_col in pdf.columns and int(pdf[hl_col].sum()) >= 1
    # the client-side toggle payload keys on peak_id and exposes the all-signal set.
    sh = vd["selectiveHighlight"]
    assert sh["idColumn"] == "peak_id"
    assert isinstance(sh["allSignalKeys"], list)
    assert sh["deconvPeaksToggle"] is True
    # with NO mass selected, nothing is highlighted (selection-driven).
    vd0 = comp._prepare_vue_data({"scan": row["scan_id"]})
    hl0 = vd0["_plotConfig"]["highlightColumn"]
    pdf0 = vd0["plotData"]
    assert hl0 not in pdf0.columns or int(pdf0[hl0].sum()) == 0


def test_deconv_spectrum_selective_highlight_wiring(mock_streamlit, temp_workspace):
    """deconv selective highlight: the SELECTED mass's stick highlights via the
    match-column path (no link frame, no z=N labels, no deconv-peaks toggle)."""
    fm = _fm(temp_workspace)
    ds = make_deconv_caches(fm)
    build_insight_caches(fm, ds, "flashdeconv")
    comp = make_builders(fm, ds, "flashdeconv")["deconv_spectrum"]()

    assert comp._highlight_selection == "mass"
    assert comp._highlight_match_column == "mass_in_scan"
    # no link frame on the deconv spectrum (match-column path) => no z=N labels.
    assert comp._highlight_link_path is None
    args = comp._get_component_args()
    assert args["selectiveHighlightEnabled"] is True
    # NO "Show Deconvolved Peaks" toggle for the deconvolved spectrum (oracle parity).
    assert args["deconvPeaksToggle"] is False
    # clicking still selects the shared "mass" slot.
    assert comp.get_interactivity_mapping() == {"mass": "mass_in_scan"}

    # round-9 finding 3-deconv-001: deconv draws the selected mass's MonoMass
    # value label (oracle mass.toFixed(2)) via the match-column value producer.
    assert comp._highlight_value_column == "mass"
    assert comp._highlight_value_template == "{:.2f}"

    # functional: selecting a mass highlights that mass's stick.
    dft = pl.read_parquet(fm.result_path(ds, "deconv_spectrum_tidy"))
    r = dft.row(0, named=True)
    vd = comp._prepare_vue_data({"scan": r["scan_id"], "mass": r["mass_in_scan"]})
    hl_col = vd["_plotConfig"]["highlightColumn"]
    pdf = vd["plotData"]
    assert hl_col in pdf.columns and bool(pdf[hl_col].any())
    # ... and draws exactly one MonoMass value label at that stick.
    anns = vd["peakAnnotations"]
    assert len(anns) == 1
    assert anns[0]["text"] == f"{r['mass']:.2f}"


def test_3d_sn_plot_dynamic_title(mock_streamlit, temp_workspace):
    """finding 3-3d-001: the 3D S/N plot has a dynamic title driven by the SAME
    scan/mass identifiers its filters use: '' (no scan) / 'Precursor signals'
    (scan, no mass) / 'Mass signals' (mass selected)."""
    fm = _fm(temp_workspace)
    ds = make_deconv_caches(fm)
    build_insight_caches(fm, ds, "flashdeconv")
    comp = make_builders(fm, ds, "flashdeconv")["3D_SN_plot"]()

    # title_selection uses the filters' identifier names ("scan"/"mass").
    assert comp._get_component_args()["titleSelection"] == {"scan": "scan", "mass": "mass"}
    assert comp.compute_dynamic_title({}) == ""
    assert comp.compute_dynamic_title({"scan": 0}) == "Precursor signals"
    assert comp.compute_dynamic_title({"scan": 0, "mass": 1}) == "Mass signals"


def test_scan_table_resets_mass_on_scan_change(mock_streamlit, temp_workspace):
    """finding 3-cascade-002: a scan-table click resets the mass selection to the
    new scan's FIRST mass (oracle updateSelectedScan -> updateSelectedMass(0)). The
    scan_table cascade-clears "mass"; the mass_table (default_row=0) then re-defaults
    to mass_in_scan 0 of the selected scan via the bridge _auto_selection, so a stale
    per-scan ordinal cannot carry across a scan switch."""
    fm = _fm(temp_workspace)
    ds = make_deconv_caches(fm)
    build_insight_caches(fm, ds, "flashdeconv")
    b = make_builders(fm, ds, "flashdeconv")

    # scan click cascade-clears the dependent mass selection.
    assert b["scan_table"]()._get_component_args()["clearsSelections"] == ["mass"]

    # with mass unset (the post-clear state), the mass_table auto-selects the first
    # mass (mass_in_scan 0) of the selected scan -> equals the oracle's mass=0 reset.
    mt = b["mass_table"]()
    vd = mt._prepare_vue_data({"scan": 1})
    assert vd.get("_auto_selection", {}).get("mass") == 0


def test_quant_traces_3d_per_trace_break(mock_streamlit, temp_workspace):
    """finding 3-quant-005: the quant 3D breaks its polyline per ACTUAL trace
    (series_column="trace_in_feature"), keeping per-charge color/legend."""
    fm = _fm(temp_workspace)
    ds = make_quant_caches(fm)
    build_insight_caches(fm, ds, "flashquant")
    args = make_builders(fm, ds, "flashquant")["quant_traces_3d"]()._get_component_args()
    assert args["seriesColumn"] == "trace_in_feature"
    assert args["categoryColumn"] == "charge"
    assert args["categoryNameTemplate"] == "Charge: {}"
    # the per-trace id is present in the traces frame so the break is real.
    traces = pl.read_parquet(fm.result_path(ds, "quant_traces"))
    assert "trace_in_feature" in traces.columns


def test_table_go_to_fields_match_oracle(mock_streamlit, temp_workspace):
    """finding 3-tables-003: each Table passes the oracle's explicit goToFields so
    auto-detect never exposes internal carrier columns (scan_id-as-mass_in_scan,
    protein_id, tag_id, etc.). The FLASHQuant feature table disables go-to ([])."""
    fm = _fm(temp_workspace)
    tnt = make_tnt_caches(fm, ds="gtf_tnt")
    build_insight_caches(fm, tnt, "flashtnt")
    b = make_builders(fm, tnt, "flashtnt")

    # scan/mass: oracle ['id','Scan'] / ['id'] -> schema id columns scan_id/mass_id.
    assert b["scan_table"]()._get_component_args()["goToFields"] == ["scan_id", "Scan"]
    assert b["mass_table"]()._get_component_args()["goToFields"] == ["mass_id"]
    # protein/tag: oracle lists verbatim; carriers (protein_id/tag_id) excluded.
    assert b["protein_table"]()._get_component_args()["goToFields"] == ["Scan", "accession"]
    assert b["tag_table"]()._get_component_args()["goToFields"] == [
        "Scan", "StartPos", "EndPos", "TagSequence",
    ]
    # carriers are not exposed as go-to fields.
    for name, carriers in (
        ("scan_table", {"mass_in_scan"}),
        ("mass_table", {"mass_in_scan", "scan_id"}),
        ("protein_table", {"protein_id", "scan_id"}),
        ("tag_table", {"tag_id", "scan_id"}),
    ):
        gtf = set(b[name]()._get_component_args()["goToFields"])
        assert not (gtf & carriers), name

    # FLASHQuant feature table: oracle had no go-to-fields -> disabled with [] (so
    # goToFields is NOT emitted to Vue, vs auto-detect exposing feature_id etc.).
    qfm = _fm(temp_workspace)
    qds = make_quant_caches(qfm, ds="gtf_quant")
    build_insight_caches(qfm, qds, "flashquant")
    qargs = make_builders(qfm, qds, "flashquant")["quant_visualization"]()._get_component_args()
    assert "goToFields" not in qargs


def test_fdr_plots_oracle_title_and_trace_labels(mock_streamlit, temp_workspace):
    """findings 3-fdr-001/002: both FDR density plots use title "FDR Plot" and the
    oracle trace legend names "Target QScores" / "Decoy QScores"."""
    fm = _fm(temp_workspace)
    # flashdeconv -> fdr_plot
    dds = make_deconv_caches(fm, ds="fdr_d")
    build_insight_caches(fm, dds, "flashdeconv")
    fdr = make_builders(fm, dds, "flashdeconv")["fdr_plot"]()
    fargs = fdr._get_component_args()
    assert fargs["title"] == "FDR Plot"
    assert fargs["targetLabel"] == "Target QScores"
    assert fargs["decoyLabel"] == "Decoy QScores"

    # flashtnt -> id_fdr_plot
    tds = make_tnt_caches(_fm(temp_workspace), ds="fdr_t")
    fm2 = _fm(temp_workspace)
    build_insight_caches(fm2, "fdr_t", "flashtnt")
    idfdr = make_builders(fm2, "fdr_t", "flashtnt")["id_fdr_plot"]()
    iargs = idfdr._get_component_args()
    assert iargs["title"] == "FDR Plot"
    assert iargs["targetLabel"] == "Target QScores"
    assert iargs["decoyLabel"] == "Decoy QScores"


def test_quant_feature_table_title_feature_groups(mock_streamlit, temp_workspace):
    """finding 3-feat-001: the FLASHQuant feature table title is "Feature groups"
    (oracle FLASHQuantView), not "Features"."""
    fm = _fm(temp_workspace)
    ds = make_quant_caches(fm)
    build_insight_caches(fm, ds, "flashquant")
    args = make_builders(fm, ds, "flashquant")["quant_visualization"]()._get_component_args()
    assert args["title"] == "Feature groups"


def _multi_proteoform_tnt(fm, ds):
    """Build tnt caches whose protein frame has TWO proteoforms on ONE Scan (so the
    best-per-spectrum flag actually distinguishes them) + one on another Scan."""
    make_tnt_caches(fm, ds=ds)
    # Scan 10: proteoforms with Score 5 and 9 (best = 9); Scan 20: a single one.
    protein_df = pd.DataFrame({
        "index": [0, 1, 2], "accession": ["P1", "P1b", "P2"],
        "description": ["d", "d", "d"],
        "sequence": ["PEPTIDEK", "PEPTIDEK", "ACDEFGHK"],
        "length": [8, 8, 8], "ProteoformMass": [900.4, 900.4, 800.3],
        "MatchingFragments": [12, 3, 8], "Coverage(%)": [55.0, 10.0, 40.0],
        "ModCount": [0, 0, 1], "TagCount": [2, 1, 1], "Score": [5.0, 9.0, 6.0],
        "ProteoformLevelQvalue": [0.01, 0.02, 0.5], "Scan": [10, 10, 20],
    })
    fm.store_data(ds, "protein_dfs", protein_df)
    build_insight_caches(fm, ds, "flashtnt", regenerate=True)


def test_protein_best_per_spectrum_toggle(mock_streamlit, temp_workspace):
    """finding 3-tables-002: best_per_spectrum=True sources the is_best_per_scan==1
    subset under a DISTINCT cache_id (so the toggle reliably swaps the cached row
    set); False sources the full table under the normal cache_id. Column chrome /
    interactivity / index_field / initial_sort stay identical across both."""
    fm = _fm(temp_workspace)
    _multi_proteoform_tnt(fm, "bps")

    best = make_builders(fm, "bps", "flashtnt", best_per_spectrum=True)["protein_table"]()
    allp = make_builders(fm, "bps", "flashtnt", best_per_spectrum=False)["protein_table"]()

    # DISTINCT cache_ids so the two row sets cache independently (toggle swap).
    assert best._cache_id == "flashtnt__bps__protein_table_best"
    assert allp._cache_id == "flashtnt__bps__protein_table"
    assert best._cache_id != allp._cache_id

    # filtered (best) shows one row per Scan (the highest Score); full shows all 3.
    best_rows = best._prepare_vue_data({})["tableData"]
    all_rows = allp._prepare_vue_data({})["tableData"]
    assert len(best_rows) == 2  # Scan 10 (best proteoform) + Scan 20
    assert len(all_rows) == 3
    # the kept Scan-10 proteoform is the higher-Score one (protein_id 1, Score 9).
    assert sorted(best_rows["protein_id"].tolist()) == [1, 2]

    # column chrome / interactivity / index / initial_sort are IDENTICAL.
    bargs, aargs = best._get_component_args(), allp._get_component_args()
    assert bargs["columnDefinitions"] == aargs["columnDefinitions"]
    assert bargs["interactivity"] == aargs["interactivity"] == {
        "protein": "protein_id", "scan": "scan_id",
    }
    assert bargs["tableIndexField"] == aargs["tableIndexField"] == "protein_id"
    assert bargs["initialSort"] == aargs["initialSort"] == [{"column": "Score", "dir": "desc"}]
    assert bargs["goToFields"] == aargs["goToFields"] == ["Scan", "accession"]

    # default wiring (no kwarg) is best-per-spectrum (oracle default ON).
    default = make_builders(fm, "bps", "flashtnt")["protein_table"]()
    assert default._cache_id == "flashtnt__bps__protein_table_best"
    assert len(default._prepare_vue_data({})["tableData"]) == 2


def test_best_per_spectrum_preserves_scan_cross_link(mock_streamlit, temp_workspace):
    """Both protein-table row sets carry scan_id, so the downstream scan-keyed
    panels (tag table / sequence view / augmented spectrum) cross-link unchanged
    regardless of the toggle."""
    fm = _fm(temp_workspace)
    _multi_proteoform_tnt(fm, "bps2")
    for flag in (True, False):
        rows = make_builders(
            fm, "bps2", "flashtnt", best_per_spectrum=flag
        )["protein_table"]()._prepare_vue_data({})["tableData"]
        # scan_id carrier present (drives the protein->scan cross-link) in both sets.
        assert "scan_id" in rows.columns
        assert rows["scan_id"].notna().all()
