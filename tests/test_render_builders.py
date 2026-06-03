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
    assert builders["sequence_view"]()._residue_identifier == "aa"

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
    assert args["seriesColumn"] == "isotope"  # break line between isotopes
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
    for h in ("ms1_deconv_heat_map", "ms2_deconv_heat_map",
              "ms1_raw_heatmap", "ms2_raw_heatmap"):
        a = b[h]()._get_component_args()
        assert a["xLabel"] == "Retention Time", h
        assert a["yLabel"] == "Monoisotopic Mass", h


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
