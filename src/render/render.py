"""FLASHApp's OpenMS-Insight builder factory (post Phase-3 migration).

This module is repurposed from the old bespoke-Vue grid-render loop
(``render_grid`` / ``render_component`` + ``StateTracker``) to a thin **builder
factory**. The grid itself now comes from the frozen, tool-agnostic template
module ``src.view.grid`` (``render_linked_grid`` + ``LayoutManager``); the viewer
pages import that and feed it the builders produced here.

``make_builders(file_manager, dataset_id, tool, settings=None)`` returns a
``{comp_name: () -> BaseComponent}`` map. Each zero-arg factory closes over
``dataset_id`` + ``file_manager`` + an Insight cache dir and uses
``file_manager.result_path(...)`` (the tidy parquet written by
``src.render.schema.build_insight_caches``) to feed ``data_path=``. ``cache_id``
is ``f"{tool}__{dataset_id}__{comp_name}"`` so component caches are per-dataset
-- this is the oracle's "dataset changed -> reset" guarantee expressed through
``cache_id`` (the StateManager is likewise scoped per ``(tool, experiment)`` via
``state_key`` inside ``render_linked_grid``).

The OLD index-based selection maps to value-based ``filters`` / ``interactivity``
(see ``migration/specs/PHASE3_PLAN.md`` 5.3 and the deleted ``update.py``):

==========================  ============================================
oracle (index-based)        insight (value-based)
==========================  ============================================
``scanIndex`` / iloc        selection ``scan`` = ``scan_id``; ``filters={"scan":"scan_id"}``
``massIndex`` / ``[idx]``    selection ``mass`` = ``mass_in_scan`` (per-scan ordinal;
                            the table/deconv-spectrum/3D all share this slot)
``proteinIndex`` + scan_map protein-row click sets ``protein`` = ``protein_id`` AND
                            ``scan`` = ``scan_id`` (denormalized deconv_index); the
                            scan-keyed panels (tag table, augmented spectrum,
                            sequence-view peaks) follow via ``filters={"scan":...}``
heatmap ``xRange/yRange``    Heatmap internal zoom (per-instance ``zoom_identifier``)
``StateTracker``            ``StateManager(session_key=state_key)``
==========================  ============================================
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
from openms_insight import Heatmap, LinePlot, Plot3D, SequenceView, Table


def _insight_cache_dir(file_manager) -> str:
    """Keep Insight's own disk caches under the workspace cache dir."""
    return str(Path(file_manager.cache_path, "insight"))


# --------------------------------------------------------------------------- #
# Oracle Tabulator column chrome (titles + formatters + sorters + initialSort)
# --------------------------------------------------------------------------- #
# Ported verbatim from the oracle Tabulator{Scan,Mass,Protein,Tag}Table.vue and
# FLASHQuantView.vue ``columnDefinitions`` arrays so the migrated Insight Tables
# show the SAME curated subset of columns with the SAME human titles, number
# formatting and per-table initial sort -- instead of the auto-generated raw
# column names + internal carrier columns. The Table renders ONLY these columns
# (carriers like scan_id / mzs / ProteinIndex stay in the data for
# filters/interactivity/index but are not listed, hence not shown).
#
# Formatter mapping (see OpenMS-Insight Table.with_fixed_format / with_placeholder
# and tabulator-formatters.ts):
#   oracle ``toFixedFormatter()``  -> {"formatter": "fixed",
#                                      "formatterParams": {"precision": 4,
#                                                          "minLength": 4}}
#     (guarded toFixed: only reformats when the value's string length exceeds
#      minLength, matching ``value.toString().length > 4 ? value.toFixed(4) :
#      value``).
#   oracle inline ``value == -1 ? '-' : value`` -> {"formatter": "placeholder",
#       "formatterParams": {"sentinels": [-1], "text": "-", "loose": True}}.
#     None of the oracle -1->"-" columns ALSO toFixed (they return the raw value
#     otherwise), so a plain placeholder is an exact match (no combine nuance).
#
# Field-name mapping (oracle field -> schema column, from src/render/schema.py):
#   * oracle ``id`` ("Index") -> the schema id column (scan_id / mass_id /
#     feature_id); the oracle set row.id = row.index client-side.
#   * FLASHQuant ``StartRetentionTime(FWHM)`` / ``EndRetentionTime(FWHM)`` ->
#     schema ``StartRT`` / ``EndRT`` (renamed by schema._QUANT_SCALAR_RENAME).
#   * all other oracle fields keep their name in the corresponding tidy frame
#     (verified against the real protein.tsv / tags.tsv FLASHTagger headers).
_FIXED_FMT = {"formatter": "fixed", "formatterParams": {"precision": 4, "minLength": 4}}
_DASH_FMT = {
    "formatter": "placeholder",
    "formatterParams": {"sentinels": [-1], "text": "-", "loose": True},
}

# Scan Table (TabulatorScanTable.vue) -- no initialSort.
_SCAN_COLUMN_DEFS = [
    {"field": "scan_id", "title": "Index", "sorter": "number",
     "headerTooltip": "The sequential index of the spectrum in the dataset."},
    {"field": "Scan", "title": "Scan Number", "sorter": "number",
     "headerTooltip": "The identifier of the mass spectrometry scan."},
    {"field": "MSLevel", "title": "MS Level", "sorter": "number",
     "headerTooltip": "The level of mass spectrometry analysis (e.g., MS1 or MS2)."},
    {"field": "RT", "title": "Retention time", "sorter": "number", **_FIXED_FMT,
     "headerTooltip": "The time at which the spectrum was detected during the "
                      "chromatographic separation in seconds."},
    {"field": "PrecursorMass", "title": "Precursor Mass", "sorter": "number", **_FIXED_FMT,
     "headerTooltip": "The mass of the precursor ion selected for fragmentation "
                      "in Daltons."},
    {"field": "#Masses", "title": "#Masses", "sorter": "number",
     "headerTooltip": "The number of detected masses in the spectrum."},
]

# Mass Table (TabulatorMassTable.vue) -- no initialSort.
_MASS_COLUMN_DEFS = [
    {"field": "mass_id", "title": "Index", "sorter": "number",
     "headerTooltip": "The sequential index of the mass entry in the dataset."},
    {"field": "MonoMass", "title": "Monoisotopic mass", "sorter": "number", **_FIXED_FMT,
     "headerTooltip": "The monoisotopic mass of the detected ion in Daltons."},
    {"field": "SumIntensity", "title": "Sum intensity", "sorter": "number", **_FIXED_FMT,
     "headerTooltip": "The total intensity of the detected mass across all "
                      "isotopic peaks and charges."},
    {"field": "MinCharges", "title": "Min charge", "sorter": "number",
     "headerTooltip": "The minimum charge state detected for the mass."},
    {"field": "MaxCharges", "title": "Max charge", "sorter": "number",
     "headerTooltip": "The maximum charge state detected for the mass."},
    {"field": "MinIsotopes", "title": "Min isotope", "sorter": "number",
     "headerTooltip": "The smallest observed isotopic shift, expressed as a "
                      "multiple of the average isotopic mass difference at 55kDA."},
    {"field": "MaxIsotopes", "title": "Max isotope", "sorter": "number",
     "headerTooltip": "The largest observed isotopic shift, expressed as a "
                      "multiple of the average isotopic mass difference at 55kDA."},
    {"field": "CosineScore", "title": "Cosine score", "sorter": "number", **_FIXED_FMT,
     "headerTooltip": "The cosine similarity score comparing the observed and "
                      "theoretical isotopic patterns."},
    {"field": "SNR", "title": "SNR", "sorter": "number", **_FIXED_FMT,
     "headerTooltip": "The signal-to-noise ratio for the detected mass."},
    {"field": "QScore", "title": "QScore", "sorter": "number", **_FIXED_FMT,
     "headerTooltip": "The quality score indicating the confidence of the mass "
                      "detection (higher is better)."},
]

# Protein Table (TabulatorProteinTable.vue) -- initialSort Score desc.
# Coverage(%) is COMMENTED OUT in the oracle, so it is intentionally omitted
# here (all other oracle ProteinTable fields exist in the real protein.tsv).
_PROTEIN_COLUMN_DEFS = [
    {"field": "Scan", "title": "Scan No.", "sorter": "number",
     "headerTooltip": "The identifier of the mass spectrometry scan associated "
                      "with the identified proteoform."},
    {"field": "accession", "title": "Accession",
     "headerTooltip": "The unique identifier for the protein in the reference "
                      "database."},
    {"field": "description", "title": "Description", "responsive": 10},
    {"field": "length", "title": "Length", "responsive": 6, "sorter": "number",
     "headerTooltip": "The total number of amino acids in the matched protein."},
    {"field": "ProteoformMass", "title": "Mass", "responsive": 8, "sorter": "number",
     **_DASH_FMT,
     "headerTooltip": "The calculated mass of the proteoform in Daltons."},
    {"field": "MatchingFragments", "title": "No. of Matched Fragments", "sorter": "number",
     "headerTooltip": "The number of fragment ions that match the protein sequence."},
    {"field": "ModCount", "title": "No. of Modifications", "sorter": "number",
     "headerTooltip": "The number of modifications identified in the protein."},
    {"field": "TagCount", "title": "No. of Tags", "sorter": "number",
     "headerTooltip": "The number of sequence tags associated with the proteoform "
                      "match."},
    {"field": "Score", "title": "Score", "sorter": "number",
     "headerTooltip": "A score indicating the confidence of the protein match "
                      "(higher is better)."},
    {"field": "ProteoformLevelQvalue", "title": "Q-Value (Proteoform Level)",
     "sorter": "number", **_DASH_FMT,
     "headerTooltip": "The confidence value of the protein match at the proteoform "
                      "level."},
]
_PROTEIN_INITIAL_SORT = [{"column": "Score", "dir": "desc"}]

# Tag Table (TabulatorTagTable.vue) -- initialSort Score desc.
_TAG_COLUMN_DEFS = [
    {"field": "Scan", "title": "Scan Number", "sorter": "number",
     "headerTooltip": "The identifier of the mass spectrometry scan containing the "
                      "sequence tag."},
    {"field": "StartPos", "title": "Start Position", "sorter": "number",
     "headerTooltip": "The position in the protein sequence where the sequence tag "
                      "begins."},
    {"field": "EndPos", "title": "End Position", "sorter": "number",
     "headerTooltip": "The position in the protein sequence where the sequence tag "
                      "ends."},
    {"field": "TagSequence", "title": "Sequence", "sorter": "number",
     "headerTooltip": "The amino acid sequence of the identified tag."},
    {"field": "Length", "title": "Length", "sorter": "number",
     "headerTooltip": "The number of amino acids in the sequence tag."},
    {"field": "Score", "title": "Tag Score", "sorter": "number",
     "headerTooltip": "A score indicating the confidence of the sequence tag "
                      "identification (higher is better)."},
    {"field": "Nmass", "title": "N mass", "sorter": "number", **_DASH_FMT,
     "headerTooltip": "The N-terminal mass offset from the start of the sequence "
                      "tag in Daltons."},
    {"field": "Cmass", "title": "C mass", "sorter": "number", **_DASH_FMT,
     "headerTooltip": "The C-terminal mass offset from the end of the sequence tag "
                      "in Daltons."},
    {"field": "DeltaMass", "title": "Δ mass", "sorter": "number",
     "headerTooltip": "Delta mass is the difference between the tag flanking mass "
                      "and the (partial) proteoform mass, from its terminal to the "
                      "tag boundary."},
]
_TAG_INITIAL_SORT = [{"column": "Score", "dir": "desc"}]

# FLASHQuant feature table (FLASHQuantView.vue featureGroupTableColumnDefinitions)
# -- no initialSort, no formatters. The oracle listed "Feature Group Quantity"
# twice (a copy-paste bug); we keep a single definition. StartRetentionTime(FWHM)
# / EndRetentionTime(FWHM) map to the schema's renamed StartRT / EndRT.
_QUANT_COLUMN_DEFS = [
    {"field": "feature_id", "title": "Index", "sorter": "number"},
    {"field": "MonoisotopicMass", "title": "Monoisotopic Mass", "sorter": "number"},
    {"field": "AverageMass", "title": "Average Mass", "sorter": "number"},
    {"field": "StartRT", "title": "Start Retention Time (FWHM)", "sorter": "number"},
    {"field": "EndRT", "title": "End Retention Time (FWHM)", "sorter": "number"},
    {"field": "FeatureGroupQuantity", "title": "Feature Group Quantity",
     "sorter": "number"},
    {"field": "MinCharge", "title": "Min Charge", "sorter": "number"},
    {"field": "MaxCharge", "title": "Max Charge", "sorter": "number"},
    {"field": "MostAbundantFeatureCharge", "title": "Most Abundant Charge",
     "sorter": "number"},
    {"field": "IsotopeCosineScore", "title": "Isotope Cosine Score", "sorter": "number"},
]


def _sequence_view(file_manager, dataset_id, tool, cid, cache, p, settings):
    """Build the SequenceView wired for the tool (deconv global vs tnt per-proteoform).

    deconv: a single global sequence (``seq_deconv``) filtered by scan; peaks are
    the deconv-spectrum long frame (neutral masses -> ``deconvolved=True``).
    tnt: per-proteoform (``seq_tnt``) filtered by protein, with coverage +
    proteoform terminal columns; ``annotation_config`` (ion types / tolerance)
    is read from the oracle ``settings`` cache when available.
    """
    if tool == "flashtnt":
        anno_cfg = None
        if settings:
            anno_cfg = {
                "ion_types": settings.get("ion_types", ["b", "y"]),
                "tolerance": settings.get("tolerance", 20.0),
            }
        return SequenceView(
            cache_id=cid("sequence_view"),
            sequence_data_path=p("seq_tnt"),
            peaks_data_path=p("deconv_spectrum_tidy"),
            cache_path=cache,
            # protein selects the proteoform's sequence (seq_tnt has protein_id);
            # scan selects that proteoform's deconv peaks (deconv_spectrum_tidy has
            # scan_id, not protein_id) -- each filter applies only where its column
            # exists, reproducing the oracle's proteoform -> scan peak resolution.
            filters={"protein": "protein_id", "scan": "scan_id"},
            interactivity={"mass": "mass_in_scan"},
            # residue clicks publish the 0-based residue index as "aa" so the
            # augmented (tagger) spectrum can derive the tag-relative selectedAA.
            residue_identifier="aa",
            deconvolved=True,
            coverage_column="coverage",
            proteoform_start_column="proteoform_start",
            proteoform_end_column="proteoform_end",
            annotation_config=anno_cfg,
            title="Sequence View",
        )
    # flashdeconv: single global sequence
    return SequenceView(
        cache_id=cid("sequence_view"),
        sequence_data_path=p("seq_deconv"),
        peaks_data_path=p("deconv_spectrum_tidy"),
        cache_path=cache,
        filters={"scan": "scan_id"},
        interactivity={"mass": "mass_in_scan"},
        deconvolved=True,
        title="Sequence View",
    )


def make_builders(file_manager, dataset_id, tool, settings=None,
                  best_per_spectrum=True):
    """Return ``{comp_name: () -> BaseComponent}`` for one ``(tool, dataset)``.

    Args:
        file_manager: FLASHApp FileManager (provides ``result_path`` + ``cache_path``).
        dataset_id: the experiment id whose tidy caches were built by
            ``build_insight_caches``.
        tool: ``"flashdeconv"`` | ``"flashtnt"`` | ``"flashquant"`` (used for the
            sequence-view wiring and cache namespacing).
        settings: optional oracle ``settings`` dict (ion types / tolerance) for the
            FLASHTnT SequenceView.
        best_per_spectrum: round-8 finding 3-tables-002. When True (the oracle
            ProteinTable default), the ``protein_table`` builder shows only the
            single best-``Score`` proteoform per ``Scan`` (sourcing the
            ``is_best_per_scan == 1`` subset under a DISTINCT cache_id so toggling
            reliably swaps the cached row set); when False it shows all proteoforms.
            The FLASHTnT viewer wires this to a per-experiment "Best per spectrum"
            checkbox above its grid.

    Returns:
        A dict mapping every supported ``comp_name`` to a zero-arg factory. The
        grid lazily calls only the factories its layout references, so building
        this full dict is cheap (no Insight component is constructed here).
    """
    p = lambda tag: file_manager.result_path(dataset_id, tag)  # parquet path
    # Plot3D does not forward its x/y/z column config through the data_path=
    # subprocess (upstream limitation), so feed it the same on-disk tidy parquet
    # via data=scan_parquet(path) (in-process). These frames are per-scan /
    # per-feature small, so the memory tradeoff is negligible.
    scan = lambda tag: pl.scan_parquet(file_manager.result_path(dataset_id, tag))
    cid = lambda name: f"{tool}__{dataset_id}__{name}"
    cache = _insight_cache_dir(file_manager)

    # Tagger tag-payload resolution is only meaningful when a tags frame exists
    # (FLASHTnT). In FLASHDeconv the augmented spectrum has no tag overlay, so the
    # resolve kwargs are omitted (the tag selection simply never fires).
    tagger_tag_kwargs = (
        dict(
            tag_data_path=p("tags"), tag_id_column="tag_id",
            tag_sequence_column="TagSequence", tag_masses_column="mzs",
            tag_start_column="StartPos", selected_aa_identifier="aa",
        )
        if file_manager.result_exists(dataset_id, "tags")
        else {}
    )

    B = {
        # ---- FLASHDeconv / shared panels ----
        "scan_table": lambda: Table(
            cache_id=cid("scan_table"), data_path=p("scans"), cache_path=cache,
            interactivity={"scan": "scan_id"}, index_field="scan_id",
            default_row=0, title="Scan Table",
            # oracle Tabulator chrome: curated titles + guarded toFixed on RT /
            # PrecursorMass; shows ONLY these columns (no initialSort in the oracle).
            column_definitions=_SCAN_COLUMN_DEFS,
            # round-8 finding 3-tables-003: oracle TabulatorScanTable.vue
            # go-to-fields ['id','Scan'] -> schema id column is scan_id. Passing it
            # explicitly stops Table auto-detect from exposing the internal
            # mass_in_scan carrier as a go-to field.
            go_to_fields=["scan_id", "Scan"],
        ),
        "mass_table": lambda: Table(
            cache_id=cid("mass_table"), data_path=p("masses"), cache_path=cache,
            # mass selection == per-scan ordinal (the oracle massIndex), which the
            # 3D S/N plot consumes as SignalPeaks[mass_in_scan]; index_field stays
            # the global mass_id for row identity / go-to navigation.
            filters={"scan": "scan_id"}, interactivity={"mass": "mass_in_scan"},
            index_field="mass_id", title="Mass Table",
            # oracle chrome: toFixed on MonoMass/SumIntensity/CosineScore/SNR/QScore;
            # mass_in_scan stays in the data (interactivity) but is not displayed.
            column_definitions=_MASS_COLUMN_DEFS,
            # round-8 finding 3-tables-003: oracle TabulatorMassTable.vue
            # go-to-fields ['id'] -> schema id column is mass_id. Explicit list keeps
            # auto-detect from exposing the mass_in_scan / scan_id carriers.
            go_to_fields=["mass_id"],
        ),
        "deconv_spectrum": lambda: LinePlot(
            cache_id=cid("deconv_spectrum"), data_path=p("deconv_spectrum_tidy"),
            cache_path=cache, filters={"scan": "scan_id"},
            # clicking a deconvolved peak selects its mass (oracle onPlotClick
            # matched x against MonoMass and emitted the per-scan index).
            interactivity={"mass": "mass_in_scan"},
            x_column="mass", y_column="SumIntensity",
            # round-8 finding (deconv selective highlight): when a mass is selected
            # ("mass"), highlight the SELECTED mass's stick. The deconv base frame
            # carries one deconvolved mass per row (mass_in_scan), so the
            # match-column path lights up base rows where mass_in_scan == sel
            # directly (no link frame). No z=N charge labels and no
            # deconv_peaks_toggle for the deconvolved spectrum (oracle parity).
            # NOTE: the match-column highlight path
            # (lineplot._compute_selective_highlight) returns no charge/value
            # descriptors, so it draws NO selected-mass MonoMass value label. The
            # priority per the finding is the selected-stick highlight, which this
            # delivers; surfacing the MonoMass value as a label would require a new
            # match-column label producer in the LinePlot (not available today).
            highlight_selection="mass",
            highlight_match_column="mass_in_scan",
            # oracle axis titles (PlotlyLineplot.vue): deconvolved x="Monoisotopic
            # Mass", y="Intensity". Without these the axes show the raw column names.
            x_label="Monoisotopic Mass", y_label="Intensity",
            title="Deconvolved Spectrum",
        ),
        "anno_spectrum": lambda: LinePlot(
            cache_id=cid("anno_spectrum"), data_path=p("anno_spectrum_tidy"),
            cache_path=cache, filters={"scan": "scan_id"},
            # Clicking a raw m/z peak must NOT drive the shared "mass" slot (the
            # oracle onPlotClick matched against the deconvolved MonoMass array -- a
            # raw m/z never matches, so a click selected nothing; driving the shared
            # mass slot from here was a parity bug). BUT the selective-highlight LINK
            # path keys its highlight set off the FIRST interactivity column
            # (lineplot._compute_link_highlight / _attach_selective_highlight read
            # ``list(interactivity.values())[0]`` as the base ``id_column``), so the
            # annotated peaks MUST carry ``peak_id`` as their interactivity/index key
            # for the highlight-link key-set to map onto the drawn peaks. We publish
            # the click to a PRIVATE "anno_peak" slot (NOT consumed by any other
            # panel), keeping the parity-bug fix while exposing peak_id as id_column.
            interactivity={"anno_peak": "peak_id"},
            x_column="mz", y_column="intensity",
            # round-8 finding 3-anno-001: SELECTION-driven highlight. Drop the static
            # is_signal highlight; instead, when a deconvolved mass is selected
            # ("mass"), highlight that mass's SIGNAL peaks via the highlight LINK
            # frame (anno_highlight_link, 1:many peak->mass), with per-peak z=N charge
            # labels and the "Show Deconvolved Peaks" modebar toggle (oracle parity).
            highlight_selection="mass",
            highlight_link_path=p("anno_highlight_link"),
            highlight_link_key_column="peak_id",
            highlight_link_match_column="mass_in_scan",
            highlight_charge_column="charge",
            highlight_annotation_template="z={}",
            deconv_peaks_toggle=True,
            # oracle annotated-spectrum axis titles: x="m/z", y="Intensity".
            x_label="m/z", y_label="Intensity",
            title="Annotated Spectrum",
        ),
        "combined_spectrum": lambda: LinePlot.tagger(
            cache_id=cid("combined_spectrum"), data_path=p("combined_tagger"),
            cache_path=cache, filters={"scan": "scan_id"},
            interactivity={"tagger_mass": "peak_id"},
            x_column="MonoMass", y_column="SumIntensity",
            signal_peaks_column="SignalPeaks", mz_column="Mzs",
            mz_intensity_column="MzIntensities", tag_identifier="tag",
            # The tag table emits a scalar tag_id; resolve it to the tag's fragment
            # masses + sequence via the tags frame (mzs is a comma-string). A residue
            # click in the SequenceView sets "aa" -> tag-relative selectedAA (gold),
            # the value-based form of the oracle selectedAApos - startPos. Only wired
            # for FLASHTnT (where a tags frame exists); see tagger_tag_kwargs above.
            **tagger_tag_kwargs,
            title="Augmented Deconvolved Spectrum",
        ),
        "3D_SN_plot": lambda: Plot3D(
            cache_id=cid("3D_SN_plot"), data=scan("precursor_signals"),
            cache_path=cache,
            # Both scan AND mass are REQUIRED filters (no default for mass): the 3D
            # is empty until a mass is selected, matching the oracle. update.py
            # filters per_scan_data to the one selected scan, so the oracle frontend
            # getPrecursorSignal's precursor-scan lookup always fails when no mass is
            # set -> empty; only SignalPeaks[mass_index] is drawn once a mass is
            # chosen. (Do NOT make mass optional -- that would show all the scan's
            # peaks, which the oracle never did.)
            filters={"scan": "scan_id", "mass": "mass_in_scan"},
            filter_defaults={"scan": -1},
            # x-axis is the oracle "Mass" = mz*charge (precomputed in schema), NOT
            # raw m/z; Plot3D's default x_label "Mass" matches the oracle axis title.
            x_column="mass", y_column="charge", z_column="intensity",
            category_column="series",
            category_colors={"Signal": "#3366CC", "Noise": "#DC3912"},
            # round-8 finding 3-3d-001: DYNAMIC title (oracle Plotly3Dplot.vue). The
            # keys are the fixed scan/mass roles; the values are the SAME selection
            # identifiers this plot's ``filters`` read ("scan" / "mass"). Plot3D
            # computes the title from the live selection: '' when no scan is set,
            # 'Precursor signals' once a scan is selected (no mass), 'Mass signals'
            # once a mass is selected. The static ``title`` is the no-title fallback.
            title_selection={"scan": "scan", "mass": "mass"},
            title="Precursor Signals",
        ),
        # ---- heatmaps: reuse the existing full-resolution oracle caches as-is ----
        # oracle PlotlyHeatmap axis titles: x="Retention Time", y="Monoisotopic Mass".
        "ms1_deconv_heat_map": lambda: Heatmap(
            cache_id=cid("ms1_deconv_heat_map"), data_path=p("ms1_deconv_heatmap"),
            cache_path=cache, x_column="rt", y_column="mass",
            intensity_column="intensity",
            x_label="Retention Time", y_label="Monoisotopic Mass",
            title="Deconvolved MS1 Heatmap",
        ),
        "ms2_deconv_heat_map": lambda: Heatmap(
            cache_id=cid("ms2_deconv_heat_map"), data_path=p("ms2_deconv_heatmap"),
            cache_path=cache, x_column="rt", y_column="mass",
            intensity_column="intensity",
            x_label="Retention Time", y_label="Monoisotopic Mass",
            title="Deconvolved MS2 Heatmap",
        ),
        "ms1_raw_heatmap": lambda: Heatmap(
            cache_id=cid("ms1_raw_heatmap"), data_path=p("ms1_raw_heatmap"),
            cache_path=cache, x_column="rt", y_column="mass",
            intensity_column="intensity",
            x_label="Retention Time", y_label="Monoisotopic Mass",
            title="Raw MS1 Heatmap",
        ),
        "ms2_raw_heatmap": lambda: Heatmap(
            cache_id=cid("ms2_raw_heatmap"), data_path=p("ms2_raw_heatmap"),
            cache_path=cache, x_column="rt", y_column="mass",
            intensity_column="intensity",
            x_label="Retention Time", y_label="Monoisotopic Mass",
            title="Raw MS2 Heatmap",
        ),
        "fdr_plot": lambda: LinePlot.density(
            cache_id=cid("fdr_plot"), data_path=p("qscore_density"),
            cache_path=cache, x_column="x", y_column="y", category_column="group",
            target_value="target", decoy_value="decoy",
            # round-8 findings 3-fdr-001/002: oracle title "FDR Plot" (FDR_plotly.vue
            # args.title) and explicit trace legend names "Target QScores" /
            # "Decoy QScores" (FDR_plotly.vue trace ``name``s). targetLabel/decoyLabel
            # flow through ``config`` -> _plot_config -> _get_component_args_density.
            title="FDR Plot",
            config={"targetLabel": "Target QScores", "decoyLabel": "Decoy QScores"},
        ),
        "id_fdr_plot": lambda: LinePlot.density(
            cache_id=cid("id_fdr_plot"), data_path=p("qscore_density_id"),
            cache_path=cache, x_column="x", y_column="y", category_column="group",
            target_value="target", decoy_value="decoy",
            # round-8 findings 3-fdr-001/002: same as fdr_plot (oracle FDR_plotly.vue).
            title="FDR Plot",
            config={"targetLabel": "Target QScores", "decoyLabel": "Decoy QScores"},
        ),
        # ---- FLASHTnT panels ----
        # round-8 finding 3-tables-002: the oracle ProteinTable defaults to showing
        # only the best-Score proteoform per Scan (``bestPerSpectrumOnly: true``),
        # with a toggle to show all. We reproduce that server-side: when
        # ``best_per_spectrum`` is True the builder sources the
        # ``is_best_per_scan == 1`` subset (minted in schema._build_proteins:
        # exactly one row per Scan, highest Score, ties -> first-seen, matching the
        # oracle ``>`` keep-first) under a DISTINCT cache_id ("..protein_table_best")
        # so flipping the viewer toggle reliably swaps the cached row set; when False
        # it sources the full table under the normal cache_id. column_definitions /
        # interactivity / index_field / initial_sort are IDENTICAL in both branches.
        # Downstream cross-links (tag table, sequence view, augmented spectrum) key
        # off ``scan`` -- both row sets carry scan_id, so they are unaffected.
        "protein_table": lambda: Table(
            cache_id=cid("protein_table_best") if best_per_spectrum
            else cid("protein_table"),
            data=(
                pl.scan_parquet(p("proteins")).filter(
                    pl.col("is_best_per_scan") == 1
                )
                if best_per_spectrum
                else None
            ),
            data_path=None if best_per_spectrum else p("proteins"),
            cache_path=cache,
            # a protein-row click resolves to its scan (value-based
            # proteoform_scan_map): it sets BOTH the protein and the scan
            # selection, so the augmented spectrum / sequence-view peaks / tag
            # table all follow the selected proteoform to its scan.
            interactivity={"protein": "protein_id", "scan": "scan_id"},
            index_field="protein_id", default_row=0, title="Protein Table",
            # oracle chrome: curated titles, -1->"-" on Mass/Q-Value, initialSort
            # by Score desc. protein_id/scan_id carriers stay for index/cross-link
            # but are not displayed (no "Index" column in the oracle protein table).
            column_definitions=_PROTEIN_COLUMN_DEFS,
            initial_sort=_PROTEIN_INITIAL_SORT,
            # round-8 finding 3-tables-003: oracle TabulatorProteinTable.vue
            # go-to-fields ['Scan','accession']. Explicit list keeps auto-detect from
            # exposing the protein_id / scan_id carriers as go-to fields.
            go_to_fields=["Scan", "accession"],
        ),
        "tag_table": lambda: Table(
            cache_id=cid("tag_table"), data_path=p("tags"), cache_path=cache,
            # tags are scan data: show every tag on the selected proteoform's scan
            # (oracle filtered by Scan), driven by the protein->scan selection.
            filters={"scan": "scan_id"}, interactivity={"tag": "tag_id"},
            # oracle secondary filter: when a sequence residue is clicked, narrow to
            # tags spanning it (StartPos <= aa <= EndPos); shows all when no residue
            # is selected. The "aa" selection is published by the SequenceView.
            interval_filters={"aa": ("StartPos", "EndPos")},
            index_field="tag_id", title="Tag Table",
            # oracle chrome: curated titles, -1->"-" on N mass / C mass, initialSort
            # by Score desc. tag_id / mzs carriers stay for index/payload resolution
            # but are not displayed; StartPos/EndPos ARE displayed AND drive the
            # residue interval_filter.
            column_definitions=_TAG_COLUMN_DEFS,
            initial_sort=_TAG_INITIAL_SORT,
            # round-8 finding 3-tables-003: oracle TabulatorTagTable.vue go-to-fields
            # ['Scan','StartPos','EndPos','TagSequence']. Explicit list keeps
            # auto-detect from exposing the tag_id / scan_id carriers as go-to fields.
            go_to_fields=["Scan", "StartPos", "EndPos", "TagSequence"],
        ),
        "sequence_view": lambda: _sequence_view(
            file_manager, dataset_id, tool, cid, cache, p, settings
        ),
        # ---- FLASHQuant panels ----
        "quant_visualization": lambda: Table(
            cache_id=cid("quant_features"), data_path=p("quant_features"),
            cache_path=cache, interactivity={"feature": "feature_id"},
            # round-8 finding 3-feat-001: oracle FLASHQuantView TabulatorTable
            # title="Feature groups" (was "Features").
            index_field="feature_id", default_row=0, title="Feature groups",
            # oracle FLASHQuantView featureGroupTableColumnDefinitions: curated
            # titles (Index/Monoisotopic Mass/.../Isotope Cosine Score), no
            # formatters, no initialSort. StartRetentionTime(FWHM)/EndRetentionTime
            # (FWHM) -> schema StartRT/EndRT.
            column_definitions=_QUANT_COLUMN_DEFS,
            # round-8 finding 3-tables-003: the oracle FLASHQuantView TabulatorTable
            # passes NO go-to-fields, so its go-to UI never rendered. Pass [] to
            # DISABLE go-to (vs None, which would auto-detect and expose feature_id
            # etc. as a go-to field the oracle never had).
            go_to_fields=[],
        ),
        "quant_traces_3d": lambda: Plot3D(
            cache_id=cid("quant_traces"), data=scan("quant_traces"),
            cache_path=cache, filters={"feature": "feature_id"},
            filter_defaults={"feature": -1},
            # oracle FLASHQuantView: x = m/z, y = retention time, z = intensity
            # (Plot3D's defaults are precursor-flavored "Mass"/"Charge", so pass
            # explicit labels for the quant recipe).
            x_column="mz", y_column="rt", z_column="intensity",
            x_label="m/z", y_label="retention time", z_label="intensity",
            category_column="charge",
            # oracle builds one trace per charge but BREAKS the polyline between
            # EVERY trace within that charge (it pushes a -1000 z sentinel
            # before/after each trace's points). round-8 finding 3-quant-005:
            # (charge, isotope) is NOT unique -- two traces of one feature can share
            # it -- so keying the break on "isotope" would merge those two traces
            # into one connected polyline. series_column="trace_in_feature" (a stable
            # per-feature running trace id minted in schema._build_quant) breaks the
            # line per ACTUAL trace, while the legend/color stay per-charge.
            series_column="trace_in_feature",
            # oracle legend label is `Charge: ${charge}` (name: `Charge: 2`).
            category_name_template="Charge: {}",
            # oracle FLASHQuantView draws ONE connected elution line per charge
            # (mode:lines), not per-point stems; category_column already splits the
            # charges into separate traces, so disable the precursor-style stems.
            stem=False,
            title="Feature group signals",
        ),
    }
    return B
