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


def make_builders(file_manager, dataset_id, tool, settings=None):
    """Return ``{comp_name: () -> BaseComponent}`` for one ``(tool, dataset)``.

    Args:
        file_manager: FLASHApp FileManager (provides ``result_path`` + ``cache_path``).
        dataset_id: the experiment id whose tidy caches were built by
            ``build_insight_caches``.
        tool: ``"flashdeconv"`` | ``"flashtnt"`` | ``"flashquant"`` (used for the
            sequence-view wiring and cache namespacing).
        settings: optional oracle ``settings`` dict (ion types / tolerance) for the
            FLASHTnT SequenceView.

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

    B = {
        # ---- FLASHDeconv / shared panels ----
        "scan_table": lambda: Table(
            cache_id=cid("scan_table"), data_path=p("scans"), cache_path=cache,
            interactivity={"scan": "scan_id"}, index_field="scan_id",
            default_row=0, title="Scan Table",
        ),
        "mass_table": lambda: Table(
            cache_id=cid("mass_table"), data_path=p("masses"), cache_path=cache,
            # mass selection == per-scan ordinal (the oracle massIndex), which the
            # 3D S/N plot consumes as SignalPeaks[mass_in_scan]; index_field stays
            # the global mass_id for row identity / go-to navigation.
            filters={"scan": "scan_id"}, interactivity={"mass": "mass_in_scan"},
            index_field="mass_id", title="Mass Table",
        ),
        "deconv_spectrum": lambda: LinePlot(
            cache_id=cid("deconv_spectrum"), data_path=p("deconv_spectrum_tidy"),
            cache_path=cache, filters={"scan": "scan_id"},
            # clicking a deconvolved peak selects its mass (oracle onPlotClick
            # matched x against MonoMass and emitted the per-scan index).
            interactivity={"mass": "mass_in_scan"},
            x_column="mass", y_column="SumIntensity",
            title="Deconvolved Spectrum",
        ),
        "anno_spectrum": lambda: LinePlot(
            cache_id=cid("anno_spectrum"), data_path=p("anno_spectrum_tidy"),
            cache_path=cache, filters={"scan": "scan_id"},
            # NO mass interactivity: the annotated (raw m/z) spectrum's x is m/z,
            # but the oracle onPlotClick matched the click against the deconvolved
            # MonoMass array -- a raw m/z never matches, so clicking it selected
            # nothing. (Driving the shared mass slot from here was a parity bug.)
            x_column="mz", y_column="intensity", highlight_column="is_signal",
            title="Annotated Spectrum",
        ),
        "combined_spectrum": lambda: LinePlot.tagger(
            cache_id=cid("combined_spectrum"), data_path=p("combined_tagger"),
            cache_path=cache, filters={"scan": "scan_id"},
            interactivity={"tagger_mass": "peak_id"},
            x_column="MonoMass", y_column="SumIntensity",
            signal_peaks_column="SignalPeaks", mz_column="Mzs",
            mz_intensity_column="MzIntensities", tag_identifier="tag",
            title="Augmented Deconvolved Spectrum",
        ),
        "3D_SN_plot": lambda: Plot3D(
            cache_id=cid("3D_SN_plot"), data=scan("precursor_signals"),
            cache_path=cache,
            filters={"scan": "scan_id", "mass": "mass_in_scan"},
            filter_defaults={"scan": -1},
            x_column="mz", y_column="charge", z_column="intensity",
            category_column="series",
            category_colors={"Signal": "#3366CC", "Noise": "#DC3912"},
            title="Precursor Signals",
        ),
        # ---- heatmaps: reuse the existing full-resolution oracle caches as-is ----
        "ms1_deconv_heat_map": lambda: Heatmap(
            cache_id=cid("ms1_deconv_heat_map"), data_path=p("ms1_deconv_heatmap"),
            cache_path=cache, x_column="rt", y_column="mass",
            intensity_column="intensity", title="Deconvolved MS1 Heatmap",
        ),
        "ms2_deconv_heat_map": lambda: Heatmap(
            cache_id=cid("ms2_deconv_heat_map"), data_path=p("ms2_deconv_heatmap"),
            cache_path=cache, x_column="rt", y_column="mass",
            intensity_column="intensity", title="Deconvolved MS2 Heatmap",
        ),
        "ms1_raw_heatmap": lambda: Heatmap(
            cache_id=cid("ms1_raw_heatmap"), data_path=p("ms1_raw_heatmap"),
            cache_path=cache, x_column="rt", y_column="mass",
            intensity_column="intensity", title="Raw MS1 Heatmap",
        ),
        "ms2_raw_heatmap": lambda: Heatmap(
            cache_id=cid("ms2_raw_heatmap"), data_path=p("ms2_raw_heatmap"),
            cache_path=cache, x_column="rt", y_column="mass",
            intensity_column="intensity", title="Raw MS2 Heatmap",
        ),
        "fdr_plot": lambda: LinePlot.density(
            cache_id=cid("fdr_plot"), data_path=p("qscore_density"),
            cache_path=cache, x_column="x", y_column="y", category_column="group",
            target_value="target", decoy_value="decoy",
            title="Score Distribution",
        ),
        "id_fdr_plot": lambda: LinePlot.density(
            cache_id=cid("id_fdr_plot"), data_path=p("qscore_density_id"),
            cache_path=cache, x_column="x", y_column="y", category_column="group",
            target_value="target", decoy_value="decoy",
            title="Score Distribution",
        ),
        # ---- FLASHTnT panels ----
        "protein_table": lambda: Table(
            cache_id=cid("protein_table"), data_path=p("proteins"),
            cache_path=cache,
            # a protein-row click resolves to its scan (value-based
            # proteoform_scan_map): it sets BOTH the protein and the scan
            # selection, so the augmented spectrum / sequence-view peaks / tag
            # table all follow the selected proteoform to its scan.
            interactivity={"protein": "protein_id", "scan": "scan_id"},
            index_field="protein_id", default_row=0, title="Protein Table",
        ),
        "tag_table": lambda: Table(
            cache_id=cid("tag_table"), data_path=p("tags"), cache_path=cache,
            # tags are scan data: show every tag on the selected proteoform's scan
            # (oracle filtered by Scan), driven by the protein->scan selection.
            filters={"scan": "scan_id"}, interactivity={"tag": "tag_id"},
            index_field="tag_id", title="Tag Table",
        ),
        "sequence_view": lambda: _sequence_view(
            file_manager, dataset_id, tool, cid, cache, p, settings
        ),
        # ---- FLASHQuant panels ----
        "quant_visualization": lambda: Table(
            cache_id=cid("quant_features"), data_path=p("quant_features"),
            cache_path=cache, interactivity={"feature": "feature_id"},
            index_field="feature_id", default_row=0, title="Features",
        ),
        "quant_traces_3d": lambda: Plot3D(
            cache_id=cid("quant_traces"), data=scan("quant_traces"),
            cache_path=cache, filters={"feature": "feature_id"},
            filter_defaults={"feature": -1},
            x_column="rt", y_column="mz", z_column="intensity",
            category_column="charge", title="Feature Traces",
        ),
    }
    return B
