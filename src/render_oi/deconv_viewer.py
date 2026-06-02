"""OpenMS-Insight rendering engine for FLASHDeconv (migration Phase 1).

This is the replacement for ``src/render/render.py``'s ``render_grid`` that draws
each cell of the ``[experiment][row][col]`` layout with an individual
OpenMS-Insight component instead of the monolithic ``flash_viewer_grid`` Vue
component. It is additive: the old engine stays in place so the rollout can be
toggled per workflow.

Design:
- ``build_component(comp_name, ...)`` is a factory mapping each
  ``FLASHDeconvLayoutManager.COMPONENT_NAMES`` entry to an OpenMS-Insight
  component, loading the existing ``.pq`` caches through the long-format
  adapters in :mod:`src.parse.long_format`.
- Each experiment panel gets its OWN :class:`StateManager` (distinct
  ``session_key``) so selections never leak across side-by-side panels.
- Cross-linking uses the identifier→column model:
    scanIndex  : scan-table click → spectra / mass-table / 3D / sequence filter
    massIndex  : mass-table click → 3D plot optional isolation
- OI component caches live under ``{workspace}/cache/oi_cache/{dataset_id}/``;
  they are rebuilt only when missing (``_ensure_*`` helpers).

The Streamlit cross-link wiring (one shared StateManager per panel, components
composed with ``st.columns``) is performed by :func:`render_experiment`.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

import polars as pl

logger = logging.getLogger(__name__)

# State identifiers (FLASHApp StateTracker key → OpenMS-Insight identifier)
SCAN = "scanIndex"
MASS = "massIndex"

# Heatmap cache names → (title, MS level descriptor)
_HEATMAP_SPEC = {
    "ms1_deconv_heat_map": ("Deconvolved MS1 Heatmap", "ms1_deconv_heatmap"),
    "ms2_deconv_heat_map": ("Deconvolved MS2 Heatmap", "ms2_deconv_heatmap"),
    "ms1_raw_heatmap": ("Raw MS1 Heatmap", "ms1_raw_heatmap"),
    "ms2_raw_heatmap": ("Raw MS2 Heatmap", "ms2_raw_heatmap"),
}


def _oi_cache_dir(file_manager, dataset_id: str) -> str:
    """Per-dataset cache directory for OpenMS-Insight component caches."""
    base = Path(file_manager.cache_path) / "oi_cache" / dataset_id
    base.mkdir(parents=True, exist_ok=True)
    return str(base)


def _load_polars(file_manager, dataset_id: str, name: str) -> pl.LazyFrame:
    """Load a cached parquet result as a Polars LazyFrame."""
    res = file_manager.get_results(dataset_id, [name], use_polars=True)
    data = res[name]
    return data if isinstance(data, pl.LazyFrame) else pl.LazyFrame(data)


def _load_pandas(file_manager, dataset_id: str, name: str):
    """Load a cached parquet result as a pandas DataFrame (via Path)."""
    res = file_manager.get_results(dataset_id, [name])
    p = res[name]
    import pandas as pd

    return pd.read_parquet(p) if not isinstance(p, pd.DataFrame) else p


# --------------------------------------------------------------------------
# Component builders. Each returns a zero-arg callable that, when invoked,
# renders the component with the given StateManager + Streamlit key.
# --------------------------------------------------------------------------
def build_component(
    comp_name: str,
    dataset_id: str,
    file_manager,
    state_manager,
    key_prefix: str,
    has_sequence: bool = False,
) -> Optional[Callable[[], Any]]:
    """Build a render callable for one layout cell.

    Args:
        comp_name: A value from ``FLASHDeconvLayoutManager.COMPONENT_NAMES``.
        dataset_id: The selected experiment id.
        file_manager: FLASHApp FileManager for cache access.
        state_manager: The per-experiment OpenMS-Insight StateManager.
        key_prefix: Unique Streamlit key prefix for this panel (keeps
            side-by-side panels' component keys disjoint).
        has_sequence: Whether a sequence was submitted (enables sequence_view /
            internal_fragment_map).

    Returns:
        A zero-arg callable that renders the component, or None if the
        component name is unknown / unavailable.
    """
    from openms_insight import (
        DensityPlot,
        Heatmap,
        InternalFragmentMap,
        LinePlot,
        Scatter3D,
        SequenceView,
        Table,
    )
    from src.parse.long_format import (
        density_series_long,
        explode_combined_spectrum_long,
        explode_signal_peaks_long,
        explode_spectrum_long,
    )

    cache_dir = _oi_cache_dir(file_manager, dataset_id)
    cid = lambda name: f"{dataset_id}__{name}"  # noqa: E731
    skey = lambda name: f"{key_prefix}_{name}"  # noqa: E731

    # ---- Heatmaps ----
    if comp_name in _HEATMAP_SPEC:
        title, cache_name = _HEATMAP_SPEC[comp_name]
        data = _load_polars(file_manager, dataset_id, cache_name)
        is_deconv = "deconv" in comp_name
        hm = Heatmap(
            cache_id=cid(comp_name),
            data=data,
            x_column="rt",
            y_column="mass",
            intensity_column="intensity",
            title=title,
            x_label="Retention time",
            y_label="Monoisotopic mass" if is_deconv else "m/z",
            # Click a point -> scanIndex (all heatmaps) plus massIndex (deconv only),
            # restoring the legacy heatmap cross-links into the spectra/mass/3D panels.
            interactivity=(
                {SCAN: "scan_idx", MASS: "mass_idx"}
                if is_deconv
                else {SCAN: "scan_idx"}
            ),
            zoom_identifier=f"{comp_name}_zoom",
            cache_path=cache_dir,
        )
        return lambda: hm(key=skey(comp_name), state_manager=state_manager)

    # ---- Scan table (master; click sets scanIndex) ----
    if comp_name == "scan_table":
        data = _load_polars(file_manager, dataset_id, "scan_table")
        tbl = Table(
            cache_id=cid("scan_table"),
            data=data,
            interactivity={SCAN: "index"},
            index_field="index",
            title="Scan Table",
            # Legacy flash_viewer_grid column titles / tooltips / numeric
            # formatting (recovered from the built bundle's TabulatorScanTable
            # columnDefinitions). Built only from the columns actually present.
            column_definitions=_scan_table_column_definitions(data),
            cache_path=cache_dir,
        )
        return lambda: tbl(key=skey("scan_table"), state_manager=state_manager)

    # ---- Mass table (filtered by scanIndex; click sets massIndex) ----
    if comp_name == "mass_table":
        per_scan = _load_polars(file_manager, dataset_id, "mass_table")
        # Explode arrays-per-scan into one row per mass with mass_id.
        long = _explode_mass_table(per_scan)
        tbl = Table(
            cache_id=cid("mass_table"),
            data=long,
            filters={SCAN: "index"},
            interactivity={MASS: "mass_id"},
            index_field="mass_id",
            title="Mass Table",
            # Legacy flash_viewer_grid TabulatorMassTable column titles /
            # tooltips / numeric formatting (recovered from the built bundle),
            # built only from the exploded columns actually present.
            column_definitions=_mass_table_column_definitions(long),
            cache_path=cache_dir,
        )
        return lambda: tbl(key=skey("mass_table"), state_manager=state_manager)

    # ---- Deconvolved spectrum (LinePlot, filtered by scanIndex) ----
    if comp_name == "deconv_spectrum":
        # Build from combined_spectrum when available: it carries the SAME
        # deconvolved sticks (MonoMass/SumIntensity) PLUS the per-mass SignalPeaks
        # ([binIdx, mz, intensity, charge]), letting us derive a per-peak charge
        # label (z=N) locally — without touching long_format.py. Fall back to the
        # plain deconv_spectrum cache (no charge data) when combined is absent.
        long, ann_col = _deconv_spectrum_with_charge(
            file_manager, dataset_id, explode_spectrum_long
        )
        lp = LinePlot(
            cache_id=cid("deconv_spectrum"),
            data=long,
            filters={SCAN: "index"},
            x_column="mass",
            y_column="intensity",
            # Per-peak charge label on each deconvolved stick (legacy showed the
            # charge state next to peaks); None when no charge data is available.
            annotation_column=ann_col,
            title="Deconvolved Spectrum",
            x_label="Monoisotopic Mass",
            y_label="Intensity",
            cache_path=cache_dir,
        )
        return lambda: lp(key=skey("deconv_spectrum"), state_manager=state_manager)

    # ---- Annotated/raw spectrum (LinePlot over annotated peaks) ----
    if comp_name == "anno_spectrum":
        per_scan = _load_polars(file_manager, dataset_id, "combined_spectrum")
        _deconv_long, anno_long = explode_combined_spectrum_long(per_scan)
        # No per-peak charge label here: the annotated series (MonoMass_Anno /
        # SumIntensity_Anno) carries no charge information in these caches (only
        # the deconvolved series has the SignalPeaks charge constituents), so we
        # have no charge to annotate — left unannotated, matching the data.
        lp = LinePlot(
            cache_id=cid("anno_spectrum"),
            data=anno_long,
            filters={SCAN: "index"},
            x_column="mass",
            y_column="intensity",
            title="Annotated Spectrum",
            x_label="m/z",
            y_label="Intensity",
            cache_path=cache_dir,
        )
        return lambda: lp(key=skey("anno_spectrum"), state_manager=state_manager)

    # ---- 3D S/N plot (Scatter3D; scanIndex required, massIndex optional) ----
    if comp_name == "3D_SN_plot":
        per_scan = _load_polars(file_manager, dataset_id, "threedim_SN_plot")
        # x-axis is the deconvoluted mass (mz * charge), matching the legacy 3D plot
        # (it plotted peak[1]*peak[3]); the long format keeps mz and charge separate.
        long = explode_signal_peaks_long(per_scan).with_columns(
            (pl.col("mz") * pl.col("charge")).alias("mass")
        )
        s3 = Scatter3D(
            cache_id=cid("3D_SN_plot"),
            data=long,
            filters={SCAN: "index"},
            optional_filters={MASS: "mass_id"},
            mz_column="mass",
            title="Precursor Signals",
            cache_path=cache_dir,
        )

        def _render_3d():
            # Reflect the selection state in the title: the default view shows the
            # scan's full precursor S/N peaks ("Precursor Signals"); once a single
            # mass is isolated (massIndex set via the mass-table / deconv-heatmap
            # click) it shows that mass's signal/noisy peaks ("Mass Signals").
            # Only the displayed title arg changes — the cached/filtered data is
            # untouched — so this does not invalidate the component cache.
            mass_selected = (
                state_manager is not None
                and state_manager.get_selection(MASS) is not None
            )
            s3._title = "Mass Signals" if mass_selected else "Precursor Signals"
            return s3(key=skey("3D_SN_plot"), state_manager=state_manager)

        return _render_3d

    # ---- FDR / score-distribution plot (DensityPlot, precomputed curves) ----
    if comp_name == "fdr_plot":
        target = _load_pandas_pl(file_manager, dataset_id, "density_target")
        decoy = _load_pandas_pl(file_manager, dataset_id, "density_decoy")
        long = density_series_long(target, decoy)
        dp = DensityPlot(
            cache_id=cid("fdr_plot"),
            data=long.lazy(),
            precomputed=True,
            series_column="series",
            series_config={
                "Target": {"label": "Target QScores", "color": "green"},
                "Decoy": {"label": "Decoy QScores", "color": "red"},
            },
            title="Score Distribution",
            x_label="QScore",
            cache_path=cache_dir,
        )
        return lambda: dp(key=skey("fdr_plot"), state_manager=state_manager)

    # ---- Sequence view (only when a sequence is provided) ----
    if comp_name == "sequence_view" and has_sequence:
        builder = _build_sequence_view(
            dataset_id, file_manager, cache_dir, cid, skey, state_manager
        )
        if builder is not None:
            return builder

    # ---- Internal fragment map (only when a sequence is provided) ----
    if comp_name == "internal_fragment_map" and has_sequence:
        builder = _build_internal_fragment_map(
            dataset_id, file_manager, cache_dir, cid, skey, state_manager
        )
        if builder is not None:
            return builder

    logger.warning("Unknown / unavailable component: %s", comp_name)
    return None


def _load_pandas_pl(file_manager, dataset_id: str, name: str):
    """Load a parquet cache as a Polars DataFrame (eager) for density curves."""
    import pandas as pd

    res = file_manager.get_results(dataset_id, [name])
    p = res[name]
    pdf = pd.read_parquet(p) if not isinstance(p, pd.DataFrame) else p
    return pl.from_pandas(pdf)


def _explode_mass_table(per_scan: pl.LazyFrame) -> pl.LazyFrame:
    """Explode the arrays-per-scan mass_table into one row per mass.

    Columns: index (scan), mass_id, MonoMass, SumIntensity, charges/isotopes/
    scores — all the per-mass scalars the original Mass Table showed.
    """
    list_cols = [
        "MonoMass",
        "SumIntensity",
        "MinCharges",
        "MaxCharges",
        "MinIsotopes",
        "MaxIsotopes",
        "CosineScore",
        "SNR",
        "QScore",
    ]
    schema = per_scan.collect_schema().names()
    present = [c for c in list_cols if c in schema]
    lf = per_scan.select(["index", *present]).explode(present)
    lf = lf.with_columns(pl.int_range(pl.len()).over("index").alias("mass_id"))
    return lf.sort(["index", "mass_id"])


# --------------------------------------------------------------------------
# Tabulator column definitions, recovered from the legacy flash_viewer_grid
# bundle (TabulatorScanTable / TabulatorMassTable columnDefinitions). The
# legacy numeric formatter was ``v => v.toString().length > 4 ? v.toFixed(4)
# : v`` (a JS function that cannot be JSON-serialized into the OI cache); the
# closest portable equivalent is Tabulator's built-in ``money`` formatter with
# a fixed precision (no symbol), which OpenMS-Insight passes straight through.
# --------------------------------------------------------------------------
_FLOAT_FMT = {"formatter": "money", "formatterParams": {"precision": 4, "symbol": ""}}

# field -> (title, headerTooltip, is_float) for the scan table. Order follows
# the legacy column order; ``index`` is mapped to the displayed "Index" column.
_SCAN_TABLE_COLUMNS = [
    ("index", "Index", "The sequential index of the spectrum in the dataset.", False),
    ("Scan", "Scan Number", "The identifier of the mass spectrometry scan.", False),
    (
        "MSLevel",
        "MS Level",
        "The level of mass spectrometry analysis (e.g., MS1 or MS2).",
        False,
    ),
    (
        "RT",
        "Retention time",
        "The time at which the spectrum was detected during the chromatographic "
        "separation in seconds.",
        True,
    ),
    (
        "PrecursorMass",
        "Precursor Mass",
        "The mass of the precursor ion selected for fragmentation in Daltons.",
        True,
    ),
    ("#Masses", "#Masses", "The number of detected masses in the spectrum.", False),
]

# field -> (title, headerTooltip, is_float) for the exploded mass table. The
# ``mass_id`` index column is shown as "Index" (legacy "Index" column).
_MASS_TABLE_COLUMNS = [
    (
        "mass_id",
        "Index",
        "The sequential index of the mass entry in the dataset.",
        False,
    ),
    (
        "MonoMass",
        "Monoisotopic mass",
        "The monoisotopic mass of the detected ion in Daltons.",
        True,
    ),
    (
        "SumIntensity",
        "Sum intensity",
        "The total intensity of the detected mass across all isotopic peaks and "
        "charges.",
        True,
    ),
    (
        "MinCharges",
        "Min charge",
        "The minimum charge state detected for the mass.",
        False,
    ),
    (
        "MaxCharges",
        "Max charge",
        "The maximum charge state detected for the mass.",
        False,
    ),
    (
        "MinIsotopes",
        "Min isotope",
        "The smallest observed isotopic shift, expressed as a multiple of the "
        "average isotopic mass difference at 55kDA.",
        False,
    ),
    (
        "MaxIsotopes",
        "Max isotope",
        "The largest observed isotopic shift, expressed as a multiple of the "
        "average isotopic mass difference at 55kDA.",
        False,
    ),
    (
        "CosineScore",
        "Cosine score",
        "The cosine similarity score comparing the observed and theoretical "
        "isotopic patterns.",
        True,
    ),
    ("SNR", "SNR", "The signal-to-noise ratio for the detected mass.", True),
    (
        "QScore",
        "QScore",
        "The quality score indicating the confidence of the mass detection "
        "(higher is better).",
        True,
    ),
]


def _column_definitions(
    present_fields, spec
) -> List[Dict[str, Any]]:
    """Build Tabulator column_definitions from a (field,title,tooltip,float) spec.

    Only fields that are actually present in the data are emitted (so no column
    the data lacks is referenced, and — combined with always covering every real
    column — no existing column is dropped). Numeric columns get a ``number``
    sorter; float columns additionally get the fixed-precision ``money``
    formatter that stands in for the legacy ``toFixed(4)`` display.
    """
    present = set(present_fields)
    defs: List[Dict[str, Any]] = []
    for field, title, tooltip, is_float in spec:
        if field not in present:
            continue
        col: Dict[str, Any] = {
            "title": title,
            "field": field,
            "headerTooltip": tooltip,
            "sorter": "number",
        }
        if is_float:
            col.update(_FLOAT_FMT)
        defs.append(col)
    return defs


def _scan_table_column_definitions(data: pl.LazyFrame) -> List[Dict[str, Any]]:
    """Legacy scan-table column titles/tooltips/formatters for the real fields."""
    fields = data.collect_schema().names()
    return _column_definitions(fields, _SCAN_TABLE_COLUMNS)


def _mass_table_column_definitions(long: pl.LazyFrame) -> List[Dict[str, Any]]:
    """Legacy mass-table column titles/tooltips/formatters for the exploded fields."""
    fields = long.collect_schema().names()
    return _column_definitions(fields, _MASS_TABLE_COLUMNS)


def _deconv_spectrum_with_charge(file_manager, dataset_id, explode_spectrum_long):
    """Deconvolved-spectrum long format + a per-peak ``charge_label`` (``z=N``).

    Returns ``(long_frame, annotation_column_or_None)``.

    The plain ``deconv_spectrum`` cache holds only ``MonoMass``/``SumIntensity``
    (no charge). ``combined_spectrum`` carries the SAME deconvolved sticks PLUS
    each mass's constituent ``SignalPeaks`` (``[binIdx, mz, intensity, charge]``),
    so we derive a representative charge per deconvolved mass (the charge of its
    most intense constituent peak) and format it as ``z=N``. This is built
    locally via Polars expressions — ``long_format.py`` is untouched.

    Falls back to the plain ``deconv_spectrum`` (no annotation) when
    ``combined_spectrum`` / ``SignalPeaks`` are unavailable, in which case the
    annotation column is ``None``.
    """
    try:
        combined = _load_polars(file_manager, dataset_id, "combined_spectrum")
        schema = combined.collect_schema().names()
        if not {"MonoMass", "SumIntensity", "SignalPeaks"} <= set(schema):
            raise KeyError("combined_spectrum lacks SignalPeaks")
        long = (
            combined.select(["index", "MonoMass", "SumIntensity", "SignalPeaks"])
            .explode(["MonoMass", "SumIntensity", "SignalPeaks"])
            .rename({"MonoMass": "mass", "SumIntensity": "intensity"})
            # mass_id matches explode_spectrum_long's per-scan peak order so the
            # scan cross-link row counts stay identical.
            .with_columns(pl.int_range(pl.len()).over("index").alias("mass_id"))
        )
        # Representative charge = charge (peak field 3) of the constituent peak
        # with the maximum intensity (peak field 2). Null when a mass has no
        # constituent peaks; then no label is shown for that stick.
        long = long.with_columns(
            pl.col("SignalPeaks")
            .list.eval(pl.element().list.get(3))
            .list.get(
                pl.col("SignalPeaks")
                .list.eval(pl.element().list.get(2))
                .list.arg_max()
            )
            .cast(pl.Int64, strict=False)
            .alias("_charge")
        ).with_columns(
            pl.when(pl.col("_charge").is_not_null())
            .then(pl.format("z={}", pl.col("_charge")))
            .otherwise(pl.lit(""))
            .alias("charge_label")
        )
        long = long.select(
            ["index", "mass_id", "mass", "intensity", "charge_label"]
        ).sort(["index", "mass_id"])
        return long, "charge_label"
    except Exception:  # pragma: no cover - defensive fallback to plain cache
        logger.info(
            "deconv_spectrum charge labels unavailable (no SignalPeaks); "
            "rendering without charge annotation"
        )
        per_scan = _load_polars(file_manager, dataset_id, "deconv_spectrum")
        return explode_spectrum_long(per_scan), None


def _build_sequence_view(
    dataset_id, file_manager, cache_dir, cid, skey, state_manager
):
    """Build a SequenceView for FLASHDeconv from the submitted sequence.

    FLASHDeconv stores a single submitted sequence (not per-proteoform); the
    sequence view is filtered by scanIndex only to show the deconvolved peaks
    of the selected scan against that fixed sequence.
    """
    from openms_insight import SequenceView
    from src.parse.long_format import explode_spectrum_long

    if not file_manager.result_exists("sequence", "sequence"):
        return None
    seq = file_manager.get_results("sequence", "sequence")["sequence"]
    sequence_str = seq.get("input_sequence")
    if not sequence_str:
        return None

    # Deconvolved peaks (neutral masses) for matching, filtered by scan.
    per_scan = _load_polars(file_manager, dataset_id, "deconv_spectrum")
    peaks_long = (
        explode_spectrum_long(per_scan)
        .with_columns(pl.int_range(pl.len()).over("index").alias("peak_id"))
    )

    sv = SequenceView(
        cache_id=cid("sequence_view"),
        sequence_data=sequence_str,  # static sequence string
        peaks_data=peaks_long,
        filters={SCAN: "index"},
        deconvolved=True,
        fixed_modifications=_fixed_mods_from_sequence(seq),
        cache_path=cache_dir,
    )
    return lambda: sv(key=skey("sequence_view"), state_manager=state_manager)


def _build_internal_fragment_map(
    dataset_id, file_manager, cache_dir, cid, skey, state_manager
):
    """Build an InternalFragmentMap from the submitted sequence + scan peaks."""
    from openms_insight import InternalFragmentMap
    from src.parse.long_format import explode_spectrum_long

    if not file_manager.result_exists("sequence", "sequence"):
        return None
    seq = file_manager.get_results("sequence", "sequence")["sequence"]
    sequence_str = seq.get("input_sequence")
    if not sequence_str:
        return None

    per_scan = _load_polars(file_manager, dataset_id, "deconv_spectrum")
    peaks_long = explode_spectrum_long(per_scan)

    ifm = InternalFragmentMap(
        cache_id=cid("internal_fragment_map"),
        sequence_data=sequence_str,
        peaks_data=peaks_long,
        mass_column="mass",
        filters={SCAN: "index"},
        title="Internal Fragment Map",
        cache_path=cache_dir,
    )
    return lambda: ifm(key=skey("internal_fragment_map"), state_manager=state_manager)


def _fixed_mods_from_sequence(seq: Dict[str, Any]) -> List[str]:
    """Derive fixed-mod residue letters from the stored sequence settings."""
    mods = []
    if seq.get("fixed_mod_cysteine"):
        mods.append("C")
    if seq.get("fixed_mod_methionine"):
        mods.append("M")
    return mods


def render_experiment(
    dataset_id: str,
    layout_rows: List[List[str]],
    file_manager,
    panel_key: str,
    has_sequence: bool = False,
) -> None:
    """Render one experiment panel's [row][col] grid with OpenMS-Insight.

    Args:
        dataset_id: The selected experiment id.
        layout_rows: The experiment's layout — a list of rows, each a list of
            ``COMPONENT_NAMES`` strings (≤3 columns per row).
        file_manager: FLASHApp FileManager for cache access.
        panel_key: Unique key for this panel; also seeds the StateManager
            ``session_key`` so SIDE-BY-SIDE PANELS HAVE INDEPENDENT SELECTIONS
            (a distinct session_key per experiment prevents cross-contamination).
        has_sequence: Whether a sequence was submitted (enables sequence_view /
            internal_fragment_map).
    """
    import streamlit as st
    from openms_insight import StateManager

    # Per-experiment StateManager — distinct session_key keeps selections from
    # leaking across side-by-side panels (Risks/watch-items in the plan).
    state_manager = StateManager(session_key=f"oi_state_{panel_key}")

    for row_index, row in enumerate(layout_rows):
        if not row:
            continue
        cols = st.columns(len(row))
        for col_index, comp_name in enumerate(row):
            with cols[col_index]:
                try:
                    render = build_component(
                        comp_name,
                        dataset_id,
                        file_manager,
                        state_manager,
                        key_prefix=f"{panel_key}_{row_index}_{col_index}",
                        has_sequence=has_sequence,
                    )
                    if render is not None:
                        render()
                    else:
                        st.warning(f"Component unavailable: {comp_name}")
                except Exception as exc:  # pragma: no cover - defensive UI guard
                    logger.exception("Failed to render %s", comp_name)
                    st.error(f"Error rendering {comp_name}: {exc}")
