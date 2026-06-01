"""Long-format adapters for the OpenMS-Insight migration.

FLASHApp's deconvolution caches historically store spectra as *arrays per scan*
(one row per scan, with ``MonoMass``/``SumIntensity`` list columns) and the old
``flash_viewer_grid`` filtered them by **row index** (``iloc[scanIndex]``).

OpenMS-Insight components filter by **column value** and expect **long format**
(one row per peak) with explicit identifier columns. These helpers explode the
existing per-scan frames into that long format so they can be fed directly to
``LinePlot`` / ``Table`` / ``Scatter3D`` with
``filters={'scanIndex': 'index'}`` (and ``massIndex`` where a per-scan peak
index is needed).

The functions are intentionally pure (Polars in, Polars out) and additive: they
do not touch the existing parse pipeline, so the old and new rendering paths can
coexist during the phased rollout.
"""

from typing import Optional

import polars as pl


def explode_spectrum_long(
    per_scan: pl.LazyFrame,
    *,
    index_column: str = "index",
    mass_array_column: str = "MonoMass",
    intensity_array_column: str = "SumIntensity",
    mass_out: str = "mass",
    intensity_out: str = "intensity",
    mass_id_out: str = "mass_id",
    drop_nonpositive_intensity: bool = False,
) -> pl.LazyFrame:
    """Explode an arrays-per-scan spectrum frame into long format.

    Each input row holds one scan with parallel ``MonoMass`` / ``SumIntensity``
    list columns. The output has one row per peak with:

        - ``index``     : the originating scan's row index (filter target for
                          ``scanIndex``); carried through verbatim.
        - ``mass``      : the peak mass / m/z.
        - ``intensity`` : the peak intensity.
        - ``mass_id``   : 0-based peak index within the scan (filter target for
                          ``massIndex``), assigned BEFORE any intensity filter so
                          it matches the original array position.

    Args:
        per_scan: LazyFrame with one row per scan and list-valued mass/intensity
            columns plus an ``index`` column.
        index_column: Name of the per-scan row-index column. Default "index".
        mass_array_column: List column of masses. Default "MonoMass".
        intensity_array_column: List column of intensities. Default "SumIntensity".
        mass_out: Output mass column name. Default "mass".
        intensity_out: Output intensity column name. Default "intensity".
        mass_id_out: Output per-scan peak-index column name. Default "mass_id".
        drop_nonpositive_intensity: If True, drop peaks with null/<=0 intensity
            AFTER mass_id assignment (default False — keep every peak so positions
            stay aligned with the original arrays).

    Returns:
        Long-format LazyFrame, sorted by ``index`` then ``mass_id``.
    """
    lf = per_scan.select(
        [
            pl.col(index_column).alias("index"),
            pl.col(mass_array_column).alias(mass_out),
            pl.col(intensity_array_column).alias(intensity_out),
        ]
    ).explode([mass_out, intensity_out])

    # Assign per-scan peak index over the original (pre-filter) order so it
    # matches the position in the source arrays — this is what massIndex selects.
    lf = lf.with_columns(
        pl.int_range(pl.len()).over("index").alias(mass_id_out)
    )

    if drop_nonpositive_intensity:
        lf = lf.filter(
            pl.col(intensity_out).is_not_null() & (pl.col(intensity_out) > 0)
        )

    return lf.sort(["index", mass_id_out])


def explode_combined_spectrum_long(
    per_scan: pl.LazyFrame,
    *,
    index_column: str = "index",
    deconv_mass_column: str = "MonoMass",
    deconv_intensity_column: str = "SumIntensity",
    anno_mass_column: str = "MonoMass_Anno",
    anno_intensity_column: str = "SumIntensity_Anno",
) -> tuple[pl.LazyFrame, pl.LazyFrame]:
    """Explode a combined (deconv + annotated) spectrum into two long frames.

    The FLASHApp ``combined_spectrum`` cache holds, per scan, both the
    deconvolved sticks (``MonoMass``/``SumIntensity``) and the raw/annotated
    peaks (``MonoMass_Anno``/``SumIntensity_Anno``). The augmented "tagger"
    spectrum overlays the latter on the former.

    Returns a ``(deconv_long, anno_long)`` pair, each in the
    :func:`explode_spectrum_long` schema, suitable for a ``LinePlot`` primary
    series + overlay series respectively (both filtered by ``scanIndex`` →
    ``index``).
    """
    deconv_long = explode_spectrum_long(
        per_scan,
        index_column=index_column,
        mass_array_column=deconv_mass_column,
        intensity_array_column=deconv_intensity_column,
    )
    anno_long = explode_spectrum_long(
        per_scan,
        index_column=index_column,
        mass_array_column=anno_mass_column,
        intensity_array_column=anno_intensity_column,
    )
    return deconv_long, anno_long


def explode_signal_peaks_long(
    per_scan: pl.LazyFrame,
    *,
    index_column: str = "index",
    signal_column: str = "SignalPeaks",
    noisy_column: str = "NoisyPeaks",
    signal_value: str = "signal",
    noise_value: str = "noise",
) -> pl.LazyFrame:
    """Explode per-scan signal/noisy peak arrays into long format for Scatter3D.

    FLASHApp stores ``SignalPeaks`` / ``NoisyPeaks`` as, per scan, a list over
    masses of a list of peaks, where each peak is ``[?, mz, intensity, charge]``
    (the 3D plot uses index 1 = m/z, 2 = intensity, 3 = charge; x is rendered as
    ``mz * charge``). This explodes both into one row per peak with:

        - ``index``     : scan row index (filter target for ``scanIndex``)
        - ``mass_id``   : mass index within the scan (filter target for
                          ``massIndex`` — isolates one mass's peaks)
        - ``mz``        : x (mass = mz * charge)
        - ``charge``    : y
        - ``intensity`` : z
        - ``kind``      : "signal" or "noise"

    Args:
        per_scan: LazyFrame with one row per scan and nested-list
            ``SignalPeaks``/``NoisyPeaks`` columns plus an ``index`` column.
        index_column: Per-scan row-index column. Default "index".
        signal_column: Nested signal-peaks column. Default "SignalPeaks".
        noisy_column: Nested noisy-peaks column. Default "NoisyPeaks".
        signal_value: ``kind`` value for signal peaks. Default "signal".
        noise_value: ``kind`` value for noise peaks. Default "noise".

    Returns:
        Long-format LazyFrame with columns
        index, mass_id, mz, charge, intensity, kind — Scatter3D-ready.
    """

    empty_schema = {
        "index": pl.Int64,
        "mass_id": pl.Int64,
        "mz": pl.Float64,
        "charge": pl.Float64,
        "intensity": pl.Float64,
        "kind": pl.Utf8,
    }

    def _one_kind(column: str, kind: str) -> pl.LazyFrame:
        # Level 1: list over masses -> add mass_id; Level 2: list over peaks.
        lf = per_scan.select(
            [
                pl.col(index_column).alias("index"),
                pl.col(column).alias("_peaks_by_mass"),
            ]
        )
        # Explode masses, then assign mass_id within each scan.
        lf = lf.explode("_peaks_by_mass").with_columns(
            pl.int_range(pl.len()).over("index").alias("mass_id")
        )
        # Now _peaks_by_mass is a list of peaks (each a list of floats).
        lf = lf.rename({"_peaks_by_mass": "_peak"}).explode("_peak")

        # Defensive: if the peak column carries no type information (e.g. an
        # all-empty column from untyped input), `.list.get()` would raise. In
        # real data the pyarrow schema keeps the list<list<double>> type, so
        # this only guards degenerate inputs — emit an empty typed frame.
        peak_dtype = lf.collect_schema().get("_peak")
        if not isinstance(peak_dtype, (pl.List, pl.Array)):
            return pl.LazyFrame(schema=empty_schema)

        lf = lf.filter(pl.col("_peak").is_not_null())
        lf = lf.with_columns(
            [
                pl.col("_peak").list.get(1).cast(pl.Float64).alias("mz"),
                pl.col("_peak").list.get(2).cast(pl.Float64).alias("intensity"),
                pl.col("_peak").list.get(3).cast(pl.Float64).alias("charge"),
                pl.lit(kind).alias("kind"),
            ]
        ).select(["index", "mass_id", "mz", "charge", "intensity", "kind"])
        return lf

    signal_lf = _one_kind(signal_column, signal_value)
    noise_lf = _one_kind(noisy_column, noise_value)
    return pl.concat([signal_lf, noise_lf]).sort(["index", "mass_id", "kind"])


def density_series_long(
    target_density: pl.DataFrame,
    decoy_density: Optional[pl.DataFrame] = None,
    *,
    target_label: str = "Target",
    decoy_label: str = "Decoy",
    x_column: str = "x",
    y_column: str = "y",
) -> pl.DataFrame:
    """Stack precomputed target/decoy density curves into one long frame.

    FLASHApp's FDR plot precomputes ``density_target`` / ``density_decoy`` as
    separate ``{x, y}`` frames. The OpenMS-Insight ``DensityPlot`` instead takes
    a single long frame with a ``series`` column (then computes the KDE itself);
    but when curves are already computed, this helper produces the equivalent
    long ``{series, x, y}`` frame directly for a thin pass-through path.

    Empty/absent decoy frames are handled (no Decoy rows emitted).
    """
    frames = []
    if target_density is not None and target_density.height > 0:
        frames.append(
            target_density.select(
                [
                    pl.lit(target_label).alias("series"),
                    pl.col(x_column).alias("x"),
                    pl.col(y_column).alias("y"),
                ]
            )
        )
    if decoy_density is not None and decoy_density.height > 0:
        frames.append(
            decoy_density.select(
                [
                    pl.lit(decoy_label).alias("series"),
                    pl.col(x_column).alias("x"),
                    pl.col(y_column).alias("y"),
                ]
            )
        )
    if not frames:
        return pl.DataFrame(
            schema={"series": pl.Utf8, "x": pl.Float64, "y": pl.Float64}
        )
    return pl.concat(frames)
