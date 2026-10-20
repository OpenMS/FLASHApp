import pandas as pd
import polars as pl
import numpy as np

from src.parse.masstable import parseFLASHDeconvOutput, getMSSignalDF, getSpectraTableDF
from src.render.compression import downsample_heatmap, compute_compression_levels
from scipy.stats import gaussian_kde

# One row per scan with heavy array-valued cells; small row groups so per-scan
# pushdown reads only the matching group(s) instead of the whole file.
SPECTRA_ROW_GROUP_SIZE = 64

# Long-format (one-row-per-peak / one-row-per-mass) frames are consumed by the
# OpenMS-Insight components (LinePlot / Table / Scatter3D), which filter by COLUMN
# VALUE (filters={'scanIndex':'index'}) rather than by the old `iloc[scanIndex]`
# ROW-INDEX path in src/render/update.py. These producers are ADDITIVE: the legacy
# array-per-scan frames (deconv_spectrum, anno_spectrum, combined_spectrum,
# mass_table) are still stored unchanged for the old render path; the long frames
# are stored under separate `*_long` tags so both can coexist during Stage B.
#
# Exploded rows are sorted/grouped by `index` so a value-filter `index == k`
# reads exactly the rows the old `iloc[k]` array slice produced (the legacy frames
# are built with with_row_index + sort('index'), so row position == index value).
# Long peak frames use a slightly larger row group (predicate pushdown is by
# value, and per-scan peak counts are modest) than the array frames.
LONG_ROW_GROUP_SIZE = 10_000


def _explode_long_by_position(indexed_lf, id_col, value_exprs):
    """Explode parallel per-scan list columns into one row per position.

    This reproduces, exactly, the legacy FLASHApp Vue per-column expansion
    (``TabulatorMassTable.vue`` ``tableData`` and the lineplot/3D consumers):
    each per-scan list column is laid out by POSITION ``0..L-1`` independently;
    the number of rows for a scan is the MAXIMUM list length across the supplied
    columns, and any column shorter than that maximum yields ``null`` for the
    missing trailing positions (the JS ``undefined``). Columns are therefore
    ALIGNED BY POSITION, never lock-stepped — important because in real
    FLASHDeconv output ``mz_array``/``intensity_array`` (the full spectrum) can be
    LONGER than the per-mass ``MinCharges``/``SignalPeaks`` axis, and the legacy
    UI pads the short columns with blanks rather than truncating.

    Args:
        indexed_lf: a polars LazyFrame carrying an integer ``index`` column.
        id_col: name of the per-scan position id column to emit (``peak_id`` or
            ``mass_id``) — the 0-based row position within the scan.
        value_exprs: list of ``(out_name, list_expr)``. ``list_expr`` is a polars
            expression evaluating to a per-scan list; its element at position
            ``id_col`` becomes the scalar ``out_name`` (null when out of range).

    Returns:
        LazyFrame with columns ``index``, ``id_col``, then each ``out_name``,
        sorted by ``index`` then ``id_col``. Scans whose columns are all empty
        contribute 0 rows (matching the old ``iloc[k]`` empty-array slice).
    """
    out_names = [name for name, _ in value_exprs]

    # Per-scan max length across all contributing list columns → number of
    # positions to emit. (max of list.len() over the columns.)
    max_len = value_exprs[0][1].list.len()
    for _, expr in value_exprs[1:]:
        max_len = pl.max_horizontal(max_len, expr.list.len())

    # Pad every value list to the per-scan ``max_len`` with nulls — gathering by
    # the position range with ``null_on_oob`` reproduces the legacy blank-tail
    # cell — then zip-explode the id column and all value lists together so each
    # output row is exactly one position.
    #
    # This stays O(total output rows). The earlier approach exploded only the id
    # column while every row still carried the full per-scan value lists, then
    # gathered the scalar — i.e. O(rows × max_len), which materialises the lists
    # `max_len` times and OOMs on real spectra (e.g. an 865k-row annotated frame
    # with multi-thousand-length lists ≈ tens of GB). Exploding the lists in
    # lock-step avoids the duplication entirely.
    positions = pl.int_ranges(0, max_len)
    lf = (
        indexed_lf
        .select(
            [pl.col("index")]
            + [expr.list.gather(positions, null_on_oob=True).alias(name)
               for name, expr in value_exprs]
            + [positions.alias(id_col)]
        )
        .explode([id_col] + out_names)
        # Empty scans explode to a single null-id row; drop so they contribute 0 rows.
        .filter(pl.col(id_col).is_not_null())
        .sort(["index", id_col])
    )
    return lf.select(["index", id_col] + out_names)


# Each per-mass entry of the nested ``SignalPeaks`` column is a list of matched
# signal peaks, and each signal peak is a 4-tuple ``[peak_index, mz, intensity,
# charge]`` (verified against the example caches under example-data/workspaces/**:
# ``SignalPeaks`` has dtype ``List(List(List(Float64)))`` — scan → mass → peak →
# quadruple). The charge drill-down ("Augmented Annotated Spectrum") needs, for
# each deconvolved mass, the per-mass mz / charge / intensity of its signal peaks
# (not just the ``is_signal`` boolean). These expressions pull those three axes
# out of the nested column WITHOUT flattening across masses: each evaluates to a
# per-scan ``List(List(...))`` whose outer position is the mass axis (aligned to
# ``MonoMass`` / ``peak_id``) and whose inner list is that one mass's small signal
# arrays. _explode_long_by_position then gathers the inner list by mass position,
# so each output row carries exactly the signal peaks belonging to that peak.
def _signal_mzs_expr():
    # mz is element index 1 of each [idx, mz, intensity, charge] signal peak.
    return pl.col("SignalPeaks").list.eval(
        pl.element().list.eval(pl.element().list.get(1))
    )


def _signal_intensities_expr():
    # intensity is element index 2 of each signal-peak quadruple.
    return pl.col("SignalPeaks").list.eval(
        pl.element().list.eval(pl.element().list.get(2))
    )


def _signal_charges_expr():
    # charge is element index 3; stored as float in the nested array, cast to i64.
    return pl.col("SignalPeaks").list.eval(
        pl.element().list.eval(pl.element().list.get(3).cast(pl.Int64))
    )


def _signal_is_signal_expr():
    # Per-mass boolean: True where that mass has >=1 matched signal peak.
    return pl.col("SignalPeaks").list.eval(pl.element().list.len() > 0)


def _fill_empty_signal_lists(lf):
    """Replace ``null`` signal/charge/intensity cells with empty lists.

    Positions past the per-mass ``SignalPeaks`` axis (e.g. a ragged scan whose
    full-spectrum ``mz_array`` is longer than ``SignalPeaks``) explode to ``null``
    list cells. Coerce them to empty lists so consumers always read a list (parity
    with ``is_signal`` being filled to ``False``), and every signal-flagged peak
    has equal-length ``signal_mzs`` / ``signal_charges`` / ``signal_intensities``.
    """
    return lf.with_columns([
        pl.col("signal_mzs").fill_null([]).alias("signal_mzs"),
        pl.col("signal_charges").fill_null([]).alias("signal_charges"),
        pl.col("signal_intensities").fill_null([]).alias("signal_intensities"),
    ])


def deconv_spectrum_long(pl_deconv_indexed):
    """One row per deconvolved peak with the per-mass signal-peak arrays.

    Columns: index, peak_id, MonoMass, SumIntensity, signal_mzs (list[f64]),
    signal_charges (list[i64]), signal_intensities (list[f64]).

    Long-format replacement for the array-valued ``deconv_spectrum`` frame,
    consumed by ``LinePlot(filters={'scanIndex':'index'}, x_column='MonoMass',
    y_column='SumIntensity')``.

    The three ``signal_*`` list columns carry, for the deconvolved mass at this
    ``peak_id`` position, the mz / charge / intensity of each of its matched
    signal peaks (the per-mass ``SignalPeaks[peak_id]`` axis). They are aligned to
    each other (same length, one entry per signal peak of this mass) and to the
    ``peak_id``/``MonoMass`` row, and back the "Augmented Annotated Spectrum"
    charge drill-down. A peak with no matched signal (or a ragged past-end
    position) carries empty lists.
    """
    lf = _explode_long_by_position(
        pl_deconv_indexed,
        "peak_id",
        [("MonoMass", pl.col("mz_array")),
         ("SumIntensity", pl.col("intensity_array")),
         ("signal_mzs", _signal_mzs_expr()),
         ("signal_charges", _signal_charges_expr()),
         ("signal_intensities", _signal_intensities_expr())],
    )
    return _fill_empty_signal_lists(lf)


def anno_spectrum_long(pl_anno_indexed):
    """One row per annotated/raw peak: index, peak_id, MonoMass_Anno,
    SumIntensity_Anno.

    Long-format replacement for the array-valued ``anno_spectrum`` frame,
    consumed by ``LinePlot(filters={'scanIndex':'index'},
    x_column='MonoMass_Anno', y_column='SumIntensity_Anno')``.
    """
    return _explode_long_by_position(
        pl_anno_indexed,
        "peak_id",
        [("MonoMass_Anno", pl.col("mz_array")),
         ("SumIntensity_Anno", pl.col("intensity_array"))],
    )


def combined_spectrum_long(pl_deconv_indexed):
    """One row per deconvolved peak with a signal-membership flag and the
    per-mass signal-peak arrays.

    Columns: index, peak_id, MonoMass, SumIntensity, is_signal (bool),
    signal_mzs (list[f64]), signal_charges (list[i64]),
    signal_intensities (list[f64]).

    ``is_signal`` is True when the corresponding per-mass entry of the nested
    ``SignalPeaks`` column is non-empty, i.e. the deconvolved mass at that
    position has at least one matched signal peak (mirrors the per-mass alignment
    the 3D plot uses: ``Plotly3Dplot.vue`` indexes ``SignalPeaks[massIndex]`` by
    the same position). ``SignalPeaks`` is the per-mass axis and in real output
    can be SHORTER than ``mz_array``; positions beyond its length therefore have
    no signal entry and are flagged ``False`` (parity with the JS ``undefined``
    → no-signal).

    ``signal_mzs`` / ``signal_charges`` / ``signal_intensities`` carry, for the
    mass at this ``peak_id`` position, the mz / charge / intensity of each matched
    signal peak (the contents of ``SignalPeaks[peak_id]``). The three lists are
    mutually aligned (one entry per signal peak of this mass, equal length) and
    aligned to the ``peak_id``/``MonoMass`` row; they back the legacy "Augmented
    Annotated Spectrum" charge drill-down. When ``is_signal`` is False (no matched
    signal, or a ragged past-end position) all three lists are empty.

    This is the long-format counterpart of the array-valued ``combined_spectrum``
    deconv side; the annotated overlay is provided separately by
    ``anno_spectrum_long`` (the OpenMS-Insight LinePlot reads the 2nd series from
    its own ``x2_column``/``y2_column`` frame).
    """
    # Per-mass list columns, all aligned to the SignalPeaks (per-mass) axis.
    # _explode_long_by_position gathers each by the same position id as MonoMass
    # and yields null past its end (coerced below to False / empty lists).
    lf = _explode_long_by_position(
        pl_deconv_indexed,
        "peak_id",
        [("MonoMass", pl.col("mz_array")),
         ("SumIntensity", pl.col("intensity_array")),
         ("is_signal", _signal_is_signal_expr()),
         ("signal_mzs", _signal_mzs_expr()),
         ("signal_charges", _signal_charges_expr()),
         ("signal_intensities", _signal_intensities_expr())],
    )
    lf = lf.with_columns(pl.col("is_signal").fill_null(False))
    return _fill_empty_signal_lists(lf)


def mass_table_long(pl_deconv_indexed):
    """One row per mass: index, mass_id, plus scalar mass-table fields.

    Long-format replacement for the array-valued ``mass_table`` frame. Each row is
    one deconvolved mass within a scan; ``MonoMass``/``SumIntensity`` and the
    per-mass charge/isotope/score columns become scalars.

    Consumed by ``Table(interactivity={'massIndex':'mass_id'},
    filters={'scanIndex':'index'})``: clicking a row sets ``massIndex`` to the
    row's ``mass_id``, and the table is filtered to the selected scan via
    ``index``. ``mass_id`` is the 0-based position of the mass within its scan,
    matching the array-subscript semantics the 3D plot uses for ``massIndex``.

    Columns are aligned BY POSITION (not lock-stepped): the legacy
    ``TabulatorMassTable.vue`` builds one row per position up to the MAX array
    length across the required columns, leaving blanks where a column is shorter.
    In real FLASHDeconv output ``MonoMass``/``SumIntensity`` (the full spectrum
    ``mz_array``/``intensity_array``) may be LONGER than the per-mass charge/
    isotope/score arrays; those trailing rows therefore carry the mass/intensity
    with ``null`` charge/isotope/score cells, exactly as the old UI rendered them.
    """
    value_exprs = [
        ("MonoMass", pl.col("mz_array")),
        ("SumIntensity", pl.col("intensity_array")),
        ("MinCharges", pl.col("MinCharges")),
        ("MaxCharges", pl.col("MaxCharges")),
        ("MinIsotopes", pl.col("MinIsotopes")),
        ("MaxIsotopes", pl.col("MaxIsotopes")),
        ("CosineScore", pl.col("cos")),
        ("SNR", pl.col("snr")),
        ("QScore", pl.col("qscore")),
    ]
    return _explode_long_by_position(pl_deconv_indexed, "mass_id", value_exprs)

def threedim_SN_plot(pl_deconv_indexed):
    """3D signal/noise scatter frame with precursor-lookup keys.

    Columns: index, Scan (i64), PrecursorScan (f64), PrecursorMass (f64),
    MonoMass (list[f64], == ``mz_array``), SignalPeaks, NoisyPeaks.

    Carries the precursor-lookup keys the OpenMS-Insight Scatter3D needs to match
    a fragment scan's precursor back to the mass that generated it: each (MS2) scan
    row has its own ``Scan`` id, the ``PrecursorScan`` it was isolated from, the
    scalar ``PrecursorMass``, and the per-mass ``MonoMass`` array of that scan. The
    viewer locates the precursor scan's row (``Scan == PrecursorScan``) and finds
    the index into that scan's ``MonoMass`` array whose value matches
    ``PrecursorMass``, using the same per-mass position (``SignalPeaks[massIndex]``)
    as the rest of the 3D plot. ``Scan`` / ``PrecursorScan`` are the join keys,
    ``PrecursorMass`` is the scalar to match, ``MonoMass`` is the per-mass axis to
    search.
    """
    return (
        pl_deconv_indexed
        .select([
            pl.col('index'),
            pl.col('Scan'),
            pl.col('PrecursorScan'),
            pl.col('PrecursorMass'),
            pl.col('mz_array').alias('MonoMass'),
            pl.col('SignalPeaks'),
            pl.col('NoisyPeaks')
        ])
    )


def parseDeconv(
        file_manager, dataset_id, out_deconv_mzML, anno_annotated_mzML, 
        spec1_tsv=None, spec2_tsv=None, logger=None
):
    logger.log("Progress of 'processing FLASHDeconv results':", level=2)
    logger.log("0.0 %", level=2)

    # Parse input files
    tolerance = parseFLASHDeconvOutput(
        anno_annotated_mzML, out_deconv_mzML,
        file_manager, dataset_id, logger=logger,
    )
    file_manager.store_data(dataset_id, 'deconv_tolerance', float(tolerance))
    
    # Immediately reload as polars LazyFrames for efficient processing
    results = file_manager.get_results(dataset_id, ['anno_dfs', 'deconv_dfs'], use_polars=True)
    pl_anno = results['anno_dfs']
    pl_deconv = results['deconv_dfs']
    
    logger.log("10.0 %", level=2)

    # Preprocess data for the heatmaps
    for df, descriptor in zip([pl_deconv, pl_anno], ['deconv', 'raw']):

        # Create full sized version - returns polars LazyFrame
        heatmap_lazy = getMSSignalDF(df)

        for ms_level in [1, 2]:
            
            # Filter for specific MS level using polars operations
            relevant_heatmap_lazy = (
                heatmap_lazy
                .filter(pl.col('MSLevel') == ms_level)
                .drop('MSLevel')
            )

            # Collect here as this is the data we are operating on
            relevant_heatmap_lazy = relevant_heatmap_lazy.collect().lazy()

            # Get count for compression level calculation
            heatmap_count = relevant_heatmap_lazy.select(pl.len()).collect().item()

            # Store full sized version
            file_manager.store_data(
                dataset_id, f'ms{ms_level}_{descriptor}_heatmap',
                relevant_heatmap_lazy
            )

            # Store compressed versions
            compression_levels = compute_compression_levels(20000, heatmap_count, logger=logger)
            current_heatmap_lazy = relevant_heatmap_lazy
            
            for size in reversed(compression_levels):
                # Downsample iteratively using polars-optimized function
                current_heatmap_lazy = downsample_heatmap(current_heatmap_lazy, max_datapoints=size)
                
                # Store compressed version - convert to pandas only at storage
                file_manager.store_data(
                    dataset_id, f'ms{ms_level}_{descriptor}_heatmap_{size}',
                    current_heatmap_lazy
                )
    
    logger.log("20.0 %", level=2)
        
    # scan_table - using native polars operations
    spectra_lazy = (
        pl_deconv
        .with_row_index("index")
        .with_columns([
            pl.col('MinCharges').list.len().alias('#Masses')
        ])
        .select([
            pl.col('index'),
            pl.col('Scan'),
            pl.col('MSLevel'),
            pl.col('rt').alias('RT'),
            pl.col('PrecursorMass'),
            pl.col('#Masses')
        ])
        .sort("index")
    )
    file_manager.store_data(dataset_id, 'scan_table', spectra_lazy)

    logger.log("30.0 %", level=2)

    # Add row indices for joining operations
    pl_deconv_indexed = pl_deconv.with_row_index("index")
    pl_anno_indexed = pl_anno.with_row_index("index")

    # anno_spectrum - using native polars LazyFrame operations
    anno_spectrum_lazy = (
        pl_anno_indexed
        .select([
            pl.col('index'),
            pl.col('mz_array').alias('MonoMass_Anno'),
            pl.col('intensity_array').alias('SumIntensity_Anno')
        ])
        .sort("index")
    )
    file_manager.store_data(dataset_id, 'anno_spectrum', anno_spectrum_lazy, row_group_size=SPECTRA_ROW_GROUP_SIZE)

    # anno_spectrum_long - long-format (one row per peak) for OpenMS-Insight LinePlot
    file_manager.store_data(
        dataset_id, 'anno_spectrum_long',
        anno_spectrum_long(pl_anno_indexed),
        row_group_size=LONG_ROW_GROUP_SIZE,
    )

    logger.log("40.0 %", level=2)

    # mass_table - using native polars LazyFrame operations
    mass_table_lazy = (
        pl_deconv_indexed
        .select([
            pl.col('index'),
            pl.col('mz_array').alias('MonoMass'),
            pl.col('intensity_array').alias('SumIntensity'),
            pl.col('MinCharges'),
            pl.col('MaxCharges'),
            pl.col('MinIsotopes'),
            pl.col('MaxIsotopes'),
            pl.col('cos').alias('CosineScore'),
            pl.col('snr').alias('SNR'),
            pl.col('qscore').alias('QScore')
        ])
        .sort("index")
    )
    file_manager.store_data(dataset_id, 'mass_table', mass_table_lazy, row_group_size=SPECTRA_ROW_GROUP_SIZE)

    # mass_table_long - long-format (one row per mass) for OpenMS-Insight Table
    file_manager.store_data(
        dataset_id, 'mass_table_long',
        mass_table_long(pl_deconv_indexed),
        row_group_size=LONG_ROW_GROUP_SIZE,
    )

    logger.log("50.0 %", level=2)

    # sequence_view - using native polars LazyFrame operations
    sequence_view_lazy = (
        pl_deconv_indexed
        .select([
            pl.col('index'),
            pl.col('mz_array').alias('MonoMass'),
            pl.col('PrecursorMass')
        ])
        .sort("index")
    )
    file_manager.store_data(dataset_id, 'sequence_view', sequence_view_lazy, row_group_size=SPECTRA_ROW_GROUP_SIZE)

    logger.log("60.0 %", level=2)

    # deconv_spectrum - using native polars LazyFrame operations
    deconv_spectrum_lazy = (
        pl_deconv_indexed
        .select([
            pl.col('index'),
            pl.col('mz_array').alias('MonoMass'),
            pl.col('intensity_array').alias('SumIntensity')
        ])
        .sort("index")
    )
    file_manager.store_data(dataset_id, 'deconv_spectrum', deconv_spectrum_lazy, row_group_size=SPECTRA_ROW_GROUP_SIZE)

    # deconv_spectrum_long - long-format (one row per peak) for OpenMS-Insight LinePlot
    file_manager.store_data(
        dataset_id, 'deconv_spectrum_long',
        deconv_spectrum_long(pl_deconv_indexed),
        row_group_size=LONG_ROW_GROUP_SIZE,
    )

    logger.log("70.0 %", level=2)

    # anno & deconv spectrum (combined_spectrum) - using native polars LazyFrame join
    combined_spectrum_lazy = (
        pl_deconv_indexed
        .select([
            pl.col('index'),
            pl.col('mz_array').alias('MonoMass'),
            pl.col('intensity_array').alias('SumIntensity'),
            pl.col('SignalPeaks')
        ])
        .join(
            pl_anno_indexed.select([
                pl.col('index'),
                pl.col('mz_array').alias('MonoMass_Anno'),
                pl.col('intensity_array').alias('SumIntensity_Anno')
            ]),
            on='index',
            how='left'
        )
        .sort("index")
    )
    file_manager.store_data(dataset_id, 'combined_spectrum', combined_spectrum_lazy, row_group_size=SPECTRA_ROW_GROUP_SIZE)

    # combined_spectrum_long - long-format deconv peaks + is_signal flag for
    # OpenMS-Insight LinePlot (primary series). The annotated overlay (2nd series)
    # is the separate anno_spectrum_long frame wired via x2_column/y2_column.
    file_manager.store_data(
        dataset_id, 'combined_spectrum_long',
        combined_spectrum_long(pl_deconv_indexed),
        row_group_size=LONG_ROW_GROUP_SIZE,
    )

    logger.log("80.0 %", level=2)

    # 3D_SN_plot - precursor-lookup keys + signal/noise peaks (see threedim_SN_plot).
    threedim_SN_plot_lazy = threedim_SN_plot(pl_deconv_indexed)
    file_manager.store_data(dataset_id, 'threedim_SN_plot', threedim_SN_plot_lazy)

    logger.log("90.0 %", level=2)

    # fdr_plot
    fdr_dfs = []
    if spec1_tsv is not None:
        fdr_dfs.append(pd.read_csv(spec1_tsv, sep='\t'))
    if spec2_tsv is not None:
        fdr_dfs.append(pd.read_csv(spec2_tsv, sep='\t'))
    if len(fdr_dfs) > 0:
        fdr_dfs = pd.concat(fdr_dfs, axis=0, ignore_index=True)
        if 'TargetDecoyType' not in fdr_dfs.columns:
            fdr_dfs['TargetDecoyType'] = 0
        density_target, density_decoy = fdr_density_distribution(fdr_dfs)
        file_manager.store_data(dataset_id, 'density_target', density_target)
        file_manager.store_data(dataset_id, 'density_decoy', density_decoy)
    
    logger.log("100.0 %", level=2)


def fdr_density_distribution(df):

    # Find density targets
    target_qscores = df[df['TargetDecoyType'] == 0]['Qscore'].dropna()
    x_target = np.linspace(target_qscores.min(), target_qscores.max(), 200)
    kde_target = gaussian_kde(target_qscores)
    density_target = pd.DataFrame({'x': x_target, 'y': kde_target(x_target)})

    # Find density decoys (if present)
    decoy_qscores = df[df['TargetDecoyType'] > 0]['Qscore'].dropna()
    if len(decoy_qscores) > 0:
        x_decoy = np.linspace(decoy_qscores.min(), decoy_qscores.max(), 200)
        kde_decoy = gaussian_kde(decoy_qscores)
        density_decoy = pd.DataFrame({'x': x_decoy, 'y': kde_decoy(x_decoy)})
    else:
        density_decoy = pd.DataFrame(columns=['x', 'y'])

    return density_target, density_decoy