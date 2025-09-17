import numpy as np
import pandas as pd
import polars as pl

from scipy.stats import binned_statistic_2d


def compute_compression_levels(minimal_size, total_amount, logger=None):
    if total_amount <= minimal_size:
        return []
    levels = np.logspace(
        int(np.log10(minimal_size)), int(np.log10(total_amount)),
        int(np.log10(total_amount)) - int(np.log10(minimal_size)) + 1,
        dtype='int'
    )*int(10**(np.log10(minimal_size)%1))
    return levels[levels<total_amount]


def downsample_heatmap(data, max_datapoints=20000, rt_bins=400, mz_bins=50, logger=None):
    """
    Downsample heatmap data using polars for efficient processing.
    
    Args:
        data: polars LazyFrame or DataFrame with columns ['mass', 'rt', 'intensity', ...]
        max_datapoints: Maximum number of points to keep
        rt_bins: Number of retention time bins
        mz_bins: Number of mass bins
        logger: Optional logger
    
    Returns:
        polars LazyFrame with downsampled data
    """
    if (rt_bins * mz_bins) > max_datapoints:
        raise ValueError("Number of bins more than maximum datapoints.")
    
    # Ensure we're working with a LazyFrame
    if isinstance(data, pl.DataFrame):
        data = data.lazy()
    elif isinstance(data, pd.DataFrame):
        data = pl.from_pandas(data).lazy()
    
    # Sort by rt and intensity for ranking
    sorted_data = (
        data
        .sort(['rt', 'intensity'], descending=[False, True])
        .with_columns([
            pl.int_range(pl.len()).over('rt').alias('rank')
        ])
        .sort(['rank', 'intensity'], descending=[False, True])
    )

    # We need to collect here because scipy requires numpy arrays
    sorted_data = sorted_data.collect()
    
    # Count peaks
    total_count = sorted_data.select(pl.count()).item()
    
    # Extract arrays for scipy
    mass_array = sorted_data['mass'].to_numpy()
    rt_array = sorted_data['rt'].to_numpy()
    intensity_array = sorted_data['intensity'].to_numpy()
    
    # Use scipy for binning (still needed for the specific binning logic)
    count, _, __, mapping = binned_statistic_2d(
        mass_array, rt_array, intensity_array, 'count',
        [mz_bins, rt_bins], expand_binnumbers=True
    )

    # Back to lazy evaluation
    sorted_data = sorted_data.lazy()
    
    # Add bin information back to polars DataFrame
    binned_data = (
        sorted_data
        .with_columns([
            pl.Series('mass_bin', mapping[0] - 1),  # scipy uses 1-based indexing
            pl.Series('rt_bin', mapping[1] - 1)
        ])
    )
    
    # Compute maximum amount of peaks per bin that does not exceed limit
    counted_peaks = 0
    max_peaks_per_bin = -1
    new_count = 0
    while ((counted_peaks + new_count) < max_datapoints):
        # commit prev result
        max_peaks_per_bin += 1
        counted_peaks += new_count
        # compute count for next value
        new_count = np.sum(count.flatten() >= (max_peaks_per_bin + 1))

        if counted_peaks >= total_count:
            break
    
    # Use polars for efficient groupby and head operations
    result = (
        binned_data
        .group_by(['mass_bin', 'rt_bin'])
        .head(max_peaks_per_bin)
        .sort('intensity')
        .drop(['rank', 'mass_bin', 'rt_bin'])
    )
    
    return result
