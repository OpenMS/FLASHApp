import numpy as np
import pandas as pd

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

    if (rt_bins * mz_bins) > max_datapoints:
        raise ValueError("Number of bins more than maximum datapoints.")
    
    data = data.sort_values(['rt', 'intensity'], ascending=[True, False])
    data['rank'] = data.groupby('rt').cumcount()
    data = data.sort_values(['rank', 'intensity'], ascending=[True, False])

    count, _, __, mapping = binned_statistic_2d(
        data['mass'], data['rt'], data['intensity'], 'count', 
        [mz_bins, rt_bins], expand_binnumbers=True
    )
    data['mass_bin'] = mapping[0]
    data['rt_bin'] = mapping[1]

    
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

        if counted_peaks >= len(data):
            break


    data = data.groupby(
        ['mass_bin', 'rt_bin'], group_keys=False, sort=False
    ).head(max_peaks_per_bin).reset_index(drop=True)

    return data.sort_values(by='intensity', ascending=True).drop(
        columns=['rank', 'mass_bin', 'rt_bin']
    )
