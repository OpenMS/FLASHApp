import pandas as pd
import numpy as np

from src.parse.masstable import parseFLASHDeconvOutput, getMSSignalDF, getSpectraTableDF
from src.render.compression import downsample_heatmap, compute_compression_levels
from scipy.stats import gaussian_kde

def parseIda(
        file_manager, dataset_id, out_simulation, logger=None
):
    
    simulation_data = pd.read_csv(out_simulation, sep='\t')
    file_manager.store_data(
        dataset_id, 'simulation_dfs', simulation_data
    )

    heatmap = simulation_data.loc[:,['monoMasses', 'rt', 'precursorIntensity']]
    heatmap = heatmap.reset_index()
    heatmap = simulation_data.rename(columns={
        'monoMasses' : 'mass', 'precursorIntensity' : 'intensity', 
        'index' : 'scan_idx'
    })

    # Store full sized version
    file_manager.store_data(
        dataset_id, f'ms1_deconv_heatmap', heatmap
    )

    # Store compressed versions
    for size in reversed(compute_compression_levels(20000, len(heatmap), logger=logger)):
        
        
        # Downsample iteratively
        heatmap = downsample_heatmap(heatmap, max_datapoints=size)
        # Store compressed version
        file_manager.store_data(
            dataset_id, f'ms1_deconv_heatmap_{size}', heatmap
        )