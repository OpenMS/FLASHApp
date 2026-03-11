import pandas as pd
import polars as pl
import numpy as np

from src.parse.masstable import parseFLASHDeconvOutput, getMSSignalDF, getSpectraTableDF
from src.render.compression import downsample_heatmap, compute_compression_levels
from scipy.stats import gaussian_kde

def parseDeconv(
        file_manager, dataset_id, out_deconv_mzML, anno_annotated_mzML, 
        spec1_tsv=None, spec2_tsv=None, logger=None
):
    logger.log("Progress of 'processing FLASHDeconv results':", level=2)
    logger.log("0.0 %", level=2)

    # Parse input files
    deconv_df, anno_df, _, _, _ = parseFLASHDeconvOutput(
        anno_annotated_mzML, out_deconv_mzML, logger=logger
    )
    file_manager.store_data(dataset_id, 'anno_dfs', anno_df)
    file_manager.store_data(dataset_id, 'deconv_dfs', deconv_df)
    del deconv_df
    del anno_df
    
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
    file_manager.store_data(dataset_id, 'anno_spectrum', anno_spectrum_lazy)

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
    file_manager.store_data(dataset_id, 'mass_table', mass_table_lazy)

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
    file_manager.store_data(dataset_id, 'sequence_view', sequence_view_lazy)

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
    file_manager.store_data(dataset_id, 'deconv_spectrum', deconv_spectrum_lazy)

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
    file_manager.store_data(dataset_id, 'combined_spectrum', combined_spectrum_lazy)

    logger.log("80.0 %", level=2)

    # 3D_SN_plot - using native polars LazyFrame operations
    threedim_SN_plot_lazy = (
        pl_deconv_indexed
        .select([
            pl.col('index'),
            pl.col('PrecursorScan'),
            pl.col('SignalPeaks'),
            pl.col('NoisyPeaks')
        ])
    )
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