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
    
    logger.log("10.0 %", level=2)

    # Preprocess data for the heatmaps
    for df, descriptor in zip([deconv_df, anno_df], ['deconv', 'raw']):

        # Create full sized version - returns polars LazyFrame
        heatmap_lazy = getMSSignalDF(df)

        for ms_level in [1, 2]:
            
            # Filter for specific MS level using polars operations
            relevant_heatmap_lazy = (
                heatmap_lazy
                .filter(pl.col('MSLevel') == ms_level)
                .drop('MSLevel')
            )

            # Get count for compression level calculation
            heatmap_count = relevant_heatmap_lazy.select(pl.len()).collect().item()

            # Store full sized version - convert to pandas only at storage
            file_manager.store_data(
                dataset_id, f'ms{ms_level}_{descriptor}_heatmap',
                relevant_heatmap_lazy.collect().to_pandas()
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
                    current_heatmap_lazy.collect().to_pandas()
                )
    
    logger.log("20.0 %", level=2)
        
    spectra_df = getSpectraTableDF(deconv_df)

    # scan_table
    scan_table = spectra_df.loc[
        :,['index', 'Scan', 'MSLevel', 'RT', 'PrecursorMass', '#Masses']
    ]
    file_manager.store_data(dataset_id, 'scan_table', scan_table)

    logger.log("30.0 %", level=2)

    # Convert to polars for efficient operations
    pl_deconv = pl.from_pandas(deconv_df.reset_index(names=['index']))
    pl_anno = pl.from_pandas(anno_df.reset_index(names=['index']))
    
    # Create index-only table for joining
    index_df = pl_deconv.select(pl.col('index'))

    # anno_spectrum
    anno_spectrum = (
        pl_anno
        .select([
            pl.col('index'),
            pl.col('mzarray').alias('MonoMass_Anno'),
            pl.col('intarray').alias('SumIntensity_Anno')
        ])
        .to_pandas()
    )
    file_manager.store_data(dataset_id, 'anno_spectrum', anno_spectrum)

    logger.log("40.0 %", level=2)

    # mass_table
    mass_table = (
        pl_deconv
        .select([
            pl.col('index'),
            pl.col('mzarray').alias('MonoMass'),
            pl.col('intarray').alias('SumIntensity'),
            pl.col('MinCharges'),
            pl.col('MaxCharges'),
            pl.col('MinIsotopes'),
            pl.col('MaxIsotopes'),
            pl.col('cos').alias('CosineScore'),
            pl.col('snr').alias('SNR'),
            pl.col('qscore').alias('QScore')
        ])
        .to_pandas()
    )
    file_manager.store_data(dataset_id, 'mass_table', mass_table)

    logger.log("50.0 %", level=2)

    # sequence_view
    sequence_view = (
        pl_deconv
        .select([
            pl.col('index'),
            pl.col('mzarray').alias('MonoMass'),
            pl.col('PrecursorMass')
        ])
        .to_pandas()
    )
    file_manager.store_data(dataset_id, 'sequence_view', sequence_view)

    logger.log("60.0 %", level=2)

    # deconv_spectrum
    deconv_spectrum = (
        pl_deconv
        .select([
            pl.col('index'),
            pl.col('mzarray').alias('MonoMass'),
            pl.col('intarray').alias('SumIntensity')
        ])
        .to_pandas()
    )
    file_manager.store_data(dataset_id, 'deconv_spectrum', deconv_spectrum)

    logger.log("70.0 %", level=2)

    # anno & deconv spectrum (combined_spectrum)
    combined_spectrum = (
        pl_deconv
        .select([
            pl.col('index'),
            pl.col('mzarray').alias('MonoMass'),
            pl.col('intarray').alias('SumIntensity'),
            pl.col('SignalPeaks')
        ])
        .join(
            pl_anno.select([
                pl.col('index'),
                pl.col('mzarray').alias('MonoMass_Anno'),
                pl.col('intarray').alias('SumIntensity_Anno')
            ]),
            on='index',
            how='left'
        )
        .to_pandas()
    )
    file_manager.store_data(dataset_id, 'combined_spectrum', combined_spectrum)

    logger.log("80.0 %", level=2)

    # 3D_SN_plot
    threedim_SN_plot = (
        pl_deconv
        .select([
            pl.col('index'),
            pl.col('PrecursorScan'),
            pl.col('SignalPeaks'),
            pl.col('NoisyPeaks')
        ])
        .to_pandas()
    )
    file_manager.store_data(dataset_id, 'threedim_SN_plot', threedim_SN_plot)

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