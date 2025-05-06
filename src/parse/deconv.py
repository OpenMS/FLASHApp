import pandas as pd
import numpy as np

from src.parse.masstable import parseFLASHDeconvOutput, getMSSignalDF, getSpectraTableDF
from src.render.compression import downsample_heatmap, compute_compression_levels


def parseDeconv(
        file_manager, dataset_id, out_deconv_mzML, anno_annotated_mzML, 
        spec1_tsv=None, spec2_tsv=None, logger=None
):
    # Parse input files
    deconv_df, anno_df, _, _, _ = parseFLASHDeconvOutput(
        anno_annotated_mzML, out_deconv_mzML, logger=logger
    )

    file_manager.store_data(dataset_id, 'anno_dfs', anno_df)
    file_manager.store_data(dataset_id, 'deconv_dfs', deconv_df)
    # Preprocess data for the heatmaps
    for df, descriptor in zip([deconv_df, anno_df], ['deconv', 'raw']):

        # Create full sized version
        heatmap = getMSSignalDF(df)

        # Store full sized version
        file_manager.store_data(
            dataset_id, f'ms1_{descriptor}_heatmap', heatmap
        )
        # Store compressed versions
        for size in reversed(compute_compression_levels(20000, len(heatmap), logger=logger)):
            
            
            # Downsample iteratively
            heatmap = downsample_heatmap(heatmap, max_datapoints=size)
            # Store compressed version
            file_manager.store_data(
                dataset_id, f'ms1_{descriptor}_heatmap_{size}', heatmap
            )
    
    spectra_df = getSpectraTableDF(deconv_df)

    # scan_table
    scan_table = spectra_df.loc[
        :,['index', 'Scan', 'MSLevel', 'RT', 'PrecursorMass', '#Masses']
    ]
    file_manager.store_data(dataset_id, 'scan_table', scan_table)

    # Subsequent tables only share index
    scan_table = scan_table.loc[:, ['index']]

    # anno_spectrum
    anno_spectrum = anno_df.loc[:,['mzarray', 'intarray']]
    anno_spectrum.rename(columns={'mzarray': 'MonoMass_Anno', 'intarray': 'SumIntensity_Anno'},
                            inplace=True)
    anno_spectrum = pd.concat([scan_table, anno_spectrum], axis=1)
    file_manager.store_data(dataset_id, 'anno_spectrum', anno_spectrum)

    # mass_table
    mass_table = deconv_df.loc[
        :,['mzarray', 'intarray', 'MinCharges', 'MaxCharges', 'MinIsotopes', 'MaxIsotopes', 'cos', 'snr', 'qscore']
    ]
    mass_table.rename(columns={'mzarray': 'MonoMass', 'intarray': 'SumIntensity', 'cos': 'CosineScore',
                                    'snr': 'SNR', 'qscore': 'QScore'},
                            inplace=True)
    mass_table = pd.concat([scan_table, mass_table], axis=1)
    file_manager.store_data(dataset_id, 'mass_table', mass_table)

    # sequence_view
    sequence_view = deconv_df.loc[:, ['mzarray', 'PrecursorMass']]
    sequence_view.rename(columns={'mzarray': 'MonoMass'}, inplace=True)
    sequence_view = pd.concat([scan_table, sequence_view], axis=1)
    file_manager.store_data(dataset_id, 'sequence_view', sequence_view)

    # deconv_spectrum
    deconv_spectrum = deconv_df.loc[
        :,['mzarray', 'intarray']
    ]
    deconv_spectrum.rename(columns={'mzarray': 'MonoMass', 'intarray': 'SumIntensity'},
                            inplace=True)
    deconv_spectrum = pd.concat([scan_table, deconv_spectrum], axis=1)
    file_manager.store_data(dataset_id, 'deconv_spectrum', deconv_spectrum)

    # anno & deconv spectrum
    combined_spectrum = pd.concat(
        [deconv_spectrum, anno_spectrum.drop(columns=['index']), 
         deconv_df.loc[:, ['SignalPeaks']]],
        axis=1
    )
    file_manager.store_data(dataset_id, 'combined_spectrum', combined_spectrum)

    # 3D_SN_plot
    threedim_SN_plot = deconv_df.loc[
        :, ['PrecursorScan', 'SignalPeaks', 'NoisyPeaks']
    ]
    threedim_SN_plot = pd.concat([scan_table, threedim_SN_plot], axis=1)
    file_manager.store_data(dataset_id, 'threedim_SN_plot', threedim_SN_plot)

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
        ecdf_target, ecdf_decoy = ecdf(fdr_dfs)
        file_manager.store_data(dataset_id, 'ecdf_target', ecdf_target)
        file_manager.store_data(dataset_id, 'ecdf_decoy', ecdf_decoy)

    
def ecdf(df):
    target_qscores = df[df['TargetDecoyType'] == 0]['Qscore']
    decoy_qscores = df[df['TargetDecoyType'] > 0]['Qscore']

    ecdf_target = pd.DataFrame({
        'x' : np.sort(target_qscores),
        'y' : np.arange(1, len(target_qscores) + 1) / len(target_qscores)
    })
    ecdf_decoy = pd.DataFrame({
        'x' : np.sort(decoy_qscores),
        'y' : np.arange(1, len(decoy_qscores) + 1) / len(decoy_qscores)
    })
    return ecdf_target, ecdf_decoy