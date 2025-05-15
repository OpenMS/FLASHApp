import json

import numpy as np
import pandas as pd

from io import StringIO
from pyopenms import AASequence
from scipy.stats import gaussian_kde

from src.parse.masstable import parseFLASHDeconvOutput, parseFLASHTaggerOutput
from src.render.sequence import (
    remove_ambigious, getFragmentDataFromSeq, getInternalFragmentDataFromSeq
)


def parseTnT(file_manager, dataset_id, deconv_mzML, anno_mzML, tag_tsv, protein_tsv, logger=None):

    deconv_df, _, tolerance, _, _,  = parseFLASHDeconvOutput(
        anno_mzML, deconv_mzML
    )
    tag_df, protein_df = parseFLASHTaggerOutput(tag_tsv, protein_tsv)
    
    # protein_table
    protein_df['length'] = protein_df['DatabaseSequence'].apply(lambda x : len(x))
    protein_df = protein_df.rename(
        columns={
            'ProteoformIndex' : 'index',
            'ProteinAccession' : 'accession',
            'ProteinDescription' : 'description',
            'DatabaseSequence' : 'sequence'
        }
    )
    protein_df['description'] = protein_df['description'].apply(
        lambda x: x[:50] + '...' if len(x) > 50 else x
    )
    file_manager.store_data(dataset_id, 'protein_dfs', protein_df)

    # tag_table

    # Process tag df into a linear data format
    new_tag_df = {c : [] for c in tag_df.columns}
    for i, row in tag_df.iterrows():
        # No splitting if it is not recognized as string
        if pd.isna(row['ProteoformIndex']):
            row['ProteoformIndex'] = -1
        if isinstance(row['ProteoformIndex'], str) and (';' in row['ProteoformIndex']):
            no_items = row['ProteoformIndex'].count(';') + 1
            for c in new_tag_df.keys():
                if (isinstance(row[c], str)) and (';' in row[c]):
                    new_tag_df[c] += row[c].split(';')
                else:
                    new_tag_df[c] += [row[c]]*no_items
        else:
            for c in new_tag_df.keys():
                new_tag_df[c].append(row[c])
    tag_df = pd.DataFrame(new_tag_df)

    tsv_buffer = StringIO()
    tag_df.to_csv(tsv_buffer, sep='\t', index=False)
    tsv_buffer.seek(0)
    tag_df = pd.read_csv(tsv_buffer, sep='\t')

    # Complete df
    tag_df['StartPosition'] = tag_df['StartPosition'] - 1
    tag_df['EndPos'] = tag_df['StartPosition'] + tag_df['Length'] - 1
    tag_df = tag_df.rename(
        columns={
            'ProteoformIndex' : 'ProteinIndex',
            'DeNovoScore' : 'Score',
            'Masses' : 'mzs',
            'StartPosition' : 'StartPos' 
        }
    )
    file_manager.store_data(dataset_id, 'tag_dfs', tag_df)
    # sequence_view & internal_fragment_map
    sequence_data = {}
    internal_fragment_data = {}
    # Compute coverage
    for i, row in protein_df.iterrows():
        pid = row['index']
        sequence = row['sequence']
        coverage = np.zeros(len(sequence), dtype='float')
        for i in range(len(sequence)):
            coverage[i] = np.sum(
                (tag_df['ProteinIndex'] == pid) &
                (tag_df['StartPos'] <= i) &
                (tag_df['EndPos'] >= i)
            )
        p_cov = np.zeros(len(coverage))
        if np.max(coverage) > 0:
            p_cov = coverage/np.max(coverage)

        proteoform_start = row['StartPosition']
        proteoform_end = row['EndPosition']
        start_index = 0 if proteoform_start <= 0 else proteoform_start - 1
        end_index = len(sequence) - 1 if proteoform_end <= 0 else proteoform_end - 1


        if row['ModCount'] > 0:
            mod_masses = [float(m) for m in str(row['ModMass']).split(';')]
            mod_starts = [int(float(s)) for s in str(row['ModStart']).split(';')]
            mod_ends = [int(float(s)) for s in str(row['ModEnd']).split(';')]
            if pd.isna(row['ModID']):
                mod_labels = [''] * row['ModCount']
            else:
                mod_labels = [s[:-1].replace(',', '; ') for s in str(row['ModID']).split(';')]
        else:
            mod_masses = []
            mod_starts = []
            mod_ends = []
            mod_labels = []
        modifications = []
        for s, e, m in zip(mod_starts, mod_ends, mod_masses):
            modifications.append((s-start_index, e-start_index, m))
        
        sequence = str(sequence)
        sequence_data[pid] = getFragmentDataFromSeq(
            str(sequence)[start_index:end_index+1], p_cov, np.max(coverage), 
            modifications
        )

        sequence_data[pid]['sequence'] = list(sequence)
        sequence_data[pid]['proteoform_start'] = proteoform_start - 1
        sequence_data[pid]['proteoform_end'] = proteoform_end - 1
        sequence_data[pid]['computed_mass'] = row['ProteoformMass']
        sequence_data[pid]['theoretical_mass'] = remove_ambigious(AASequence.fromString(sequence)).getMonoWeight()
        sequence_data[pid]['modifications'] = [
            {
                # Modfications are zero based
                'start' : s - 1,
                'end' : e - 1,
                'mass_diff' : m,
                'labels' : l
            } for s, e, m, l in zip(mod_starts, mod_ends, mod_masses, mod_labels)
        ]

        internal_fragment_data[pid] = getInternalFragmentDataFromSeq(
            str(sequence)[start_index:end_index+1], modifications
        )  

    file_manager.store_data(dataset_id, 'sequence_data', sequence_data)
    file_manager.store_data(
        dataset_id, 'internal_fragment_data', internal_fragment_data
    )

    fragments = ['b', 'y']
    if file_manager.result_exists(dataset_id, 'FTnT_parameters_json'):
        tnt_settings_file = file_manager.get_results(
            dataset_id, ['FTnT_parameters_json']
        )['FTnT_parameters_json']
        with open(tnt_settings_file, 'r') as f:
            tnt_settings = json.load(f)
        if 'ion_type' in tnt_settings:
            fragments = tnt_settings['ion_type'].split('\n')
    settings = {
        'tolerance' : tolerance,
        'ion_types' : fragments
    }
    file_manager.store_data(
        dataset_id, 'settings', settings
    )

    density_target, density_decoy = fdr_density_distribution(protein_df, logger=logger)
    file_manager.store_data(dataset_id, 'density_id_target', density_target)
    file_manager.store_data(dataset_id, 'density_id_decoy', density_decoy)


def fdr_density_distribution(df, logger=None):
    df = df[df['ProteoformLevelQvalue'] > 0]
    # Find density targets
    target_qscores = df[~df['accession'].str.startswith('DECOY_')]['ProteoformLevelQvalue'].dropna()
    if len(target_qscores) > 0:
        x_target = np.linspace(target_qscores.min(), target_qscores.max(), 200)
        kde_target = gaussian_kde(target_qscores)
        density_target = pd.DataFrame({'x': x_target, 'y': kde_target(x_target)})
    else:
        density_target = pd.DataFrame(columns=['x', 'y'])

    # Find density decoys (if present)
    decoy_qscores = df[df['accession'].str.startswith('DECOY_')]['ProteoformLevelQvalue'].dropna()
    if len(decoy_qscores) > 0:
        x_decoy = np.linspace(decoy_qscores.min(), decoy_qscores.max(), 200)
        kde_decoy = gaussian_kde(decoy_qscores)
        density_decoy = pd.DataFrame({'x': x_decoy, 'y': kde_decoy(x_decoy)})
    else:
        density_decoy = pd.DataFrame(columns=['x', 'y'])

    return density_target, density_decoy