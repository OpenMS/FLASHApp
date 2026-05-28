import json

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from src.render.sequence_data_store import build_table, ROW_GROUP_SIZE

from io import StringIO
from pyopenms import AASequence
from scipy.stats import gaussian_kde

from src.parse.masstable import parseFLASHTaggerOutput
from src.parse.tag_resolution import build_tagspace_to_proteoform_map
from src.render.sequence import (
    remove_ambigious, getFragmentDataFromSeq, getInternalFragmentDataFromSeq
)


def _coverage_from_ranges(starts, ends, seq_len):
    """Per-residue coverage: number of tags whose inclusive [StartPos, EndPos]
    range contains the position. Equivalent to the old per-position sum, but
    walks each tag once. StartPos/EndPos arrive as float64 (the column carries
    NaN for unplaced tags); NaN tags are skipped, matching the old comparison
    where NaN<=i and NaN>=i are both False. Integer-valued positions cast
    exactly; StartPos=-1 clamps to 0; numpy caps the slice stop."""
    coverage = np.zeros(seq_len, dtype='float')
    for s, e in zip(starts, ends):
        if s != s or e != e:          # NaN positions contribute nothing
            continue
        coverage[(int(s) if s > 0 else 0):int(e) + 1] += 1
    return coverage


def _linearize_tag_df(tag_df):
    """Expand ';'-packed multi-proteoform tag rows to one row per (tag, proteoform).

    A tag can match N proteoforms; proteoform-specific fields are ';'-packed.
    Only rows whose ProteoformIndex contains ';' expand; the rest pass through
    untouched. Both groups are written to one TSV buffer and re-read, so the
    roundtrip alone determines dtypes — matching the original code exactly
    (the regression comparator checks dtypes strictly; a pd.concat here would
    upcast int->float when the split frame is empty)."""
    is_split = tag_df['ProteoformIndex'].astype(str).str.contains(';', regex=False)

    non_split = tag_df[~is_split].copy()
    non_split['ProteoformIndex'] = non_split['ProteoformIndex'].fillna(-1)

    expanded = {c: [] for c in tag_df.columns}
    for _, row in tag_df[is_split].iterrows():
        no_items = row['ProteoformIndex'].count(';') + 1
        for c in tag_df.columns:
            if isinstance(row[c], str) and (';' in row[c]):
                expanded[c] += row[c].split(';')
            else:
                expanded[c] += [row[c]] * no_items
    split_expanded = pd.DataFrame(expanded, columns=tag_df.columns)

    buf = StringIO()
    non_split.to_csv(buf, sep='\t', index=False)
    if len(split_expanded):
        split_expanded.to_csv(buf, sep='\t', index=False, header=False)
    buf.seek(0)
    return pd.read_csv(buf, sep='\t')


def parseTnT(file_manager, dataset_id, deconv_mzML, anno_mzML, tag_tsv, protein_tsv, logger=None):
    logger.log("Progress of 'processing FLASHTnT results':", level=2)
    logger.log("0.0 %", level=2)

    tolerance = file_manager.get_results(dataset_id, ['deconv_tolerance'])['deconv_tolerance']
    tag_df, protein_df = parseFLASHTaggerOutput(tag_tsv, protein_tsv)
    # Map FLASHTagger tag-space ProteoformIndex -> protein-space index from the
    # raw frames (before protein_df is renamed and tag_df is linearized).
    tagspace_to_proteoform = build_tagspace_to_proteoform_map(tag_df, protein_df)
    logger.log("10.0 %", level=2)
    
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
    logger.log("30.0 %", level=2)

    # tag_table

    # Process tag df into a linear data format
    tag_df = _linearize_tag_df(tag_df)

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
    logger.log("50.0 %", level=2)

    # sequence_view & internal_fragment_map
    sequence_data = {}
    # internal_fragment_data = {}  # Disabled
    # Compute coverage
    # tag_df['ProteinIndex'] is tag-space; map to protein-space so coverage uses
    # each proteoform's own tags (the two enumerations diverge on large runs).
    proteoform_of_tag = tag_df['ProteinIndex'].map(
        lambda q: tagspace_to_proteoform.get(int(q), -1) if pd.notna(q) else -1
    )
    tag_groups = {
        pid: (g['StartPos'].to_numpy(), g['EndPos'].to_numpy())
        for pid, g in tag_df.groupby(proteoform_of_tag)[['StartPos', 'EndPos']]
    }
    for i, row in protein_df.iterrows():
        pid = row['index']
        sequence = row['sequence']
        L = len(sequence)
        if pid in tag_groups:
            starts, ends = tag_groups[pid]
            coverage = _coverage_from_ranges(starts, ends, L)
        else:
            coverage = np.zeros(L, dtype='float')
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

        # internal_fragment_data[pid] = getInternalFragmentDataFromSeq(
        #     str(sequence)[start_index:end_index+1], modifications
        # )  # Disabled

    sequence_data_table = build_table(sequence_data)
    with file_manager.parquet_sink(dataset_id, 'sequence_data') as sequence_data_path:
        pq.write_table(sequence_data_table, sequence_data_path, row_group_size=ROW_GROUP_SIZE)
    # file_manager.store_data(
    #     dataset_id, 'internal_fragment_data', internal_fragment_data
    # )  # Disabled
    logger.log("70.0 %", level=2)

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
    logger.log("90.0 %", level=2)

    density_target, density_decoy = fdr_density_distribution(protein_df, logger=logger)
    file_manager.store_data(dataset_id, 'density_id_target', density_target)
    file_manager.store_data(dataset_id, 'density_id_decoy', density_decoy)
    logger.log("100.0 %", level=2)


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