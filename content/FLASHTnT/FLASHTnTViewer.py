import json

import numpy as np
import pandas as pd
import streamlit as st

from io import StringIO
from pathlib import Path
from pyopenms import AASequence

from src.common.common import page_setup, save_params
from src.masstable import getSpectraTableDF
from src.components import Tabulator, SequenceView, InternalFragmentMap, \
    FlashViewerComponent, flash_viewer_grid_component, PlotlyLineplotTagger
from src.sequence import getFragmentDataFromSeq, getInternalFragmentDataFromSeq
from src.workflow.FileManager import FileManager
from src.sequence import remove_ambigious


DEFAULT_LAYOUT = [
    ['protein_table'], 
    ['sequence_view'], 
    ['tag_table'],
    ['deconv_spectrum']
]

def sendDataToJS(selected_data, layout_info_per_exp, grid_key='flash_viewer_grid'):

    # Get data
    results = file_manager.get_results(
        selected_data, ['deconv_dfs', 'anno_dfs', 'tag_dfs', 'protein_dfs']
    )
    
    # getting data from mzML files
    spec_df = results['deconv_dfs']
    anno_df = results['anno_dfs']
    tag_df = results['tag_dfs']
    protein_df = results['protein_dfs']

    # Get the ion type from configuration if exists
    fragments = ['b', 'y']
    if file_manager.result_exists(selected_data, 'FTnT_parameters_json'):
        tnt_settings_file = file_manager.get_results(
            selected_data, ['FTnT_parameters_json']
        )['FTnT_parameters_json']
        with open(tnt_settings_file, 'r') as f:
            tnt_settings = json.load(f)
        if 'ion_type' in tnt_settings:
            fragments = tnt_settings['ion_type'].split('\n')

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

    protein_df['length'] = protein_df['DatabaseSequence'].apply(lambda x : len(x))
    protein_df = protein_df.rename(
        columns={
            'ProteoformIndex' : 'index',
            'ProteinAccession' : 'accession',
            'ProteinDescription' : 'description',
            'DatabaseSequence' : 'sequence'
        }
    )

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

    components = []
    data_to_send = {}
    per_scan_contents = {'spectrum_view': False, 'sequence_view': False}
    for row in layout_info_per_exp:
        components_of_this_row = []
        for _, comp_name in enumerate(row):
            component_arguments = None

            # prepare component arguments
            if comp_name == 'deconv_spectrum':
                per_scan_contents['spectrum_view'] = True
                component_arguments = PlotlyLineplotTagger(title="Deconvolved Spectrum")
            elif comp_name == 'protein_table':
                data_to_send['protein_table'] = protein_df
                data_to_send['per_scan_data'] = getSpectraTableDF(spec_df)
                component_arguments = Tabulator('ProteinTable')
            elif comp_name == 'tag_table':
                data_to_send['tag_table'] = tag_df
                component_arguments = Tabulator('TagTable')
            elif comp_name == 'sequence_view':
                per_scan_contents['sequence_view'] = True
                component_arguments = SequenceView()
            elif comp_name == 'internal_fragment_map':
                data_to_send['internal_fragment_data'] = internal_fragment_data
                component_arguments = InternalFragmentMap()

            components_of_this_row.append(FlashViewerComponent(component_arguments))
        components.append(components_of_this_row)

    if any(per_scan_contents.values()):
        scan_table = data_to_send['per_scan_data']
        dfs = [scan_table]
        for key, exist in per_scan_contents.items():
            if not exist: continue

            if key == 'spectrum_view':
                tmp_df = spec_df[['mzarray', 'intarray', 'CombinedPeaks', 'SignalPeaks', 'NoisyPeaks']].copy()
                tmp_df.rename(columns={'mzarray': 'MonoMass', 'intarray': 'SumIntensity'}, inplace=True)
                dfs.append(tmp_df)

                tmp_df = anno_df[['mzarray', 'intarray']].copy()
                tmp_df.rename(columns={'mzarray': 'MonoMass_Anno', 'intarray': 'SumIntensity_Anno'}, inplace=True)
                dfs.append(tmp_df)                               
                
            elif key == 'sequence_view':
                # Deconvolved data
                tmp_df = spec_df[['mzarray']].copy()
                tmp_df.rename(columns={'mzarray': 'MonoMass'}, inplace=True)
                dfs.append(tmp_df)

        combined_dfs = pd.concat(dfs, axis=1)
        combined_dfs = combined_dfs[np.isin(combined_dfs['Scan'], protein_df['Scan'])]
        combined_dfs = combined_dfs.loc[:, ~combined_dfs.columns.duplicated()]
        data_to_send['per_scan_data'] = combined_dfs

    # Set sequence data
    data_to_send['sequence_data'] = sequence_data
    data_to_send['settings'] = {
        'tolerance' : spec_df['tol'].to_numpy(dtype='float')[0],
        'ion_types' : fragments
    }
    data_to_send['dataset'] = selected_data
    
    # data_to_send['internal_fragment_data'] = pd.DataFrame(internal_fragment_data)
    flash_viewer_grid_component(components=components, data=data_to_send, component_key=grid_key)


def setSequenceViewInDefaultView():
    if 'input_sequence' in st.session_state and st.session_state.input_sequence:
        global DEFAULT_LAYOUT
        DEFAULT_LAYOUT = DEFAULT_LAYOUT + [['sequence_view']] + [['internal_fragment_map']]

def select_experiment():
    st.session_state.selected_experiment0_tagger = st.session_state.selected_experiment_dropdown_tagger
    if len(layout) > 1:
        for exp_index in range(1, len(layout)):
            if st.session_state[f'selected_experiment_dropdown_{exp_index}_tagger'] is None:
                continue
            st.session_state[f"selected_experiment{exp_index}_tagger"] = st.session_state[f'selected_experiment_dropdown_{exp_index}_tagger']



params = page_setup("TaggerViewer")


# Get available results
file_manager = FileManager(
    st.session_state["workspace"],
    Path(st.session_state['workspace'], 'flashtnt', 'cache')
)
results = file_manager.get_results_list(
    ['deconv_dfs', 'anno_dfs', 'tag_dfs', 'protein_dfs']
)

if file_manager.result_exists('layout', 'layout'):
    layout = file_manager.get_results('layout', 'layout')['layout']
    side_by_side = layout['side_by_side']
    layout = layout['layout']
    
else:
    layout = [DEFAULT_LAYOUT]
    side_by_side = False

### if no input file is given, show blank page
if len(results) == 0:
    st.error('No results to show yet. Please run a workflow first!')
    st.stop()

# Map names to index
name_to_index = {n : i for i, n in enumerate(results)}

if len(layout) == 2 and side_by_side:
    c1, c2 = st.columns(2)
    with c1:
        st.selectbox(
            "choose experiment", results, 
            key="selected_experiment_dropdown_tagger", 
            index=name_to_index[st.session_state.selected_experiment0_tagger] if 'selected_experiment0_tagger' in st.session_state else None,
            on_change=select_experiment
        )
        if 'selected_experiment0_tagger' in st.session_state:
            with st.spinner('Loading component...'):
                sendDataToJS(st.session_state.selected_experiment0_tagger,  layout[0])
    with c2:
        st.selectbox(
            "choose experiment", results, 
            key=f'selected_experiment_dropdown_1_tagger',
            index = name_to_index[st.session_state[f'selected_experiment1_tagger']] if f'selected_experiment1_tagger' in st.session_state else None,
            on_change=select_experiment
        )
        if f"selected_experiment1_tagger" in st.session_state:
            with st.spinner('Loading component...'):
                sendDataToJS(st.session_state["selected_experiment1_tagger"], layout[1], 'flash_viewer_grid_1')

else:
    ### for only single experiment on one view
    st.selectbox(
        "choose experiment", results, 
        key="selected_experiment_dropdown_tagger", 
        index=name_to_index[st.session_state.selected_experiment0_tagger] if 'selected_experiment0_tagger' in st.session_state else None,
        on_change=select_experiment
    )

    if 'selected_experiment0_tagger' in st.session_state:
        with st.spinner('Loading component...'):
            sendDataToJS(st.session_state.selected_experiment0_tagger, layout[0])

    ### for multiple experiments on one view
    if len(layout) > 1:

        for exp_index, exp_layout in enumerate(layout):
            if exp_index == 0: continue  # skip the first experiment

            st.divider() # horizontal line

            st.selectbox(
                "choose experiment", results, 
                key=f'selected_experiment_dropdown_{exp_index}_tagger',
                index = name_to_index[st.session_state[f'selected_experiment{exp_index}_tagger']] if f'selected_experiment{exp_index}_tagger' in st.session_state else None,
                on_change=select_experiment
            )

            # if #experiment input files are less than #layouts, all the pre-selection will be the first experiment
            if f"selected_experiment{exp_index}_tagger" in st.session_state:
                with st.spinner('Loading component...'):
                    sendDataToJS(st.session_state["selected_experiment%d_tagger" % exp_index], layout[exp_index], 'flash_viewer_grid_%d' % exp_index)

save_params(params)