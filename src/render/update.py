import pandas as pd
import streamlit as st
import pyarrow.dataset as ds

from src.render.compression import downsample_heatmap
from src.render.sequence import getFragmentDataFromSeq, getInternalFragmentDataFromSeq


# Ignore raw data for caching, too ressource intensive
hash_funcs = {pd.DataFrame : lambda x : 1}
@st.cache_data(max_entries=4, show_spinner=False, hash_funcs=hash_funcs)
def render_heatmap(full_data, selection, dataset_name, component_name):
    if (
        (selection['xRange'][0] < 0)
        and (selection['xRange'][1] < 0)
        and (selection['yRange'][0] < 0)
        and (selection['yRange'][1] < 0)
    ):
        return downsample_heatmap(full_data[0])
    else:
        for dataset in full_data:

            relevant_data = dataset.loc[
                (dataset['rt'] >= selection['xRange'][0])
                & (dataset['rt'] <= selection['xRange'][1])
                & (dataset['mass'] >= selection['yRange'][0])
                & (dataset['mass'] <= selection['yRange'][1])
            ]
            if len(relevant_data) >= 20000:
                break
        if len(relevant_data) <= 20000:
            return relevant_data
        downsampled = downsample_heatmap(relevant_data)
        return downsampled


@st.cache_data(max_entries=1, show_spinner=False)
def render_sequence_data(sequence):
    return getFragmentDataFromSeq(sequence)


@st.cache_data(max_entries=1, show_spinner=False)
def render_internal_fragment_data(sequence):
    return getInternalFragmentDataFromSeq(sequence)


def update_data(data, out_components, additional_data, tool):
    component = out_components[0][0]['componentArgs']['title']
    if (
        (component in ['Sequence View', 'Internal Fragment Map']) 
        and (tool != 'flashtnt')
    ):
        data['sequence_data'] = {
            0: render_sequence_data(st.session_state.input_sequence)
        }
    if (component == 'Internal Fragment Map') and (tool != 'flashtnt'):
        data['internal_fragment_data'] = {
            0: render_internal_fragment_data(st.session_state.input_sequence)
        }
    
    return data  


def filter_data(data, out_components, selection_store, additional_data, tool):
    data = data.copy()
    
    # Assumption: We are only dealing with one component
    component = out_components[0][0]['componentArgs']['title']

    # Filter data if possible
    if component in [
        'Annotated Spectrum', 'Deconvolved Spectrum', 
        'Augmented Deconvolved Spectrum', 
        'Mass Table', 'Sequence View', 'Internal Fragment Map'
    ]:
        if 'scanIndex' not in selection_store:
            data['per_scan_data'] = data['per_scan_data'].iloc[0:0,:]
        else:
            data['per_scan_data'] = data['per_scan_data'].iloc[selection_store['scanIndex']:selection_store['scanIndex']+1,:]
    elif component == 'Precursor Signals':
        scan_index = selection_store.get("scanIndex")
        mass_index = selection_store.get("massIndex")
        if scan_index is None:
            data['per_scan_data'] = data['per_scan_data'].to_table(filter=(ds.field("index") == -1)).slice(0, 0)
        else:
            filtered_table = data['per_scan_data'].to_table(filter=(ds.field("index") == scan_index))
            if mass_index is not None:
                df = filtered_table.to_pandas()
                df['SignalPeaks'] = df['SignalPeaks'].apply(lambda peaks: peaks[mass_index] if len(peaks) > mass_index else None)
                df['NoisyPeaks'] = df['NoisyPeaks'].apply(lambda peaks: peaks[mass_index] if len(peaks) > mass_index else None)
                filtered_table = df
            data['per_scan_data'] = filtered_table


    elif (component == 'Deconvolved MS1 Heatmap'):
        if 'heatmap_deconv' in selection_store:
            data['deconv_heatmap_df'] = render_heatmap(
                additional_data['deconv_heatmap_df'], 
                selection_store['heatmap_deconv'],
                additional_data['dataset'], component
            )
    elif (component == 'Raw MS1 Heatmap'):
        if 'heatmap_raw' in selection_store:
            data['raw_heatmap_df'] = render_heatmap(
                additional_data['raw_heatmap_df'], 
                selection_store['heatmap_raw'],
                additional_data['dataset'], component
            )

    if (
        (component in ['Internal Fragment Map', 'Sequence View']) 
        and (tool == 'flashtnt')
    ):
        if 'proteinIndex' not in selection_store:
            data['sequence_data'] = {}
        else:
            data['sequence_data'] = {
                selection_store['proteinIndex'] : data[
                    'sequence_data'
                ][selection_store['proteinIndex']]
            }

    if (component == 'Internal Fragment Map') and (tool == 'flashtnt'):
        if 'proteinIndex' not in selection_store:
            data['internal_fragment_data'] = {}
        else:
            data['internal_fragment_data'] = {
                selection_store['proteinIndex'] : data[
                    'internal_fragment_data'
                ][selection_store['proteinIndex']]
            }

    return data