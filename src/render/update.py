import pandas as pd
import polars as pl
import streamlit as st
import pyarrow.dataset as ds

from src.render.compression import downsample_heatmap
from src.workflow.FileManager import FileManager
from src.render.sequence import getFragmentDataFromSeq, getInternalFragmentDataFromSeq
from pathlib import Path


def get_sequence(selection_store):
    if 'sequenceOut' in selection_store:
        if len(selection_store['sequenceOut']) > 0:
            return selection_store['sequenceOut'], None, None
    # Setup cache access
    file_manager = FileManager(
        st.session_state["workspace"],
        Path(st.session_state['workspace'], 'flashdeconv', 'cache')
    )

    # Check if sequence has been set
    if not file_manager.result_exists('sequence', 'sequence'):
        return None
    # fetch sequence from cache
    sequence = file_manager.get_results('sequence', 'sequence')['sequence']

    return sequence['input_sequence'], sequence['fixed_mod_cysteine'], sequence['fixed_mod_methionine']


# Ignore raw data for caching, too ressource intensive
hash_funcs = {pl.LazyFrame : lambda x : 1}
@st.cache_data(max_entries=4, show_spinner=False, hash_funcs=hash_funcs)
def render_heatmap(full_data, selection, dataset_name, component_name):
    if (
        (selection['xRange'][0] < 0)
        and (selection['xRange'][1] < 0)
        and (selection['yRange'][0] < 0)
        and (selection['yRange'][1] < 0)
    ):
        return downsample_heatmap(full_data[0]).collect(engine="streaming")

    x0, x1 = selection['xRange']
    y0, y1 = selection['yRange']

    relevant_data = None
    est_count = 0
    for lf in full_data:
        filtered = lf.filter(
            (
                (pl.col("rt") >= x0) & (pl.col("rt") <= x1)
                & (pl.col("mass") >= y0) & (pl.col("mass") <= y1)
            )
        )
        est_count = (
            filtered
            .limit(20000)
            .select(pl.len().alias("n"))
            .collect(streaming=True)["n"][0]
        )

        relevant_data = filtered
        if est_count >= 20000:
            break

    if est_count <= 20000:
        # Small enough: return the filtered data eagerly
        return relevant_data.collect(engine="streaming")

    # Large: downsample lazily, then collect
    downsampled = downsample_heatmap(relevant_data)
    return downsampled.collect(engine="streaming")


@st.cache_data(max_entries=1, show_spinner=False)
def render_sequence_data(sequence):
    return getFragmentDataFromSeq(sequence)


@st.cache_data(max_entries=1, show_spinner=False)
def render_internal_fragment_data(sequence):
    return getInternalFragmentDataFromSeq(sequence)


def update_data(data, out_components, selection_store, additional_data, tool):
    component = out_components[0][0]['componentArgs']['title']
    if (
        (component in ['Sequence View', 'Internal Fragment Map']) 
        and (tool != 'flashtnt')
    ):
        data['sequence_data'] = {
            0: render_sequence_data(get_sequence(selection_store)[0])
        }
    if (component == 'Internal Fragment Map') and (tool != 'flashtnt'):
        data['internal_fragment_data'] = {
            0: render_internal_fragment_data(get_sequence(selection_store)[0])
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

    elif (component in ['Deconvolved MS1 Heatmap', 'Deconvolved MS2 Heatmap']):
        selection = 'heatmap_deconv' if '1' in component else 'heatmap_deconv2'
        if selection not in selection_store:
            selected_data = {
                'xRange' : [-1, -1],
                'yRange' : [-1, -1]
            }
        else:
            selected_data = selection_store[selection]
        data['deconv_heatmap_df'] = render_heatmap(
            additional_data['deconv_heatmap_df'], 
            selected_data,
            additional_data['dataset'], component
        )
    elif (component == ['Raw MS1 Heatmap', 'Raw MS2 Heatmap']):
        selection = 'heatmap_raw' if '1' in component else 'heatmap_raw2'
        if selection not in selection_store:
            selected_data = {
                'xRange' : [-1, -1],
                'yRange' : [-1, -1]
            }
        else:
            selected_data = selection_store[selection]
        data['raw_heatmap_df'] = render_heatmap(
            additional_data['raw_heatmap_df'], 
            selected_data,
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