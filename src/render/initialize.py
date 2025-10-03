import polars as pl

from src.render.components import (
    PlotlyHeatmap, PlotlyLineplot, PlotlyLineplotTagger, Plotly3Dplot, 
    Tabulator, SequenceView, InternalFragmentMap, FlashViewerComponent, 
    FDRPlotly, FLASHQuant, Chromatogram
)
from src.render.compression import compute_compression_levels

def initialize_data(comp_name, selected_data, file_manager, tool):

    data_to_send = {}
    additional_data = {'dataset' : selected_data}

    if comp_name == 'ms1_deconv_heat_map':

        # Fetch full dataset
        data_full = file_manager.get_results(
            selected_data,  ['ms1_deconv_heatmap'], use_polars=True
        )['ms1_deconv_heatmap']

        # Fetch all caches
        cached_compression_levels = []
        for size in compute_compression_levels(20000, data_full.select(pl.len()).collect(engine="streaming").item()):
            cached_compression_levels.append(
                file_manager.get_results(
                    selected_data,  [f'ms1_deconv_heatmap_{size}'], use_polars=True
                )[f'ms1_deconv_heatmap_{size}']
            )
        cached_compression_levels.append(data_full)

        # Get smallest compression level
        data_to_send['deconv_heatmap_df'] = cached_compression_levels[0]

        additional_data['deconv_heatmap_df'] = cached_compression_levels
        component_arguments = PlotlyHeatmap(title="Deconvolved MS1 Heatmap")
    elif comp_name == 'ms2_deconv_heat_map':

        # Fetch full dataset
        data_full = file_manager.get_results(
            selected_data,  ['ms2_deconv_heatmap'], use_polars=True
        )['ms2_deconv_heatmap']

        # Fetch all caches
        cached_compression_levels = []
        for size in compute_compression_levels(20000, data_full.select(pl.len()).collect(engine="streaming").item()):
            cached_compression_levels.append(
                file_manager.get_results(
                    selected_data,  [f'ms2_deconv_heatmap_{size}'], use_polars=True
                )[f'ms2_deconv_heatmap_{size}']
            )
        cached_compression_levels.append(data_full)

        # Get smallest compression level
        data_to_send['deconv_heatmap_df'] = cached_compression_levels[0]

        additional_data['deconv_heatmap_df'] = cached_compression_levels
        component_arguments = PlotlyHeatmap(title="Deconvolved MS2 Heatmap")

    elif comp_name == 'ms1_raw_heatmap':

        # Fetch full dataset
        data_full = file_manager.get_results(
            selected_data,  ['ms1_raw_heatmap'], use_polars=True
        )['ms1_raw_heatmap']

        # Fetch all caches
        cached_compression_levels = []
        for size in compute_compression_levels(20000, data_full.select(pl.len()).collect(engine="streaming").item()):
            cached_compression_levels.append(
                file_manager.get_results(
                    selected_data,  [f'ms1_raw_heatmap_{size}'], use_polars=True
                )[f'ms1_raw_heatmap_{size}']
            )
        cached_compression_levels.append(data_full)

        # Get smallest compression level
        data_to_send['raw_heatmap_df'] = cached_compression_levels[0]

        additional_data['raw_heatmap_df'] = cached_compression_levels

        component_arguments = PlotlyHeatmap(title="Raw MS1 Heatmap")
    elif comp_name == 'ms2_raw_heatmap':

        # Fetch full dataset
        data_full = file_manager.get_results(
            selected_data,  ['ms2_raw_heatmap'], use_polars=True
        )['ms2_raw_heatmap']

        # Fetch all caches
        cached_compression_levels = []
        for size in compute_compression_levels(20000, data_full.select(pl.len()).collect(engine="streaming").item()):
            cached_compression_levels.append(
                file_manager.get_results(
                    selected_data,  [f'ms2_raw_heatmap_{size}'], use_polars=True
                )[f'ms2_raw_heatmap_{size}']
            )
        cached_compression_levels.append(data_full)

        # Get smallest compression level
        data_to_send['raw_heatmap_df'] = cached_compression_levels[0]

        additional_data['raw_heatmap_df'] = cached_compression_levels

        component_arguments = PlotlyHeatmap(title="Raw MS2 Heatmap")
    elif comp_name == 'scan_table':
        data = file_manager.get_results(selected_data, ['scan_table'])
        data_to_send['per_scan_data'] = data['scan_table']
        component_arguments = Tabulator('ScanTable')
    elif comp_name == 'deconv_spectrum':
        data = file_manager.get_results(selected_data, ['deconv_spectrum'])
        data_to_send['per_scan_data'] = data['deconv_spectrum']
        component_arguments = PlotlyLineplot(title="Deconvolved Spectrum")
    elif comp_name == 'combined_spectrum':
        data = file_manager.get_results(selected_data, ['combined_spectrum'])
        data_to_send['per_scan_data'] = data['combined_spectrum']
        component_arguments = PlotlyLineplotTagger(title="Augmented Deconvolved Spectrum")
    elif comp_name == 'anno_spectrum':
        data = file_manager.get_results(selected_data,  ['combined_spectrum'])
        data_to_send['per_scan_data'] = data['combined_spectrum']
        component_arguments = PlotlyLineplot(title="Annotated Spectrum")
    elif comp_name == 'mass_table':
        data = file_manager.get_results(selected_data,  ['mass_table'])
        data_to_send['per_scan_data'] = data['mass_table']
        component_arguments = Tabulator('MassTable')
    elif comp_name == '3D_SN_plot':
        data = file_manager.get_results(selected_data,  ['threedim_SN_plot'], use_pyarrow=True)
        data_to_send['per_scan_data'] = data['threedim_SN_plot']
        component_arguments = Plotly3Dplot(title="Precursor Signals")
    elif comp_name == 'sequence_view':
        data = file_manager.get_results(selected_data,  ['sequence_view'])
        data_to_send['per_scan_data'] = data['sequence_view']
        if tool == 'flashtnt':
            data = file_manager.get_results(selected_data,  ['sequence_data'])
            data_to_send['sequence_data'] = data['sequence_data']
            data = file_manager.get_results(selected_data,  ['settings'])
            data_to_send['settings'] = data['settings']
        component_arguments = SequenceView(title='Sequence View')
    elif comp_name == 'internal_fragment_map':
        data = file_manager.get_results(selected_data,  ['sequence_view'])
        data_to_send['per_scan_data'] = data['sequence_view']
        if tool == 'flashtnt':
            data = file_manager.get_results(selected_data,  ['sequence_data'])
            data_to_send['sequence_data'] = data['sequence_data']
            data = file_manager.get_results(selected_data,  ['internal_fragment_data'])
            data_to_send['internal_fragment_data'] = data['internal_fragment_data']
        component_arguments = InternalFragmentMap(title="Internal Fragment Map")    
    elif comp_name == 'fdr_plot':
        data = file_manager.get_results(selected_data,  ['density_target'])
        data_to_send['density_target'] = data['density_target']
        data = file_manager.get_results(selected_data,  ['density_decoy'])
        data_to_send['density_decoy'] = data['density_decoy']
        component_arguments = FDRPlotly(title="FDR Plot")
    elif comp_name == 'id_fdr_plot':
        data = file_manager.get_results(selected_data,  ['density_id_target'])
        data_to_send['density_target'] = data['density_id_target']
        data = file_manager.get_results(selected_data,  ['density_id_decoy'])
        data_to_send['density_decoy'] = data['density_id_decoy']
        component_arguments = FDRPlotly(title="FDR Plot")
    elif comp_name == 'protein_table':
        # TODO: Unify lookup or remove in vue
        data = file_manager.get_results(selected_data,  ['scan_table'])
        data_to_send['per_scan_data'] = data['scan_table']
        data = file_manager.get_results(selected_data,  ['protein_dfs'])
        data_to_send['protein_table'] = data['protein_dfs']
        component_arguments = Tabulator('ProteinTable')
    elif comp_name == 'tag_table':
        data = file_manager.get_results(selected_data,  ['tag_dfs'])
        data_to_send['tag_table'] = data['tag_dfs']
        component_arguments = Tabulator('TagTable')
    elif comp_name == 'quant_visualization':
        data = file_manager.get_results(selected_data,  ['quant_dfs'])
        data_to_send['quant_data'] = data['quant_dfs']
        component_arguments = FLASHQuant()
    elif comp_name == 'tic_chromatogram':
        data = file_manager.get_results(selected_data,  ['tic'])
        data_to_send['tic'] = data['tic']
        component_arguments = Chromatogram()

    components = [[FlashViewerComponent(component_arguments)]]

    return data_to_send, components, additional_data
