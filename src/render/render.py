"""Viewer grid built from openms-insight components.

Every layout cell is one Insight component. The components of one grid share a
StateManager, so a click in one (a scan, a mass, a proteoform, a tag, a residue)
filters or highlights the others. Each component keeps a preprocessed cache next
to the workflow results (cache/insight/<dataset>/), built on first view and
rebuilt when the results it was built from change.
"""

import hashlib
import json
import shutil
import threading
from pathlib import Path

import polars as pl
import streamlit as st
from streamlit.runtime.scriptrunner import get_script_run_ctx

import openms_insight
from openms_insight import (
    DeconvolvedSpectrum, Heatmap, Plot3D, ProteinSequenceView, SeriesPlot,
    StateManager, Table, defer_reruns,
)
from openms_insight.core.cache import CacheMissError

from src.render import insight_data
from src.render.sequence import getFragmentDataFromSeq


def _code_version():
    """Changes whenever the configuration of a component can: with this module,
    the frame builders or openms-insight. Part of every cache id, so caches
    are rebuilt after an update."""
    digest = hashlib.md5(openms_insight.__version__.encode())
    for module in (__file__, insight_data.__file__):
        try:
            digest.update(Path(module).read_bytes())
        except OSError:
            # Frozen builds may not ship the sources: the version alone
            digest.update(module.encode())
    return digest.hexdigest()


CODE_VERSION = _code_version()


def _fixed(precision=4):
    """FLASHViewer's number format: four decimals for values longer than four
    characters, short values as they are."""
    return {'formatter': 'fixed', 'formatterParams': {'precision': precision, 'minLength': 4}}


# -1 marks a value FLASHTnT did not determine; tags it could not place on the
# proteoform have no positions
_DASH_FOR_MISSING = {'formatter': 'fixed', 'formatterParams': {'placeholders': {'-1': '-', 'NaN': ''}, 'minLength': 1000}}
_BLANK_FOR_MISSING = {'formatter': 'fixed', 'formatterParams': {'placeholders': {'NaN': ''}, 'minLength': 1000}}

SCAN_COLUMNS = [
    {'title': 'Index', 'field': 'index', 'sorter': 'number',
     'headerTooltip': 'The sequential index of the spectrum in the dataset.'},
    {'title': 'Scan Number', 'field': 'Scan', 'sorter': 'number',
     'headerTooltip': 'The identifier of the mass spectrometry scan.'},
    {'title': 'MS Level', 'field': 'MSLevel', 'sorter': 'number',
     'headerTooltip': 'The level of mass spectrometry analysis (e.g., MS1 or MS2).'},
    {'title': 'Retention time', 'field': 'RT', 'sorter': 'number', **_fixed(),
     'headerTooltip': 'The time at which the spectrum was detected during the chromatographic separation in seconds.'},
    {'title': 'Precursor Mass', 'field': 'PrecursorMass', 'sorter': 'number', **_fixed(),
     'headerTooltip': 'The mass of the precursor ion selected for fragmentation in Daltons.'},
    {'title': '#Masses', 'field': '#Masses', 'sorter': 'number',
     'headerTooltip': 'The number of detected masses in the spectrum.'},
]

MASS_COLUMNS = [
    {'title': 'Index', 'field': 'mass_idx', 'sorter': 'number',
     'headerTooltip': 'The sequential index of the mass entry in the dataset.'},
    {'title': 'Monoisotopic mass', 'field': 'MonoMass', 'sorter': 'number', **_fixed(),
     'headerTooltip': 'The monoisotopic mass of the detected ion in Daltons.'},
    {'title': 'Sum intensity', 'field': 'SumIntensity', 'sorter': 'number', **_fixed(),
     'headerTooltip': 'The total intensity of the detected mass across all isotopic peaks and charges.'},
    {'title': 'Min charge', 'field': 'MinCharges', 'sorter': 'number',
     'headerTooltip': 'The minimum charge state detected for the mass.'},
    {'title': 'Max charge', 'field': 'MaxCharges', 'sorter': 'number',
     'headerTooltip': 'The maximum charge state detected for the mass.'},
    {'title': 'Min isotope', 'field': 'MinIsotopes', 'sorter': 'number',
     'headerTooltip': 'The smallest observed isotopic shift, expressed as a multiple of the average isotopic mass difference at 55kDA.'},
    {'title': 'Max isotope', 'field': 'MaxIsotopes', 'sorter': 'number',
     'headerTooltip': 'The largest observed isotopic shift, expressed as a multiple of the average isotopic mass difference at 55kDA.'},
    {'title': 'Cosine score', 'field': 'CosineScore', 'sorter': 'number', **_fixed(),
     'headerTooltip': 'The cosine similarity score comparing the observed and theoretical isotopic patterns.'},
    {'title': 'SNR', 'field': 'SNR', 'sorter': 'number', **_fixed(),
     'headerTooltip': 'The signal-to-noise ratio for the detected mass.'},
    {'title': 'QScore', 'field': 'QScore', 'sorter': 'number', **_fixed(),
     'headerTooltip': 'The quality score indicating the confidence of the mass detection (higher is better).'},
]

PROTEIN_COLUMNS = [
    {'title': 'Scan No.', 'field': 'Scan', 'sorter': 'number',
     'headerTooltip': 'The identifier of the mass spectrometry scan associated with the identified proteoform.'},
    {'title': 'Accession', 'field': 'accession',
     'headerTooltip': 'The unique identifier for the protein in the reference database.'},
    {'title': 'Description', 'field': 'description', 'responsive': 10},
    {'title': 'Length', 'field': 'length', 'responsive': 6, 'sorter': 'number',
     'headerTooltip': 'The total number of amino acids in the matched protein.'},
    {'title': 'Mass', 'field': 'ProteoformMass', 'responsive': 8, 'sorter': 'number', **_DASH_FOR_MISSING,
     'headerTooltip': 'The calculated mass of the proteoform in Daltons.'},
    {'title': 'No. of Matched Fragments', 'field': 'MatchingFragments', 'sorter': 'number',
     'headerTooltip': 'The number of fragment ions that match the protein sequence.'},
    {'title': 'No. of Modifications', 'field': 'ModCount', 'sorter': 'number',
     'headerTooltip': 'The number of modifications identified in the protein.'},
    {'title': 'No. of Tags', 'field': 'TagCount', 'sorter': 'number',
     'headerTooltip': 'The number of sequence tags associated with the proteoform match.'},
    {'title': 'Score', 'field': 'Score', 'sorter': 'number',
     'headerTooltip': 'A score indicating the confidence of the protein match (higher is better).'},
    {'title': 'Q-Value (Proteoform Level)', 'field': 'ProteoformLevelQvalue', 'sorter': 'number',
     **_DASH_FOR_MISSING,
     'headerTooltip': 'The confidence value of the protein match at the proteoform level.'},
]

TAG_COLUMNS = [
    {'title': 'Scan Number', 'field': 'Scan', 'sorter': 'number',
     'headerTooltip': 'The identifier of the mass spectrometry scan containing the sequence tag.'},
    {'title': 'Start Position', 'field': 'StartPos', 'sorter': 'number', **_BLANK_FOR_MISSING,
     'headerTooltip': 'The position in the protein sequence where the sequence tag begins.'},
    {'title': 'End Position', 'field': 'EndPos', 'sorter': 'number', **_BLANK_FOR_MISSING,
     'headerTooltip': 'The position in the protein sequence where the sequence tag ends.'},
    # FLASHTnT also reports long de novo tags (hundreds of residues)
    {'title': 'Sequence', 'field': 'TagSequence', 'sorter': 'string', 'maxWidth': 240,
     'tooltip': True, 'headerTooltip': 'The amino acid sequence of the identified tag.'},
    {'title': 'Length', 'field': 'Length', 'sorter': 'number',
     'headerTooltip': 'The number of amino acids in the sequence tag.'},
    {'title': 'Tag Score', 'field': 'Score', 'sorter': 'number',
     'headerTooltip': 'A score indicating the confidence of the sequence tag identification (higher is better).'},
    {'title': 'N mass', 'field': 'Nmass', 'sorter': 'number', **_DASH_FOR_MISSING,
     'headerTooltip': 'The N-terminal mass offset from the start of the sequence tag in Daltons.'},
    {'title': 'C mass', 'field': 'Cmass', 'sorter': 'number', **_DASH_FOR_MISSING,
     'headerTooltip': 'The C-terminal mass offset from the end of the sequence tag in Daltons.'},
    {'title': 'Δ mass', 'field': 'DeltaMass', 'sorter': 'number', **_BLANK_FOR_MISSING,
     'headerTooltip': 'Delta mass is the difference between the tag flanking mass and the (partial) proteoform mass, from its terminal to the tag boundary.'},
]

QUANT_COLUMNS = [
    {'title': 'Index', 'field': 'FeatureGroupIndex', 'sorter': 'number'},
    {'title': 'Monoisotopic Mass', 'field': 'MonoisotopicMass', 'sorter': 'number'},
    {'title': 'Average Mass', 'field': 'AverageMass', 'sorter': 'number'},
    {'title': 'Start Retention Time (FWHM)', 'field': 'StartRetentionTime(FWHM)', 'sorter': 'number'},
    {'title': 'End Retention Time (FWHM)', 'field': 'EndRetentionTime(FWHM)', 'sorter': 'number'},
    {'title': 'Feature Group Quantity', 'field': 'FeatureGroupQuantity', 'sorter': 'number'},
    # Twice, as in FLASHViewer (the second one is the first in the collapsed details)
    {'title': 'Feature Group Quantity', 'field': 'FeatureGroupQuantity', 'sorter': 'number'},
    {'title': 'Min Charge', 'field': 'MinCharge', 'sorter': 'number'},
    {'title': 'Max Charge', 'field': 'MaxCharge', 'sorter': 'number'},
    {'title': 'Most Abundant Charge', 'field': 'MostAbundantFeatureCharge', 'sorter': 'number'},
    {'title': 'Isotope Cosine Score', 'field': 'IsotopeCosineScore', 'sorter': 'number'},
]

# Columns of the FLASHDeconv spectrum caches (see src/parse/masstable.py)
SPECTRUM_COLUMNS = dict(
    mass_column='MonoMass', mass_intensity_column='SumIntensity',
    peak_mz_column='MonoMass_Anno', peak_intensity_column='SumIntensity_Anno',
    signal_column='SignalPeaks',
)
SPECTRUM_TITLES = {
    'deconvolved': 'Deconvolved Spectrum',
    'annotated': 'Annotated Spectrum',
    'ladder': 'Augmented Deconvolved Spectrum',
}
SEQUENCE_VIEW_OPTIONS = dict(
    mass_column='MonoMass', precursor_mass_column='PrecursorMass',
    annotation_title='FLASHAnnotator', svg_filename='FLASHViewer-sequence',
)

HEATMAPS = {
    # component: (result, title, y label, click sets the mass too)
    'ms1_raw_heatmap': ('ms1_raw_heatmap', 'Raw MS1 Heatmap', 'm/z', False),
    'ms2_raw_heatmap': ('ms2_raw_heatmap', 'Raw MS2 Heatmap', 'm/z', False),
    'ms1_deconv_heat_map': ('ms1_deconv_heatmap', 'Deconvolved MS1 Heatmap', 'Monoisotopic Mass', True),
    'ms2_deconv_heat_map': ('ms2_deconv_heatmap', 'Deconvolved MS2 Heatmap', 'Monoisotopic Mass', True),
}


# Clicking the highlighted row of a table again only removes the highlight; the
# selection stays (FLASHViewer's tables behaved so)
RECLICK = 'unhighlight'
# A selected row is scrolled to the top of its table, a clicked row only when cut
# off (not in FLASHQuant, whose table did not follow its selection)
SCROLL = 'top'
# FLASHViewer's tables: as tall as their rows, up to 320 px below a 56 px title bar
TABLE_HEIGHT = 376


class MissingResult(Exception):
    """A result the component is drawn from was not written for this dataset."""


def _stamp(file_manager, dataset_id, sources, extra):
    """Hash of what a component cache is built from: the source files as
    written (path, size, mtime) and the configuration in `extra`."""
    parts = [CODE_VERSION, extra]
    for name_tag in sources:
        path = insight_data.result_path(file_manager, dataset_id, name_tag)
        if path is None:
            parts.append([name_tag, None])
            continue
        stat = path.stat()
        parts.append([name_tag, path.name, stat.st_size, stat.st_mtime_ns])
    return hashlib.md5(json.dumps(parts, sort_keys=True, default=str).encode()).hexdigest()[:12]


def _component(cls, name, file_manager, dataset_id, sources, build, extra=None, required=None):
    """The Insight component `name` of a dataset: loaded from its cache, or
    built with the keyword arguments `build()` returns.

    `sources` are the results the cache is built from; `required` (default:
    all of them) must exist.
    """
    for name_tag in sources if required is None else required:
        if insight_data.result_path(file_manager, dataset_id, name_tag) is None:
            raise MissingResult(name_tag)
    cache_path = Path(file_manager.cache_path, 'insight', dataset_id)
    cache_id = f'{name}-{_stamp(file_manager, dataset_id, sources, extra)}'
    try:
        return cls(cache_id=cache_id, cache_path=str(cache_path))
    except CacheMissError:
        pass
    # Caches built from earlier results or configurations are not used again.
    # Not for caches that depend on the selection (`extra`, e.g. a typed
    # sequence): two grids on a page may each need their own.
    if extra is None and cache_path.exists():
        for old in cache_path.glob(f'{name}-*'):
            if old.name != cache_id:
                shutil.rmtree(old, ignore_errors=True)
    return cls(cache_id=cache_id, cache_path=str(cache_path), **build())


def _heatmap(comp_name, file_manager, dataset_id):
    result, title, y_label, sets_mass = HEATMAPS[comp_name]
    interactivity = {'scan': 'scan_idx'}
    if sets_mass:
        interactivity['mass'] = 'mass_idx'
    if insight_data.result_path(file_manager, dataset_id, result) is None:
        raise MissingResult(result)
    # The levels the parse downsampled; without them Insight builds its own
    levels = insight_data.heatmap_levels(file_manager, dataset_id, result)
    return _component(
        Heatmap, comp_name, file_manager, dataset_id, [result, *(levels or [])],
        lambda: dict(
            data=insight_data.scan(file_manager, dataset_id, result),
            levels=None if levels is None else [
                insight_data.scan(file_manager, dataset_id, name) for name in levels
            ],
            x_column='rt', y_column='mass', intensity_column='intensity',
            interactivity=interactivity, zoom_identifier=f'{comp_name}_zoom',
            title=title, x_label='Retention Time', y_label=y_label,
            colorscale='Portland', log_scale=True,
            min_points=insight_data.HEATMAP_MIN_POINTS,
            # As in FLASHViewer, a zoomed view shows every point of the first
            # level with 20,000 points in range
            downsample_zoomed=False,
            # FLASHViewer's export name; its wheel scrolled the page
            svg_filename='FLASHViewer-heatmap', scroll_zoom=False,
        ),
    )


def _fdr_plot(comp_name, file_manager, dataset_id):
    if comp_name == 'id_fdr_plot':
        target, decoy = 'density_id_target', 'density_id_decoy'
    else:
        target, decoy = 'density_target', 'density_decoy'
    return _component(
        SeriesPlot, comp_name, file_manager, dataset_id, [target, decoy],
        lambda: dict(
            data=insight_data.fdr_densities(file_manager, dataset_id, target, decoy),
            x_column='x', y_column='y', series_column='kind',
            series_names={'target': 'Target QScores', 'decoy': 'Decoy QScores'},
            series_colors={'target': 'green', 'decoy': 'red'},
            mode='lines+markers', title='FDR Plot', x_label='QScore', y_label='Density',
            fixed_y_range=True, svg_filename='FDR-plot',
            # FLASHViewer's FDR plot: Arial in black, densities unrounded
            font={'family': 'Arial', 'color': 'black'}, float_precision='float64',
        ),
        required=[target],
    )


def _deconv_sequence(file_manager, state_manager):
    """Sequence data for the FLASHDeconv sequence view: the sequence the user
    typed into the view, else the one saved on the Sequence Input page."""
    stored = None
    if file_manager.result_exists('sequence', 'sequence'):
        stored = file_manager.get_results('sequence', 'sequence')['sequence']
    if stored is None:
        return None
    sequence = state_manager.get_selection('sequence') or stored['input_sequence']
    return _sequence_data(
        sequence, stored.get('fixed_mod_cysteine'), stored.get('fixed_mod_methionine')
    )


@st.cache_data(max_entries=4, show_spinner=False)
def _sequence_data(sequence, cysteine, methionine):
    data = getFragmentDataFromSeq(
        sequence,
        fixed_mods={'fixed_mod_cysteine': cysteine, 'fixed_mod_methionine': methionine},
    )
    # Plain floats: the sequence is part of the component configuration
    return json.loads(json.dumps(data, default=float))


def _selected_tag(file_manager, dataset_id, state_manager):
    """The row of the selected sequence tag, or None."""
    tag_index = state_manager.get_selection('tag')
    if tag_index is None:
        return None
    rows = (
        insight_data.tags(file_manager, dataset_id)
        .filter(pl.col('TagIndex') == tag_index)
        .select('StartPos', 'EndPos', 'TagSequence', 'mzs')
        .head(1)
        .collect()
    )
    return rows.row(0, named=True) if rows.height else None


def _build(comp_name, file_manager, dataset_id, tool, state_manager):
    """Insight component of a layout cell, and keyword arguments for rendering it."""
    fm, ds = file_manager, dataset_id

    if comp_name in HEATMAPS:
        return _heatmap(comp_name, fm, ds), {'height': 450}

    if comp_name == 'scan_table':
        return _component(
            Table, comp_name, fm, ds, ['scan_table'],
            lambda: dict(
                data=insight_data.scan_table(fm, ds),
                column_definitions=SCAN_COLUMNS, title='Scan Table',
                index_field='index', go_to_fields=['index', 'Scan'],
                interactivity={'scan': 'index'}, on_select={'mass': 0},
                reclick=RECLICK, scroll_position=SCROLL, column_menu=False, reset_sort=True, pagination=False, float_precision='float64',
            ),
        ), {'height': TABLE_HEIGHT}

    if comp_name == 'mass_table':
        return _component(
            Table, comp_name, fm, ds, ['mass_table'],
            lambda: dict(
                data=insight_data.mass_table(fm, ds),
                column_definitions=MASS_COLUMNS, title='Mass Table',
                index_field='mass_idx', go_to_fields=['mass_idx'],
                filters={'scan': 'scan_idx'}, interactivity={'mass': 'mass_idx'},
                # As in FLASHViewer, the first mass is highlighted but not
                # selected until a mass or another scan is clicked
                default_row_selects=False,
                reclick=RECLICK, scroll_position=SCROLL, column_menu=False, reset_sort=True, pagination=False, float_precision='float64',
            ),
        ), {'height': TABLE_HEIGHT}

    if comp_name in ('deconv_spectrum', 'anno_spectrum', 'combined_spectrum'):
        if comp_name == 'deconv_spectrum':
            result, view = 'deconv_spectrum', 'deconvolved'
        elif comp_name == 'anno_spectrum':
            result, view = 'combined_spectrum', 'annotated'
        else:
            result, view = 'combined_spectrum', 'ladder'
        # FLASHTnT spectra follow the spectrum of the selected proteoform
        scan_identifier = 'pscan' if tool == 'flashtnt' else 'scan'
        component = _component(
            DeconvolvedSpectrum, f'{comp_name}-{tool}', fm, ds, [result],
            lambda: dict(
                data=insight_data.spectra(fm, ds, result),
                filters={scan_identifier: 'index'}, view=view, mass_identifier='mass',
                title=SPECTRUM_TITLES[view], svg_filename='FLASHViewer-lineplot',
                drilldown_title='Augmented Annotated Spectrum',
                **SPECTRUM_COLUMNS,
                float_precision='float64',
            ),
        )
        if view == 'ladder':
            tag = _selected_tag(fm, ds, state_manager)
            if tag is None:
                component.set_ladder(None)
            else:
                masses = [float(m) for m in str(tag['mzs']).split(',') if m.strip() and float(m) != 0]
                aa = state_manager.get_selection('aa')
                selected = None
                if (
                    aa is not None and tag['StartPos'] is not None and tag['EndPos'] is not None
                    and tag['StartPos'] <= aa <= tag['EndPos']
                ):
                    selected = int(aa - tag['StartPos'])
                component.set_ladder(
                    masses, tag['TagSequence'], selected,
                    key=state_manager.get_selection('tag'),
                )
        return component, {}

    if comp_name == '3D_SN_plot':
        return _component(
            Plot3D, comp_name, fm, ds,
            ['threedim_SN_plot', 'scan_table', 'deconv_spectrum'],
            lambda: dict(
                data=insight_data.signal_peaks(fm, ds),
                filters={'scan': 'scan_idx', 'mass': 'mass_idx'},
                filter_defaults={'mass': -1},
                x_column='mass', y_column='charge', z_column='intensity', mode='stems',
                series_column='kind',
                series_colors={'Signal': '#3366CC', 'Noise': '#DC3912'},
                # Both in the legend, as in FLASHViewer, also when one is empty
                series_names={'Signal': 'Signal', 'Noise': 'Noise'},
                title='', selection_titles={'mass': 'Mass signals', 'scan': 'Precursor signals'},
                x_label='Mass', y_label='Charge', z_label='Intensity', y_dtick=1,
                camera_eye={'x': 2.5, 'y': 0, 'z': 0.2}, svg_filename='FLASHViewer-3d-plot',
            ),
        ), {'height': 800}

    if comp_name in ('fdr_plot', 'id_fdr_plot'):
        return _fdr_plot(comp_name, fm, ds), {'height': 400}

    if comp_name == 'protein_table':
        return _component(
            Table, comp_name, fm, ds, ['protein_dfs', 'scan_table'],
            lambda: dict(
                data=insight_data.proteins(fm, ds),
                column_definitions=PROTEIN_COLUMNS, title='Protein Table',
                index_field='index', go_to_fields=['Scan', 'accession'],
                initial_sort=[{'column': 'Score', 'dir': 'desc'}],
                toggle_filters=[{
                    'label': 'Best per spectrum', 'column': 'is_best', 'group_column': 'Scan',
                    # FLASHViewer's checkbox colour
                    'color': '#6200EE',
                }],
                interactivity={'protein': 'index', 'pscan': 'scan_idx'},
                on_select={'tag': None, 'aa': None, 'mass': None},
                reclick=RECLICK, scroll_position=SCROLL, column_menu=False, reset_sort=True, pagination=False, float_precision='float64',
            ),
        ), {'height': TABLE_HEIGHT}

    if comp_name == 'tag_table':
        return _component(
            Table, comp_name, fm, ds, ['tag_dfs', 'scan_table'],
            lambda: dict(
                data=insight_data.tags(fm, ds),
                column_definitions=TAG_COLUMNS, title='Tag Table',
                index_field='TagIndex', go_to_fields=['Scan', 'StartPos', 'EndPos', 'TagSequence'],
                initial_sort=[{'column': 'Score', 'dir': 'desc'}],
                filters={'pscan': 'scan_idx'}, range_filters={'aa': ['StartPos', 'EndPos']},
                interactivity={'tag': 'TagIndex'}, on_select={'mass': None},
                reclick=RECLICK, scroll_position=SCROLL, column_menu=False, reset_sort=True,
                reset_sort_on=['protein'], pagination=False, float_precision='float64',
            ),
        ), {'height': TABLE_HEIGHT}

    if comp_name == 'sequence_view':
        if tool == 'flashtnt':
            component = _component(
                ProteinSequenceView, f'{comp_name}-{tool}', fm, ds,
                ['sequence_data', 'protein_dfs', 'scan_table', 'sequence_view', 'settings'],
                lambda: dict(
                    data=insight_data.proteoform_sequences(fm, ds),
                    filters={'protein': 'proteoform_index'},
                    mass_identifier='mass', aa_identifier='aa',
                    settings=fm.get_results(ds, ['settings'])['settings'],
                    **SEQUENCE_VIEW_OPTIONS,
                ),
            )
            tag = _selected_tag(fm, ds, state_manager)
            # Tags FLASHTnT could not place on the proteoform have no positions
            if tag is None or tag['StartPos'] is None or tag['EndPos'] is None:
                component.set_tag(None)
            else:
                component.set_tag(int(tag['StartPos']), int(tag['EndPos']))
            return component, {}

        sequence = _deconv_sequence(fm, state_manager)
        if sequence is None:
            raise MissingResult('sequence')
        return _component(
            ProteinSequenceView, f'{comp_name}-{tool}', fm, ds, ['sequence_view'],
            lambda: dict(
                data=insight_data.sequence_view_spectra(fm, ds),
                filters={'scan': 'index'}, sequence=sequence,
                mass_identifier='mass', aa_identifier='aa', sequence_identifier='sequence',
                **SEQUENCE_VIEW_OPTIONS,
            ),
            extra=sequence,
        ), {}

    if comp_name == 'quant_visualization':
        return None, {}

    raise ValueError(f'Unknown viewer component {comp_name!r}')


def _render_quant(file_manager, dataset_id, state_manager, key):
    table = _component(
        Table, 'quant_table', file_manager, dataset_id, ['quant_dfs'],
        lambda: dict(
            data=insight_data.quant_table(file_manager, dataset_id),
            column_definitions=QUANT_COLUMNS, title='Feature groups',
            index_field='FeatureGroupIndex', interactivity={'feature_group': 'FeatureGroupIndex'},
            go_to_fields=[], reclick=RECLICK, column_menu=False, reset_sort=True, pagination=False, float_precision='float64',
        ),
    )
    traces = _component(
        Plot3D, 'quant_traces', file_manager, dataset_id, ['quant_dfs'],
        lambda: dict(
            data=insight_data.quant_traces(file_manager, dataset_id),
            filters={'feature_group': 'FeatureGroupIndex'},
            x_column='mz', y_column='rt', z_column='intensity', mode='lines',
            trace_column='trace', series_column='charge', series_name_template='Charge: {}',
            color='#3366CC', title='Feature group signals',
            x_label='m/z', y_label='retention time', z_label='intensity',
            # FLASHViewer's plot: plain title, no border, Plotly's PNG download
            title_bold=False, border=False, export_format='png',
        ),
    )
    # FLASHViewer drew the table and the plot in one 1184 px panel: the table 4 px
    # from its edges, the plot 16 px from the left, 90% wide, overlapping the
    # table's last 12 px
    panel = f'{key}_panel'
    with st.container(key=panel):
        st.markdown(
            f'''<style>
            .st-key-{panel} {{ gap: 0; }}
            .st-key-{panel} > div:nth-child(1) {{ display: none; }}
            .st-key-{panel} > div:nth-child(2) {{ margin: 4px 4px 0; width: calc(100% - 8px); }}
            .st-key-{panel} > div:nth-child(3) {{
                margin: -20px 0 16px 16px; width: calc((100% - 32px) * 0.9); z-index: 1;
            }}
            </style>''',
            unsafe_allow_html=True,
        )
        table(key=f'{key}_table', state_manager=state_manager, height=TABLE_HEIGHT)
        traces(key=f'{key}_traces', state_manager=state_manager, height=800)


def _render_cell(comp_name, file_manager, dataset_id, tool, state_manager, key):
    if comp_name == 'internal_fragment_map':
        st.info('The internal fragment map is not available yet.')
        return
    if comp_name == 'quant_visualization':
        _render_quant(file_manager, dataset_id, state_manager, key)
        return
    try:
        component, render_args = _build(comp_name, file_manager, dataset_id, tool, state_manager)
    except MissingResult as e:
        if str(e) == 'sequence':
            st.info('Enter a sequence on the Sequence Input page to show the sequence view.')
        else:
            st.info(f'This dataset has no data for this view ({e}).')
        return
    component(key=key, state_manager=state_manager, **render_args)


def _state_manager(tool, identifier, dataset_id):
    """Selections shared by the components of one grid; cleared when the grid
    switches to another dataset."""
    state_manager = StateManager(session_key=f'flash_insight_state_{tool}_{identifier}')
    dataset_key = f'flash_insight_dataset_{tool}_{identifier}'
    if st.session_state.get(dataset_key) != dataset_id:
        state_manager.clear()
        st.session_state[dataset_key] = dataset_id
    return state_manager


# One lock per Streamlit session. With runner.fastReruns (the default) Streamlit
# starts a new script thread while the previous one is still executing; two threads
# rendering the grid at the same time share the session's polars/pyarrow objects and
# crashed the server (see docs/white-screen-root-cause.md). Serialising render_grid
# per session removes that overlap: a runner that has been stopped releases the lock
# as soon as its next Streamlit call raises. Keyed by session id (not stored in
# st.session_state) so two concurrent first calls cannot create two different locks.
_GRID_LOCKS: dict[str, threading.Lock] = {}
_GRID_LOCKS_GUARD = threading.Lock()


def _grid_lock() -> threading.Lock:
    ctx = get_script_run_ctx()
    key = ctx.session_id if ctx is not None else '__no_session__'
    with _GRID_LOCKS_GUARD:
        return _GRID_LOCKS.setdefault(key, threading.Lock())


def render_grid(
    selected_data, layout_info_per_exp, file_manager, tool, identifier,
    grid_key='flash_viewer_grid'
):
    with _grid_lock():
        state_manager = _state_manager(tool, identifier, selected_data)
        # A click updates every view of the grid in one rerun, not one per view
        with defer_reruns():
            for row_index, row in enumerate(layout_info_per_exp):
                columns = st.columns(len(row))
                for col_index, (col, comp_name) in enumerate(zip(columns, row)):
                    with col:
                        _render_cell(
                            comp_name, file_manager, selected_data, tool, state_manager,
                            f'{grid_key}_{row_index}_{col_index}',
                        )
