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
from functools import lru_cache
from pathlib import Path

import polars as pl
import streamlit as st
from streamlit.runtime.scriptrunner import get_script_run_ctx

import openms_insight
from openms_insight import (
    DeconvolvedSpectrum, Heatmap, LinePlot, MirrorPlot, Plot3D, ProteinSequenceView,
    ProteinStructureView, SeriesPlot, StateManager, Table, defer_reruns,
)
from openms_insight.core.cache import CacheMissError

from src.render import (
    go_terms, hierarchy_data, insight_data, modifications, prsm_data, prsm_overview, qc_data, validation_data,
)
from src.render.fragmentation import mzml_scan_methods
from src.quant import isobaric
from src.render.layout_editor.catalog import DECONV_COMPONENTS, ISOBARIC_COMPONENTS, TNT_COMPONENTS
from src.render.sequence import getFragmentDataFromSeq
from src.render.structure import (
    StructureStore, download_enabled, get_structure, structure_dir, structure_key,
)
from src.render.structure_files import structure_action_handler, upload_errors


def _code_version():
    """Changes whenever the configuration of a component can: with this module,
    the frame builders or openms-insight. Part of every cache id, so caches
    are rebuilt after an update."""
    digest = hashlib.md5(openms_insight.__version__.encode())
    for module in (
        __file__, insight_data.__file__, validation_data.__file__, qc_data.__file__,
        hierarchy_data.__file__, prsm_data.__file__, go_terms.__file__,
        str(Path(__file__).with_name('quant_cells.py')),
    ):
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

# Masses in Da; FLASHTnT's undetermined values arrive as nulls and show "-"
_MASS = {'formatter': 'fixed', 'formatterParams': {'precision': 4, 'minLength': 4, 'placeholders': {'NaN': '-'}}}
# A column that starts hidden; the column chooser shows it
_HIDDEN = {'visible': False}

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

# The protein table shows the essential columns; the rest start hidden (column
# chooser) and all of them appear in a row's overview (the button at the
# start of each row), in this order.
PROTEIN_COLUMNS = [
    {'title': 'Scan No.', 'field': 'Scan', 'sorter': 'number',
     'headerTooltip': 'The identifier of the mass spectrometry scan associated with the identified proteoform.'},
    {'title': 'Accession', 'field': 'accession',
     'headerTooltip': 'The unique identifier for the protein in the reference database.'},
    {'title': 'Description', 'field': 'description', **_HIDDEN},
    {'title': 'Proteoform', 'field': 'Proforma', 'tooltip': True, **_HIDDEN,
     'headerTooltip': 'The proteoform in ProForma notation.'},
    {'title': 'Residues', 'field': 'residues', **_HIDDEN,
     'headerTooltip': 'First and last residue of the proteoform in the protein and the protein length; "?" for a terminus that was not determined.'},
    {'title': 'Length', 'field': 'length', 'sorter': 'number', **_HIDDEN,
     'headerTooltip': 'The total number of amino acids in the matched protein.'},
    {'title': 'Characterization', 'field': 'Characterization', 'sorter': 'string',
     'formatter': 'badge', 'formatterParams': {'colorMap': {
         'Characterized': '#2E7D32', 'Partial': '#F9A825', 'Ambiguous': '#C62828',
     }},
     'headerTooltip': (
         'How completely the fragments pin down the proteoform (in the spirit of '
         "ProSight's C-score colors). Characterized: both termini determined, every "
         'mass shift localized to one residue and at least 25% of the bonds explained '
         'by fragments. Ambiguous: less than 10% of the bonds explained. Partial: '
         'anything in between.'
     )},
    {'title': 'Retention time (s)', 'field': 'RetentionTime', 'sorter': 'number', **_fixed(), **_HIDDEN},
    {'title': 'Precursor m/z', 'field': 'precursor_mz', 'sorter': 'number', **_HIDDEN, **_MASS},
    {'title': 'Precursor charge', 'field': 'precursor_charge', 'sorter': 'number', **_HIDDEN},
    {'title': 'Precursor mass', 'field': 'precursor_mass', 'sorter': 'number', **_HIDDEN, **_MASS,
     'headerTooltip': 'Monoisotopic precursor mass in Daltons.'},
    {'title': 'Mass', 'field': 'ProteoformMass', 'sorter': 'number', **_DASH_FOR_MISSING, **_HIDDEN,
     'headerTooltip': 'The calculated mass of the proteoform in Daltons.'},
    {'title': 'Δ mass', 'field': 'delta_mass', 'sorter': 'number', **_HIDDEN, **_MASS,
     'headerTooltip': 'Precursor mass minus proteoform mass, in Daltons.'},
    {'title': 'No. of Matched Fragments', 'field': 'MatchingFragments', 'sorter': 'number',
     'headerTooltip': 'The number of fragment ions that match the protein sequence.'},
    {'title': 'Coverage (%)', 'field': 'Coverage(%)', 'sorter': 'number', **_HIDDEN,
     'headerTooltip': 'Residues of the protein covered by the proteoform, in %.'},
    {'title': 'Cleavages (%)', 'field': 'Cleavages', 'sorter': 'number', **_HIDDEN,
     'headerTooltip': 'Bonds of the proteoform explained by matched terminal fragments, in %.'},
    {'title': 'No. of Modifications', 'field': 'ModCount', 'sorter': 'number', **_HIDDEN,
     'headerTooltip': 'The number of modifications identified in the protein.'},
    {'title': 'Localized Mods', 'field': 'LocalizedMods', 'sorter': 'string', **_HIDDEN,
     'headerTooltip': 'Modifications (mass shifts) localized to a single residue, of all modifications.'},
    {'title': 'Modification sites', 'field': 'modification_sites', 'tooltip': True, **_HIDDEN,
     'headerTooltip': 'Mass shift and residue of each modification, or the residue range of one that was not localized.'},
    {'title': 'No. of Tags', 'field': 'TagCount', 'sorter': 'number', **_HIDDEN,
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
    # FLASHViewer's dialog button colour
    accent_color='#6200EE',
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
# Selection the sequence view's dialog writes its fixed modifications to
FIXED_MODS_IDENTIFIER = 'fixed_mods'

RECLICK = 'unhighlight'
# A selected row is scrolled to the top of its table, a clicked row only when cut
# off
SCROLL = 'top'
# Identification levels: proteins, the proteoforms of a protein,
# the PrSMs of a proteoform (src/render/hierarchy_data.py)
_QVALUE = {'formatter': 'fixed', 'formatterParams': {'precision': 4, 'minLength': 0, 'placeholders': {'NaN': '-'}}}
# hierarchy_data turns FLASHTnT's -1 / 0 into nulls; show them as "-"
_POSITION = {'formatter': 'fixed', 'formatterParams': {'minLength': 1000, 'placeholders': {'NaN': '-'}}}
# The three tables are meant to sit side by side in one layout row, so each
# shows only a few narrow columns; the rest start hidden (column chooser).
PROTEIN_GROUP_COLUMNS = [
    {'title': '#', 'field': 'protein_group', 'sorter': 'number', 'width': 44,
     'headerTooltip': 'Proteins numbered by their best score.'},
    {'title': 'Accession', 'field': 'accession', 'widthGrow': 2, 'tooltip': True},
    {'title': 'Description', 'field': 'description', **_HIDDEN},
    {'title': 'Forms', 'field': 'proteoforms', 'sorter': 'number', 'minWidth': 92,
     'headerTooltip': 'Number of distinct proteoforms (ProForma strings) identified for the protein.'},
    {'title': 'PrSMs', 'field': 'prsms', 'sorter': 'number', 'minWidth': 92,
     'headerTooltip': 'Number of spectra whose best hit is a proteoform of the protein.'},
    {'title': 'Score', 'field': 'best_score', 'sorter': 'number', 'minWidth': 92, **_HIDDEN,
     'headerTooltip': 'Best score of its PrSMs.'},
    {'title': 'q-value', 'field': 'best_qvalue', 'sorter': 'number', 'minWidth': 96, **_QVALUE,
     'headerTooltip': 'Lowest proteoform-level q-value of its PrSMs. Only PrSMs at q ≤ 0.01 are counted.'},
]
PROTEOFORM_COLUMNS = [
    {'title': '#', 'field': 'proteoform_id', 'sorter': 'number', 'width': 44,
     'headerTooltip': 'Proteoforms numbered by protein, then by their best score.'},
    {'title': 'Proteoform', 'field': 'proteoform', 'widthGrow': 3, 'tooltip': True,
     'headerTooltip': 'The proteoform in ProForma notation, as FLASHTnT reports it.'},
    {'title': 'First residue', 'field': 'start', 'sorter': 'number', **_POSITION, **_HIDDEN,
     'headerTooltip': 'First residue of the proteoform in the protein (1-based); "-" when the N-terminus was not determined.'},
    {'title': 'Last residue', 'field': 'end', 'sorter': 'number', **_POSITION, **_HIDDEN,
     'headerTooltip': 'Last residue of the proteoform in the protein (1-based); "-" when the C-terminus was not determined.'},
    {'title': 'Mass', 'field': 'mass', 'sorter': 'number', **_MASS, **_HIDDEN,
     'headerTooltip': 'Proteoform mass of its best PrSM.'},
    {'title': 'Mods', 'field': 'modifications', 'sorter': 'number', **_HIDDEN},
    {'title': 'PrSMs', 'field': 'prsms', 'sorter': 'number', 'minWidth': 92},
    {'title': 'Score', 'field': 'best_score', 'sorter': 'number', **_HIDDEN,
     'headerTooltip': 'Best score of its PrSMs.'},
    {'title': 'q-value', 'field': 'best_qvalue', 'sorter': 'number', 'minWidth': 96, **_QVALUE,
     'headerTooltip': 'Lowest proteoform-level q-value of its PrSMs. Only PrSMs at q ≤ 0.01 are counted.'},
]
PROTEOFORM_PRSM_COLUMNS = [
    {'title': 'Scan', 'field': 'Scan', 'sorter': 'number', 'minWidth': 84},
    {'title': 'RT', 'field': 'RetentionTime', 'sorter': 'number', **_fixed(), **_HIDDEN,
     'headerTooltip': 'Retention time.'},
    {'title': 'Precursor', 'field': 'PrecursorMass', 'sorter': 'number', 'minWidth': 110, **_HIDDEN,
     'headerTooltip': 'Precursor mass.', 'formatter': 'fixed',
     'formatterParams': {'precision': 4, 'minLength': 4, 'placeholders': {'-1': '-', '-1.0': '-'}}},
    {'title': 'Proteoform mass', 'field': 'ProteoformMass', 'sorter': 'number', **_MASS, **_HIDDEN},
    {'title': 'Fragments', 'field': 'MatchingFragments', 'sorter': 'number', 'minWidth': 112,
     'headerTooltip': 'Number of matched fragments.'},
    {'title': 'Score', 'field': 'Score', 'sorter': 'number', 'minWidth': 92},
    {'title': 'q-value', 'field': 'qvalue', 'sorter': 'number', 'minWidth': 96, **_QVALUE,
     'headerTooltip': 'Proteoform-level q-value. Only PrSMs at q ≤ 0.01 are counted.'},
]


def _hierarchy_columns(columns, with_qvalues):
    """Without a decoy search there are no q-values, so their column is hidden."""
    if with_qvalues:
        return columns
    return [
        {**column, 'visible': False} if column['field'] in ('best_qvalue', 'qvalue') else column
        for column in columns
    ]


# Rows of 33 px below the title bar and header; one table never exceeds TABLE_HEIGHT
_HIERARCHY_CHROME = 100
_HIERARCHY_ROW = 33


_IDENTIFICATION_ON_SELECT = {'tag': None, 'aa': None, 'mass': None}

# Components with nothing to show for a run deconvolved only (the sequence
# view then shows the sequence set in the view)
IDENTIFICATION_ONLY = (
    {c['name'] for c in TNT_COMPONENTS} - {c['name'] for c in DECONV_COMPONENTS} - {'sequence_view'}
)

ISOBARIC_ONLY = {c['name'] for c in ISOBARIC_COMPONENTS}


def _plex(fm, ds):
    """The labels a run was quantified with: (plex, channel names)."""
    info = fm.get_results(ds, ['isobaric_plex'])['isobaric_plex'] if fm.result_exists(ds, 'isobaric_plex') else {}
    plex = info.get('plex', isobaric.DEFAULT_PLEX)
    return plex, info.get('channels') or isobaric.channel_names(plex)


def _reporter_columns(identified):
    """The reporter table's columns: the spectrum and how many labels it
    carries. The intensity of each label is in the reporter ion plot of the
    selected spectrum, not one column per label."""
    columns = [
        {'title': 'Scan', 'field': 'Scan', 'sorter': 'number'},
        {'title': 'RT (min)', 'field': 'RT', 'sorter': 'number', **_fixed(2)},
        {'title': 'Precursor mass', 'field': 'PrecursorMass', 'sorter': 'number', **_fixed(2),
         'headerTooltip': 'Monoisotopic mass of the precursor (FLASHDeconv), in Daltons.'},
        {'title': 'Labels found', 'field': 'QuantifiedChannels', 'sorter': 'number',
         'headerTooltip': 'How many labels have a reporter ion in the spectrum.'},
    ]
    if identified:
        columns.append({'title': 'Identified as', 'field': 'Identified', 'sorter': 'string',
                        'headerTooltip': 'The protein of the best PrSM of the spectrum (FLASHTnT).'})
    return columns


# FLASHViewer's tables: as tall as their rows, up to 320 px below a 56 px title bar
# PrSM mass table: every deconvolved mass of the identified spectrum
# and the fragment ion it explains
PRSM_MASS_COLUMNS = [
    {'title': 'Mass ID', 'field': 'mass_id', 'sorter': 'number',
     'headerTooltip': 'The position of the mass in the deconvolved spectrum (1-based).'},
    {'title': 'Mono mass', 'field': 'MonoMass', 'sorter': 'number', **_fixed(),
     'headerTooltip': 'The monoisotopic mass of the deconvolved mass in Daltons.'},
    {'title': 'Intensity', 'field': 'SumIntensity', 'sorter': 'number', **_fixed(), 'responsive': 6,
     'headerTooltip': 'The total intensity of the mass across all isotopic peaks and charges.'},
    {'title': 'Charge', 'field': 'charges', 'responsive': 5,
     'headerTooltip': 'The charge states the mass was detected in.'},
    {'title': 'QScore', 'field': 'QScore', 'sorter': 'number', **_fixed(), 'responsive': 8,
     'headerTooltip': 'The quality score of the mass detection (higher is better).'},
    {'title': 'Status', 'field': 'status', 'responsive': 0,
     'headerTooltip': 'Whether the mass matches a fragment ion of the proteoform. The table filter button shows only matched or only not matched masses.'},
    {'title': 'Ion', 'field': 'ions', 'responsive': 0,
     'headerTooltip': 'The fragment ion(s) the mass matches; the numbers to the right are those of the first one.'},
    {'title': 'Position', 'field': 'position', 'sorter': 'number', **_fixed(0), 'responsive': 3,
     'headerTooltip': 'The residue next to the cleaved bond (1-based, in the database sequence).'},
    {'title': 'Theoretical mass', 'field': 'theoretical_mass', 'sorter': 'number', **_fixed(), 'responsive': 4,
     'headerTooltip': 'The theoretical mass of the matched fragment ion in Daltons.'},
    {'title': 'Mass error (Da)', 'field': 'mass_error', 'sorter': 'number', **_fixed(), 'responsive': 2,
     'headerTooltip': 'Observed minus theoretical mass in Daltons.'},
    {'title': 'Error (ppm)', 'field': 'ppm_error', 'sorter': 'number', **_fixed(2), 'responsive': 0,
     'headerTooltip': 'Observed minus theoretical mass in parts per million.'},
]

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


def _scan_methods(file_manager, dataset_id):
    """Fragmentation method of each spectrum of a dataset, by spectrum index,
    read from the deconvolved mzML (FLASHDeconv keeps the uploaded scans'
    activation methods). Empty when the file is not stored."""
    path = insight_data.result_path(file_manager, dataset_id, 'out_deconv_mzML')
    return () if path is None else mzml_scan_methods(path)


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


def _proteoform_heatmap(comp_name, file_manager, dataset_id):
    """Every PrSM as a point at its retention time and proteoform mass,
    coloured by score. A click selects the PrSM in the protein, proteoform and
    PrSM tables (and the cells that follow them); the selected PrSM is ringed,
    whichever cell selected it."""
    points = hierarchy_data.proteoform_map_from(hierarchy_data.prsms(file_manager, dataset_id))
    if not points.height:
        raise MissingResult('proteoform masses')
    return _component(
        Heatmap, comp_name, file_manager, dataset_id, ['protein_dfs', 'scan_table'],
        lambda: dict(
            data=points.lazy(),
            x_column='rt', y_column='mass', intensity_column='score',
            interactivity={
                'protein': 'index', 'pscan': 'scan_idx', 'scan': 'scan_idx',
                'proteoform': 'proteoform_id', 'protein_group': 'protein_group',
            },
            on_select=_IDENTIFICATION_ON_SELECT,
            highlight={'protein': 'index'},
            # A box or lasso selects the proteins of the PrSMs in it (GO terms)
            multi_select={go_terms.RUN_SET: 'accession'},
            hover_column='label',
            zoom_identifier=f'{comp_name}_zoom',
            title=f'Proteoform Heatmap ({points.height} PrSMs)',
            x_label='Retention Time', y_label='Proteoform Mass',
            colorscale='Portland', log_scale=False, intensity_label='Score',
            min_points=insight_data.HEATMAP_MIN_POINTS, downsample_zoomed=False,
            svg_filename='FLASHViewer-proteoform-heatmap', scroll_zoom=False,
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
            # FLASHViewer's FDR plot: Arial in the theme's text colour, densities unrounded
            font={'family': 'Arial'}, float_precision='float64',
        ),
        required=[target],
    )


# Ion type colours of the fragment error plot: prefix ions blue, suffix ions red
FRAGMENT_ION_COLORS = {
    'a': '#9ecae1', 'b': '#4292c6', 'c': '#08519c',
    'x': '#fcbba1', 'y': '#ef3b2c', 'z': '#a50f15', 'z+1': '#cb181d', 'z+2': '#67000d',
    'internal': '#8E44AD',
}


def _fragment_error_plot(comp_name, file_manager, dataset_id):
    """Mass error of the matched fragments of the selected proteoform; a click
    selects the residue next to the fragment's bond."""
    return _component(
        SeriesPlot, comp_name, file_manager, dataset_id,
        ['sequence_data', 'protein_dfs', 'scan_table', 'sequence_view', 'settings'],
        lambda: dict(
            data=validation_data.fragment_errors(file_manager, dataset_id),
            filters={'protein': 'proteoform_index'}, interactivity={'aa': 'residue'},
            x_column='theoretical_mass', y_column='ppm_error', series_column='ion_type',
            series_colors=FRAGMENT_ION_COLORS, mode='markers',
            y_range_mode='normal', hover_column='fragment',
            title='Fragment Mass Errors', x_label='Fragment mass (Da)', y_label='Error (ppm)',
            svg_filename='FLASHViewer-fragment-errors',
        ),
        required=['sequence_data', 'protein_dfs', 'scan_table', 'sequence_view'],
    )


def _isotope_envelope(comp_name, file_manager, dataset_id, tool):
    """Observed isotope envelope of the selected mass (the precursor while no
    mass is selected) over the averagine envelope."""
    scan_identifier = 'pscan' if tool == 'flashtnt' else 'scan'

    def data():
        envelopes = validation_data.isotope_envelopes(file_manager, dataset_id)
        observed = pl.col('side') == 'observed'
        return envelopes.with_columns(
            pl.when(observed).then(pl.col('scan_idx')).alias('top_scan'),
            pl.when(observed).then(pl.col('mass_idx')).alias('top_mass'),
            pl.when(~observed).then(pl.col('scan_idx')).alias('bottom_scan'),
            pl.when(~observed).then(pl.col('mass_idx')).alias('bottom_mass'),
        )

    return _component(
        MirrorPlot, f'{comp_name}-{tool}', file_manager, dataset_id,
        ['mass_table', 'threedim_SN_plot', 'scan_table'],
        lambda: dict(
            data=data(),
            filters_top={scan_identifier: 'top_scan', 'mass': 'top_mass'},
            filter_defaults_top={'mass': -1},
            filters_bottom={scan_identifier: 'bottom_scan', 'mass': 'bottom_mass'},
            filter_defaults_bottom={'mass': -1},
            x_column='mass', y_column='abundance', title='Isotope Envelope Fit',
            title_top='Observed (all charges)', title_bottom='Averagine',
            x_label='Mass (Da)', y_label='Relative abundance',
            styling={'unhighlightedColor': '#3366CC'},
        ),
    )


def _precursor_xic(comp_name, file_manager, dataset_id, tool):
    """The precursor mass of the selected fragmentation spectrum traced over
    the MS1 spectra, per charge."""
    scan_identifier = 'pscan' if tool == 'flashtnt' else 'scan'
    return _component(
        SeriesPlot, f'{comp_name}-{tool}', file_manager, dataset_id,
        ['mass_table', 'threedim_SN_plot', 'scan_table'],
        lambda: dict(
            data=validation_data.precursor_xics(file_manager, dataset_id),
            filters={scan_identifier: 'scan_idx'},
            x_column='rt', y_column='intensity', series_column='series',
            series_names={'All charges': 'All charges'}, mode='lines+markers',
            title='Precursor XIC by Charge', x_label='Retention time (s)', y_label='Intensity',
            svg_filename='FLASHViewer-precursor-xic',
        ),
    )


def _chromatogram(comp_name, file_manager, dataset_id):
    """TIC, base peak and deconvolved intensity over retention time; a click
    selects the spectrum."""
    level = qc_data.chromatogram_level(file_manager, dataset_id)
    return _component(
        SeriesPlot, comp_name, file_manager, dataset_id,
        ['scan_table', 'anno_spectrum', 'deconv_spectrum'],
        lambda: dict(
            data=qc_data.chromatograms(file_manager, dataset_id, level),
            interactivity={'scan': 'index'},
            x_column='x', y_column='y', series_column='series',
            series_names={name: name for name in qc_data.CHROMATOGRAM_SERIES.values()},
            series_colors=dict(zip(
                qc_data.CHROMATOGRAM_SERIES.values(), ['#1f77b4', '#ff7f0e', '#2ca02c']
            )),
            mode='lines', y_range_mode='nonnegative',
            title=f'Chromatogram (MS{level})', x_label='Retention time (s)',
            y_label='Intensity', svg_filename='FLASHViewer-chromatogram',
        ),
        required=['scan_table', 'anno_spectrum'],
    )


def _proteins_with_characterization(file_manager, dataset_id):
    """Protein table rows with the characterization columns (empty when the
    results have no sequence data)."""
    proteins = insight_data.proteins(file_manager, dataset_id)
    if any(
        insight_data.result_path(file_manager, dataset_id, name) is None
        for name in ('sequence_data', 'sequence_view')
    ):
        return proteins.with_columns(
            pl.lit(None, dtype=pl.Float64).alias('Cleavages'),
            pl.lit(None, dtype=pl.String).alias('LocalizedMods'),
            pl.lit(None, dtype=pl.String).alias('Characterization'),
        )
    return proteins.join(
        validation_data.characterization(file_manager, dataset_id),
        left_on='index', right_on='proteoform_index', how='left', maintain_order='left',
    )


def _stored_sequence(file_manager):
    """The sequence saved for the workspace (dict), or an empty dict."""
    if file_manager.result_exists('sequence', 'sequence'):
        return file_manager.get_results('sequence', 'sequence')['sequence']
    return {}


def _deconv_sequence(file_manager, state_manager):
    """Sequence data for the FLASHDeconv sequence view: the sequence set in the
    view's "Change sequence" dialog, else the one saved earlier, else an empty
    sequence (the view then offers its dialog to set one).

    A sequence set in the view is saved for the workspace, so it stays after a
    reload and unlocks the components that need a sequence."""
    stored = _stored_sequence(file_manager)
    sequence = state_manager.get_selection('sequence') or stored.get('input_sequence') or ''
    mod_ids = list(_fixed_modification_ids(stored, state_manager))
    if sequence and (sequence != stored.get('input_sequence') or mod_ids != stored.get('fixed_modifications')):
        file_manager.store_data('sequence', 'sequence', {
            'input_sequence': sequence,
            'fixed_mod_cysteine': None,
            'fixed_mod_methionine': None,
            'fixed_modifications': mod_ids,
        })
    return _sequence_data(sequence, tuple(mod_ids))


def _fixed_modification_ids(stored, state_manager):
    """The fixed modifications of the FLASHDeconv sequence view: those picked
    in the view's "Change sequence" dialog, else those saved with the sequence."""
    return modifications.fixed_ids(stored, state_manager.get_selection(FIXED_MODS_IDENTIFIER))


@lru_cache(maxsize=1)
def _modification_picker_library():
    """The modification library as the sequence view's dialog takes it."""
    return [{'id': e['id'], 'mass': e['mass'], 'site': e['site']} for e in modifications.library()]


@st.cache_data(max_entries=4, show_spinner=False)
def _sequence_data(sequence, mod_ids):
    fixed_mods = {'fixed_modifications': list(mod_ids)}
    data = getFragmentDataFromSeq(sequence, fixed_mods=fixed_mods)
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


def _hierarchy_table(comp_name, file_manager, dataset_id):
    """One level of the protein → proteoform → PrSM tables and its height.
    Selecting a proteoform or PrSM selects its (best) PrSM in the other
    cells, as a protein table row does. The tables share one height: that of
    the level with the most rows it can show at once (all proteins, the
    proteoforms of the protein with the most, the PrSMs of the proteoform
    with the most), up to TABLE_HEIGHT."""
    hits = hierarchy_data.prsms(file_manager, dataset_id)
    with_qvalues = hierarchy_data.has_qvalues(hits)
    common = dict(
        reclick=RECLICK, scroll_position=SCROLL, reset_sort=True,
        pagination=False, float_precision='float64',
        # Columns share the cell's width, so a narrow cell never scrolls sideways
        layout='fitColumns',
        # Stay at the shared height below, whatever the rows
        fill_height=True,
    )
    if comp_name == 'protein_group_table':
        build = lambda: dict(
            data=hierarchy_data.protein_groups_from(hits).lazy(),
            column_definitions=_hierarchy_columns(PROTEIN_GROUP_COLUMNS, with_qvalues),
            title=hierarchy_data.level_title('Proteins', hierarchy_data.key_figures(hits)['proteins']),
            index_field='protein_group', go_to_fields=['accession'],
            interactivity={'protein_group': 'protein_group'},
            multi_select={go_terms.RUN_SET: 'accession'}, **common,
        )
    elif comp_name == 'proteoform_table':
        build = lambda: dict(
            data=hierarchy_data.proteoforms_from(hits).lazy(),
            column_definitions=_hierarchy_columns(PROTEOFORM_COLUMNS, with_qvalues),
            title=hierarchy_data.level_title('Proteoforms', hierarchy_data.key_figures(hits)['proteoforms']),
            index_field='proteoform_id', go_to_fields=['proteoform'],
            filters={'protein_group': 'protein_group'},
            interactivity={'proteoform': 'proteoform_id', 'protein': 'best_index', 'pscan': 'scan_idx', 'scan': 'scan_idx'},
            multi_select={go_terms.RUN_SET: 'accession'},
            on_select=_IDENTIFICATION_ON_SELECT, **common,
        )
    else:
        build = lambda: dict(
            data=hierarchy_data.prsm_rows_from(hits).lazy(),
            column_definitions=_hierarchy_columns(PROTEOFORM_PRSM_COLUMNS, with_qvalues),
            title=hierarchy_data.level_title('PrSMs', hierarchy_data.key_figures(hits)['prsms']),
            index_field='index', go_to_fields=['Scan'],
            filters={'proteoform': 'proteoform_id'},
            interactivity={'protein': 'index', 'pscan': 'scan_idx', 'scan': 'scan_idx'},
            multi_select={go_terms.RUN_SET: 'accession'},
            on_select=_IDENTIFICATION_ON_SELECT, **common,
        )
    # All three are as tall as the one with the most rows, so that tables
    # sitting side by side in a layout row line up
    rows = max(hierarchy_data.most_rows(hits).values())
    height = min(TABLE_HEIGHT, _HIERARCHY_CHROME + _HIERARCHY_ROW * max(rows, 1))
    table = _component(
        Table, comp_name, file_manager, dataset_id, ['protein_dfs', 'scan_table'], build,
    )
    return table, height


def _build(comp_name, file_manager, dataset_id, tool, state_manager):
    """Insight component of a layout cell, and keyword arguments for rendering it."""
    fm, ds = file_manager, dataset_id

    if comp_name in HEATMAPS:
        return _heatmap(comp_name, fm, ds), {'height': 450}

    if comp_name == 'proteoform_heatmap':
        return _proteoform_heatmap(comp_name, fm, ds), {'height': 450}

    if comp_name == 'scan_table':
        return _component(
            Table, comp_name, fm, ds, ['scan_table'],
            lambda: dict(
                data=insight_data.scan_table(fm, ds),
                column_definitions=SCAN_COLUMNS, title='Scan Table',
                index_field='index', go_to_fields=['index', 'Scan'],
                interactivity={'scan': 'index'}, on_select={'mass': 0},
                reclick=RECLICK, scroll_position=SCROLL, reset_sort=True, pagination=False, float_precision='float64',
            ),
        ), {'height': TABLE_HEIGHT}

    if comp_name == 'reporter_table':
        plex, _ = _plex(fm, ds)
        return _component(
            Table, comp_name, fm, ds, ['isobaric_dfs', 'scan_table'],
            lambda: dict(
                data=insight_data.reporter_table(fm, ds),
                column_definitions=_reporter_columns(fm.result_exists(ds, 'protein_dfs')),
                title=f'Reporter ions ({isobaric.label(plex)})',
                index_field='scan_idx', go_to_fields=['Scan'], interactivity={'scan': 'scan_idx'},
                reclick=RECLICK, scroll_position=SCROLL, reset_sort=True, pagination=False, float_precision='float64',
            ),
            required=['isobaric_dfs', 'scan_table'],
        ), {'height': TABLE_HEIGHT}

    if comp_name == 'reporter_ions':
        plex, channels = _plex(fm, ds)
        mzs = dict(isobaric.PLEXES.get(plex, ('', []))[1])
        return _component(
            LinePlot, comp_name, fm, ds, ['isobaric_dfs', 'scan_table'],
            lambda: dict(
                data=insight_data.reporter_ions(fm, ds, channels, [mzs.get(c, 0.0) for c in channels]),
                filters={'scan': 'scan_idx'}, x_column='mz', y_column='intensity',
                annotation_column='label', title='Reporter ions of the selected spectrum',
                x_label='m/z', y_label='Intensity',
            ),
        ), {}

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
                reclick=RECLICK, scroll_position=SCROLL, reset_sort=True, pagination=False, float_precision='float64',
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

    if comp_name == 'fragment_error_plot':
        return _fragment_error_plot(comp_name, fm, ds), {'height': 400}

    if comp_name == 'isotope_envelope':
        return _isotope_envelope(comp_name, fm, ds, tool), {'height': 400}

    if comp_name == 'precursor_xic':
        return _precursor_xic(comp_name, fm, ds, tool), {'height': 400}

    if comp_name == 'chromatogram':
        return _chromatogram(comp_name, fm, ds), {'height': 350}

    if comp_name in ('fdr_plot', 'id_fdr_plot'):
        return _fdr_plot(comp_name, fm, ds), {'height': 400}

    if comp_name in ('protein_group_table', 'proteoform_table', 'proteoform_prsm_table'):
        table, height = _hierarchy_table(comp_name, fm, ds)
        return table, {'height': height}

    if comp_name == 'protein_table':
        return _component(
            Table, comp_name, fm, ds,
            ['protein_dfs', 'scan_table', 'sequence_data', 'sequence_view', 'settings', 'spec2_tsv'],
            lambda: dict(
                data=_proteins_with_characterization(fm, ds).join(
                    prsm_overview.protein_overview(fm, ds), on='index', how='left', maintain_order='left',
                ),
                column_definitions=PROTEIN_COLUMNS, title='Protein Table',
                # A button per row opens the whole row, hidden columns included
                row_overview="Scan {Scan}: {accession}",
                index_field='index', go_to_fields=['Scan', 'accession'],
                initial_sort=[{'column': 'Score', 'dir': 'desc'}],
                toggle_filters=[{
                    'label': 'Best per spectrum', 'column': 'is_best', 'group_column': 'Scan',
                    # FLASHViewer's checkbox colour
                    'color': '#6200EE',
                }],
                interactivity={'protein': 'index', 'pscan': 'scan_idx', 'scan': 'scan_idx'},
                multi_select={go_terms.RUN_SET: 'accession'},
                on_select={'tag': None, 'aa': None, 'mass': None},
                reclick=RECLICK, scroll_position=SCROLL, reset_sort=True, pagination=False, float_precision='float64',
            ),
            required=['protein_dfs', 'scan_table'],
        ), {'height': TABLE_HEIGHT}

    if comp_name == 'prsm_mass_table':
        return _component(
            Table, comp_name, fm, ds,
            ['sequence_data', 'protein_dfs', 'scan_table', 'sequence_view', 'mass_table', 'settings'],
            lambda: dict(
                data=prsm_data.prsm_masses(fm, ds),
                column_definitions=PRSM_MASS_COLUMNS, title='PrSM Mass Table',
                index_field='mass_idx', go_to_fields=['mass_id', 'ions'],
                filters={'protein': 'proteoform_index'},
                toggle_filters=[{
                    'label': 'Matched masses only', 'column': 'matched', 'default': False,
                    'color': '#6200EE',
                }],
                # A row selects its mass: the spectra highlight it and the
                # Sequence View shows its fragment
                interactivity={'mass': 'mass_idx'}, default_row_selects=False,
                reclick=RECLICK, scroll_position=SCROLL, reset_sort=True,
                pagination=False, float_precision='float64',
            ),
            required=['sequence_data', 'protein_dfs', 'scan_table', 'sequence_view'],
        ), {'height': TABLE_HEIGHT}

    if comp_name == 'tag_table':
        return _component(
            Table, comp_name, fm, ds, ['tag_dfs', 'scan_table'],
            lambda: dict(
                # Rows come sorted by score, with no sort shown on the header
                data=insight_data.tags(fm, ds).sort(
                    'Score', descending=True, nulls_last=True, maintain_order=True
                ),
                column_definitions=TAG_COLUMNS, title='Tag Table',
                index_field='TagIndex', go_to_fields=['Scan', 'StartPos', 'EndPos', 'TagSequence'],
                filters={'scan': 'scan_idx'}, range_filters={'aa': ['StartPos', 'EndPos']},
                interactivity={'tag': 'TagIndex'}, on_select={'mass': None},
                reclick=RECLICK, scroll_position=SCROLL, reset_sort=True,
                reset_sort_on=['protein'], pagination=False, float_precision='float64',
            ),
        ), {'height': TABLE_HEIGHT}

    if comp_name == 'sequence_view':
        # The proteoform of an identified run, else the sequence set in the view
        if tool == 'flashtnt' or (tool == 'topdown' and fm.result_exists(ds, 'protein_dfs')):
            component = _component(
                ProteinSequenceView, f'{comp_name}-{tool}', fm, ds,
                ['sequence_data', 'protein_dfs', 'scan_table', 'sequence_view', 'settings',
                 'out_deconv_mzML'],
                lambda: dict(
                    # Without FLASHAnnotator fragments, the proteoform's scan
                    # selects the ion series of its fragmentation method
                    data=insight_data.with_fragmentation(
                        insight_data.proteoform_sequences(fm, ds), _scan_methods(fm, ds), 'scan_idx'
                    ),
                    filters={'protein': 'proteoform_index'},
                    mass_identifier='mass', aa_identifier='aa', fragment_identifier='fragment',
                    fragmentation_column='fragmentation',
                    settings=fm.get_results(ds, ['settings'])['settings'],
                    entry_title='Proteoform', observed_mass_label='Observed proteoform mass',
                    **SEQUENCE_VIEW_OPTIONS,
                ),
                required=['sequence_data', 'protein_dfs', 'scan_table', 'sequence_view', 'settings'],
            )
            tag = _selected_tag(fm, ds, state_manager)
            # Tags FLASHTnT could not place on the proteoform have no positions
            if tag is None or tag['StartPos'] is None or tag['EndPos'] is None:
                component.set_tag(None)
            else:
                component.set_tag(int(tag['StartPos']), int(tag['EndPos']))
            return component, {}

        sequence = _deconv_sequence(fm, state_manager)
        component = _component(
            ProteinSequenceView, f'{comp_name}-{tool}', fm, ds,
            ['sequence_view', 'out_deconv_mzML'],
            lambda: dict(
                # Each scan selects the ion series of its own fragmentation method
                data=insight_data.with_fragmentation(
                    insight_data.sequence_view_spectra(fm, ds), _scan_methods(fm, ds), 'index'
                ),
                filters={'scan': 'index'}, sequence=sequence,
                mass_identifier='mass', aa_identifier='aa', sequence_identifier='sequence',
                fragmentation_column='fragmentation',
                **SEQUENCE_VIEW_OPTIONS,
            ),
            extra=sequence, required=['sequence_view'],
        )
        # The "Change sequence" dialog also picks the fixed modifications
        if hasattr(component, 'set_modification_picker'):
            mod_ids = _fixed_modification_ids(_stored_sequence(fm), state_manager)
            component.set_modification_picker(
                _modification_picker_library(), FIXED_MODS_IDENTIFIER, mod_ids
            )
        return component, {}

    if comp_name == 'structure_view':
        component = _component(
            ProteinStructureView, f'{comp_name}-{tool}', fm, ds,
            ['sequence_data', 'protein_dfs', 'scan_table', 'sequence_view', 'settings'],
            lambda: dict(
                data=insight_data.proteoform_sequences(fm, ds),
                filters={'protein': 'proteoform_index'}, mass_column='MonoMass',
                aa_identifier='aa', fragment_identifier='fragment',
                settings=fm.get_results(ds, ['settings'])['settings'],
                empty_message='Select a proteoform in the protein table to show its structure.',
            ),
        )
        component.set_structure(_structure(fm, ds, state_manager.get_selection('protein')))
        component.on_structure_action(structure_action_handler(structure_dir(fm)))
        return component, {}

    raise ValueError(f'Unknown viewer component {comp_name!r}')


def _structure(file_manager, dataset_id, protein_index):
    """The selected proteoform's structure for the Structure View (see
    src/render/structure.py), or None while no proteoform is selected."""
    if protein_index is None:
        return None
    proteins = file_manager.get_results(dataset_id, ['protein_dfs'])['protein_dfs']
    match = proteins.loc[proteins['index'] == protein_index, 'accession']
    if len(match) == 0:
        return None
    structure = get_structure(
        StructureStore(structure_dir(file_manager)), match.iloc[0],
        download=download_enabled(st.session_state.get('settings')),
    )
    error = upload_errors().get(structure_key(match.iloc[0]))
    if error is not None:
        structure['upload_error'] = error
    if structure['status'] == 'ok':
        structure['label'] = (
            f"AlphaFold DB model {structure.get('name') or structure['accession']}"
            if structure.get('source') == 'alphafold'
            else f"Uploaded structure {structure.get('name') or ''}".strip()
        )
    return structure


RUN_QC_PLOTS = [
    ('qc_qscore', qc_data.qscore_distribution, ['mass_table', 'scan_table'],
     'QScore Distribution', 'QScore', 'Masses', 'lines', 'nonnegative'),
    ('qc_charge', qc_data.charge_distribution, ['mass_table', 'scan_table'],
     'Charge States', 'Charge', 'Masses', 'lines+markers', 'nonnegative'),
    ('qc_mass', qc_data.mass_distribution, ['mass_table', 'scan_table'],
     'Mass Distribution', 'Monoisotopic mass (Da)', 'Masses', 'lines', 'nonnegative'),
    ('qc_masses_per_spectrum', qc_data.masses_per_spectrum, ['scan_table'],
     'Masses per Spectrum', 'Retention time (s)', 'Masses', 'markers', 'nonnegative'),
]
ID_QC_PLOTS = [
    ('qc_identifications', qc_data.identifications_over_time, ['protein_dfs', 'scan_table'],
     'Identifications over Time', 'Retention time (s)', 'PrSMs per minute', 'lines+markers',
     'nonnegative'),
    ('qc_fragment_error', qc_data.fragment_error_over_time,
     ['sequence_data', 'protein_dfs', 'scan_table', 'sequence_view', 'settings'],
     'Fragment Mass Error over Time', 'Retention time (s)', 'Median error (ppm)', 'markers',
     'normal'),
]


def _render_qc(plots, file_manager, dataset_id, tool, state_manager, key):
    """A panel of run-level plots, two per row; plots whose results are
    missing are left out."""
    components = []
    for name, builder, sources, title, x_label, y_label, mode, y_range in plots:
        try:
            components.append(_component(
                SeriesPlot, f'{name}-{tool}', file_manager, dataset_id, sources,
                lambda builder=builder, name=name, title=title, x_label=x_label,
                y_label=y_label, mode=mode, y_range=y_range: dict(
                    data=builder(file_manager, dataset_id),
                    x_column='x', y_column='y', series_column='series',
                    mode=mode, title=title, x_label=x_label, y_label=y_label,
                    y_range_mode=y_range, svg_filename=f'FLASHViewer-{name}',
                ),
                required=[source for source in sources if source != 'settings'],
            ))
        except MissingResult:
            continue
    if not components:
        st.info('This dataset has no data for the run QC plots.')
        return
    for row_start in range(0, len(components), 2):
        columns = st.columns(2)
        for offset, component in enumerate(components[row_start:row_start + 2]):
            with columns[offset]:
                component(
                    key=f'{key}_{row_start + offset}', state_manager=state_manager, height=320,
                )


def _render_cell(comp_name, file_manager, dataset_id, tool, state_manager, key):
    if (
        tool == 'topdown' and comp_name in IDENTIFICATION_ONLY
        and not file_manager.result_exists(dataset_id, 'protein_dfs')
    ):
        st.info('This run was deconvolved only. Run it with **Identify proteoforms** on to fill this view.')
        return
    if comp_name in ISOBARIC_ONLY and not file_manager.result_exists(dataset_id, 'isobaric_dfs'):
        st.info('This run has no reporter ions. Run it with **Quantification: Isobaric labels** to fill this view.')
        return
    if comp_name == 'internal_fragment_map':
        st.info('The internal fragment map is not available yet.')
        return
    if comp_name == 'run_qc':
        _render_qc(RUN_QC_PLOTS, file_manager, dataset_id, tool, state_manager, key)
        return
    if comp_name == 'id_qc':
        _render_qc(ID_QC_PLOTS, file_manager, dataset_id, tool, state_manager, key)
        return
    if comp_name == 'go_enrichment':
        from src.render import go_view
        go_view.render_run(file_manager, dataset_id, state_manager, key)
        return
    try:
        component, render_args = _build(comp_name, file_manager, dataset_id, tool, state_manager)
    except MissingResult as e:
        if str(e) == 'proteoform masses':
            st.info(
                'No PrSM of this dataset has a proteoform or precursor mass, '
                'so the proteoform heatmap has nothing to show.'
            )
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
