"""Viewer grid built from openms-insight components.

Every layout cell is one Insight component. The components of one grid share a
StateManager, so a click in one (a scan, a mass, a proteoform, a tag, a residue)
filters or highlights the others. Each component keeps a preprocessed cache next
to the workflow results (cache/insight/<dataset>/), built on first view and
rebuilt when the results it was built from change.
"""

import hashlib
import inspect
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

# Where the proteoforms of a protein lie on its sequence; an openms-insight
# release without it shows a note in its place
SequenceCoverage = getattr(openms_insight, 'SequenceCoverage', None)

# The heatmaps of MS data also draw their points in 3D (one cell instead of a
# heatmap and a 3D map); an openms-insight release without it draws 2D only
PeakMap = getattr(openms_insight, 'PeakMap', None)

try:  # the comparison plots
    from openms_insight import PCAPlot, VolcanoPlot
    _QUANT_ZOOMABLE = (PCAPlot, VolcanoPlot)
except ImportError:
    _QUANT_ZOOMABLE = ()

from src.render import (
    envelope_overlay, go_terms, hierarchy_data, insight_data, linked_spectra, mass_shifts, modifications,
    prsm_data, prsm_overview, qc_data, validation_data,
)
from src.render.cache_dirs import build_lock, insight_cache_dir
from src.render.fragmentation import mzml_scan_methods
from src.quant import isobaric
from src.common.plot_style import PEAK, PICKED, spectrum_styling
from src.render import proteoform_edits, teaching
from src.render.sequence import getFragmentDataFromSeq
from src.render.structure import (
    StructureStore, download_enabled, get_structure, structure_dir, structure_key,
)
from src.render.structure_files import structure_action_handler, upload_errors
from src.viz import columns as columns_lib, plumbing
from src.viz.ms_data import chromatogram, peak_map_3d, precursor_spectrum, raw_heatmap, scan_table, spectrum
from src.viz.features import (
    deconvolved_heatmap, isotope_envelope, mass_table, precursor_xic, qscore_fdr, run_qc, signal_noise_3d,
)
from src.viz.identification import (
    _hierarchy, fragment_errors, identification_qc, mass_shift_histogram, mass_shift_table, prsm_mass_table,
    prsm_table, proteoform_heatmap, proteoform_prsms, proteoforms, protein_groups, score_distribution,
    sequence_coverage, sequence_view, structure_view, tag_table,
)
from src.viz.quantification import reporter_ions, reporter_table


def _code_version():
    """Changes whenever the configuration of a component can: with this module,
    the frame builders or openms-insight. Part of every cache id, so caches
    are rebuilt after an update."""
    digest = hashlib.md5(openms_insight.__version__.encode())
    for module in (
        __file__, insight_data.__file__, validation_data.__file__, qc_data.__file__,
        hierarchy_data.__file__, prsm_data.__file__, go_terms.__file__, mass_shifts.__file__,
        str(Path(__file__).with_name('quant_cells.py')),
        # The library panels the cells are configured with
        *sorted(str(f) for f in Path(__file__).parents[1].joinpath('viz').rglob('*.py')),
    ):
        try:
            digest.update(Path(module).read_bytes())
        except OSError:
            # Frozen builds may not ship the sources: the version alone
            digest.update(module.encode())
    return digest.hexdigest()


CODE_VERSION = _code_version()


_fixed = columns_lib.fixed
_DASH_FOR_MISSING = columns_lib.DASH_FOR_MISSING
_BLANK_FOR_MISSING = columns_lib.BLANK_FOR_MISSING
_MASS = columns_lib.MASS
_HIDDEN = columns_lib.HIDDEN

# The columns of every table come from the visualization library (src/viz)
SCAN_COLUMNS = scan_table.COLUMNS
MASS_COLUMNS = mass_table.COLUMNS
PROTEIN_COLUMNS = prsm_table.COLUMNS
TAG_COLUMNS = tag_table.COLUMNS
SPECTRUM_COLUMNS = spectrum.COLUMNS
SPECTRUM_TITLES = spectrum.TITLES
# Selections of the precursor spectrum cell: the precursor spectrum, the
# precursor mass in it, and the spectrum they were looked up for
SURVEY, SURVEY_MASS, SURVEY_FOR = 'survey', 'survey_mass', 'survey_for'
SEQUENCE_VIEW_OPTIONS = dict(
    sequence_view.OPTIONS, annotation_title='FLASHAnnotator', svg_filename='FLASHViewer-sequence',
)


# The sequence view's annotation names the tool and its download FLASHViewer
_SEQUENCE_VIEW_NAMES = dict(annotation_title='FLASHAnnotator', svg_filename='FLASHViewer-sequence')


def _svg(name):
    """FLASHApp names its downloaded plots after FLASHViewer."""
    return {'svg_filename': f'FLASHViewer-{name}'}

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
# Selection the identification sequence view writes its edits to (see
# src.render.proteoform_edits)
EDIT_IDENTIFIER = 'sequence_edit'
# The fragment ions the sequence view matched, marked on the spectra; left out
# with an openms-insight release without it
FRAGMENT_MATCHES_IDENTIFIER = 'fragment_matches'
FRAGMENT_MARKS = hasattr(DeconvolvedSpectrum, 'set_fragment_matches_identifier')
# openms-insight releases whose LinePlot marks labelled x positions (x_markers)
X_MARKERS = 'x_markers' in inspect.signature(LinePlot.__init__).parameters
# Reporter ions of neighbouring TMT channels are 0.0063 m/z apart
REPORTER_TOLERANCE = reporter_ions.TOLERANCE


def reporter_plot_args(channels, mzs, markers=X_MARKERS):
    """How the reporter ion plot shows its channels (see
    ``src.viz.quantification.reporter_ions.plot_args``)."""
    return reporter_ions.plot_args(channels, mzs, markers)


# Plots keep their zoom in the grid's state (zoom@<component>), so a view can
# store it with its selection; left out with a release without zoom_identifier
ZOOMABLE = (DeconvolvedSpectrum, Heatmap, LinePlot, MirrorPlot, SeriesPlot)


def _zoom_args(cls, name):
    """``zoom_identifier`` for a plot that keeps its zoom, named after its
    layout component (the cache name without its ``-<tool>`` suffix)."""
    if not issubclass(cls, ZOOMABLE + _QUANT_ZOOMABLE):
        return {}
    if 'zoom_identifier' not in inspect.signature(cls.__init__).parameters:
        return {}
    return {'zoom_identifier': zoom_identifier(name.split('-', 1)[0])}


def _identification_settings(file_manager, dataset_id):
    """Fragment settings of the identification sequence views: the run's
    stored settings, and the residues a proteoform's truncated termini leave
    out marked by a break line per terminus, labelled with their positions
    (count and residues on hover), which keeps a short proteoform of a long protein on a few rows.
    An openms-insight release without the option ignores it."""
    settings = dict(file_manager.get_results(dataset_id, ['settings'])['settings'] or {})
    settings['show_truncations'] = False
    return settings


def _spectrum_key(column):
    """Keyword arguments naming the column of a sequence view's data that holds
    the scan index of its spectrum, so the spectra mark its fragment ions only
    on that spectrum."""
    return {'spectrum_key_column': column} if FRAGMENT_MARKS else {}

RECLICK = 'unhighlight'
# A selected row is scrolled to the top of its table, a clicked row only when cut
# off
SCROLL = 'top'
# Identification levels: proteins, the proteoforms of a protein,
# the PrSMs of a proteoform (src/render/hierarchy_data.py)
_QVALUE = columns_lib.QVALUE
_POSITION = columns_lib.POSITION
PROTEIN_GROUP_COLUMNS = protein_groups.COLUMNS
PROTEOFORM_COLUMNS = proteoforms.COLUMNS
PROTEOFORM_PRSM_COLUMNS = proteoform_prsms.COLUMNS
# Rows below the title bar and header of the hierarchy tables
HIERARCHY_HEIGHT = _hierarchy.HEIGHT

# Components with nothing to show for a run deconvolved only, or without labels
from src.render.layout_editor.catalog import (  # noqa: E402
    IDENTIFICATION, IDENTIFICATION_ONLY, ISOBARIC, ISOBARIC_ONLY, TOPDOWN_COMPONENTS,
)
from src.layout_manager import missing_note, zoom_identifier  # noqa: E402

_COMPONENT_LABELS = {c['name']: c['label'] for c in TOPDOWN_COMPONENTS}


def _plex(fm, ds):
    """The labels a run was quantified with: (plex, channel names)."""
    info = fm.get_results(ds, ['isobaric_plex'])['isobaric_plex'] if fm.result_exists(ds, 'isobaric_plex') else {}
    plex = info.get('plex', isobaric.DEFAULT_PLEX)
    return plex, info.get('channels') or isobaric.channel_names(plex)


PRSM_MASS_COLUMNS = prsm_mass_table.COLUMNS
TABLE_HEIGHT = columns_lib.TABLE_HEIGHT


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


def _component(cls, name, file_manager, dataset_id, sources, build, extra=None, required=None,
               render=None):
    """The Insight component `name` of a dataset: loaded from its cache, or
    built with the keyword arguments `build()` returns.

    `sources` are the results the cache is built from; `required` (default:
    all of them) must exist. `render` are render-time arguments, not part of
    the cache, passed whether the cache is loaded or built.
    """
    for name_tag in sources if required is None else required:
        if insight_data.result_path(file_manager, dataset_id, name_tag) is None:
            raise MissingResult(name_tag)
    cache_path = insight_cache_dir(file_manager.cache_path, dataset_id)
    cache_id = f'{name}-{_stamp(file_manager, dataset_id, sources, extra)}'
    zoom = {**_zoom_args(cls, name), **(render or {})}
    try:
        return cls(cache_id=cache_id, cache_path=str(cache_path), **zoom)
    except CacheMissError:
        pass
    with build_lock(cache_path, cache_id):
        # Another tab may have built it while this one waited
        try:
            return cls(cache_id=cache_id, cache_path=str(cache_path), **zoom)
        except CacheMissError:
            pass
        # Caches built from earlier results or configurations are not used again.
        # Not for caches that depend on the selection (`extra`, e.g. a typed
        # sequence): two grids on a page may each need their own.
        if extra is None and cache_path.exists():
            for old in cache_path.glob(f'{name}-*'):
                if old.name != cache_id:
                    shutil.rmtree(old, ignore_errors=True)
        return cls(cache_id=cache_id, cache_path=str(cache_path), **zoom, **build())


def _heatmap(comp_name, file_manager, dataset_id, view='2d'):
    """A heatmap of MS data, which its toolbar also draws in 3D; ``view``
    is the drawing it opens with."""
    result, title, y_label, sets_mass = HEATMAPS[comp_name]
    shared = dict(min_points=insight_data.HEATMAP_MIN_POINTS, **_svg('heatmap'))
    cls = PeakMap or Heatmap
    if insight_data.result_path(file_manager, dataset_id, result) is None:
        raise MissingResult(result)
    # The levels the parse downsampled; without them Insight builds its own
    levels = insight_data.heatmap_levels(file_manager, dataset_id, result)
    return _component(
        cls, comp_name, file_manager, dataset_id, [result, *(levels or [])],
        lambda: dict(
            data=insight_data.scan(file_manager, dataset_id, result),
            levels=None if levels is None else [
                insight_data.scan(file_manager, dataset_id, name) for name in levels
            ],
            **(deconvolved_heatmap.settings(title, **shared) if sets_mass
               else raw_heatmap.settings(title, y_label, **shared)),
        ),
        render=dict(view=view) if PeakMap is not None else None,
    )


# Whether Plot3D draws the points of a selected value apart; left out with an
# openms-insight release without it
_PLOT3D_PARAMETERS = inspect.signature(Plot3D.__init__).parameters
_PLOT3D_HIGHLIGHT = 'highlight' in _PLOT3D_PARAMETERS


def _ms1_3d_map(comp_name, file_manager, dataset_id, tool, state_manager):
    """The MS1 peaks of the run in 3D (m/z, retention time, intensity), the
    strongest peaks of each spectrum. The precursor spectrum of the selected
    spectrum (an MS1 spectrum itself) is drawn in orange."""
    if insight_data.result_path(file_manager, dataset_id, 'ms1_raw_heatmap') is None:
        raise MissingResult('ms1_raw_heatmap')
    if _PLOT3D_HIGHLIGHT and insight_data.result_path(file_manager, dataset_id, 'scan_table') is not None:
        try:
            _follow_selected_spectrum(
                file_manager, dataset_id, state_manager, 'pscan' if tool == 'flashtnt' else 'scan',
            )
        except MissingResult:
            pass
    return _component(
        Plot3D, comp_name, file_manager, dataset_id, ['ms1_raw_heatmap'],
        lambda: dict(
            data=insight_data.ms1_peak_map(file_manager, dataset_id),
            # Options of a newer openms-insight (the box's shape) are left out
            **plumbing.only_accepted(Plot3D, **peak_map_3d.settings(
                insight_data.MS1_MAP_PEAKS_PER_SPECTRUM, PEAK,
                SURVEY if _PLOT3D_HIGHLIGHT else None, PICKED, **_svg('ms1-3d-map'),
            )),
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
            **proteoform_heatmap.settings(
                points.height, protein_set=go_terms.RUN_SET, min_points=insight_data.HEATMAP_MIN_POINTS,
                **_svg('proteoform-heatmap'),
            ),
        ),
    )


def _fdr_plot(comp_name, file_manager, dataset_id):
    if comp_name == 'id_fdr_plot':
        return _id_score_plot(comp_name, file_manager, dataset_id)
    target, decoy = 'density_target', 'density_decoy'
    return _component(
        SeriesPlot, comp_name, file_manager, dataset_id, [target, decoy],
        lambda: dict(
            data=insight_data.fdr_densities(file_manager, dataset_id, target, decoy),
            **qscore_fdr.settings(),
        ),
        required=[target],
    )


def _id_score_plot(comp_name, file_manager, dataset_id):
    """Score distribution of the identified proteoforms, targets and (when
    FLASHTnT kept them) decoys. The q-values of the reported hits are almost
    all 0, so their distribution left the plot empty."""
    def build():
        data = insight_data.id_score_densities(file_manager, dataset_id)
        kinds = data.select(pl.col('kind').n_unique()).collect().item()
        return dict(data=data, **score_distribution.settings(kinds, **_svg('identification-scores')))
    return _component(SeriesPlot, comp_name, file_manager, dataset_id, ['protein_dfs'], build)


# Ion type colours of the fragment error plot: prefix ions blue, suffix ions red
FRAGMENT_ION_COLORS = fragment_errors.ION_COLORS


def _fragment_error_plot(comp_name, file_manager, dataset_id):
    """Mass error of the matched fragments of the selected proteoform; a click
    selects the residue next to the fragment's bond."""
    return _component(
        SeriesPlot, comp_name, file_manager, dataset_id,
        ['sequence_data', 'protein_dfs', 'scan_table', 'sequence_view', 'settings'],
        lambda: dict(
            data=validation_data.fragment_errors(file_manager, dataset_id),
            **fragment_errors.settings(**_svg('fragment-errors')),
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
            **isotope_envelope.settings(scan_identifier, styling=spectrum_styling()),
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
            **precursor_xic.settings(scan_identifier, **_svg('precursor-xic')),
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
            **chromatogram.settings(list(qc_data.CHROMATOGRAM_SERIES.values()), level, **_svg('chromatogram')),
        ),
        required=['scan_table', 'anno_spectrum'],
    )


# Precursor spectrum links per dataset, keyed by the scan table file and its mtime
_SURVEY_LINKS = {}


def _survey_links(file_manager, dataset_id):
    path = insight_data.result_path(file_manager, dataset_id, 'scan_table')
    if path is None:
        raise MissingResult('scan_table')
    key = (str(path), path.stat().st_mtime_ns)
    if key not in _SURVEY_LINKS:
        _SURVEY_LINKS.clear()
        _SURVEY_LINKS[key] = linked_spectra.survey_links(file_manager, dataset_id)
    return _SURVEY_LINKS[key]


def _follow_selected_spectrum(file_manager, dataset_id, state_manager, scan_identifier):
    """Select the precursor spectrum and precursor mass of the selected
    spectrum, once per newly selected spectrum, so a mass clicked in the
    precursor spectrum stays selected until another spectrum is."""
    links = _survey_links(file_manager, dataset_id)
    if links['survey_idx'].is_null().all():
        raise MissingResult('no MS1 spectra')
    scan = state_manager.get_selection(scan_identifier)
    if scan == state_manager.get_selection(SURVEY_FOR) and state_manager.get_selection(SURVEY) is not None:
        return
    survey, survey_mass = linked_spectra.survey_of(links, scan)
    for identifier, value in ((SURVEY, survey), (SURVEY_MASS, survey_mass), (SURVEY_FOR, scan)):
        if value is None:
            state_manager.clear_selection(identifier)
        else:
            state_manager.set_selection(identifier, int(value))


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


# Isotope errors a sequence view tries by default when matching deconvolved
# masses: FLASHDeconv can place a monoisotopic mass one isotope off
ISOTOPE_ERRORS = 1


def _sequence_view_options(component):
    """Options of a sequence view, each when the installed openms-insight
    supports it: the residue menus offer the whole modification library, the
    shown proteoform can be copied as ProForma, the "Change sequence" dialog
    reads ProForma, fragment and precursor masses may be off by up to
    ISOTOPE_ERRORS isotopes (changeable in the view), the matched fragment
    ions are marked on the spectra, and proton loss/addition is offered only
    with electron or photo fragmentation ion series (c/z or radical series),
    where a hydrogen moves between the fragments. The observed masses no
    fragment ion explains are listed below the fragment table, each with the
    nearest fragment ion."""
    if FRAGMENT_MARKS:
        component.set_fragment_matches_identifier(FRAGMENT_MATCHES_IDENTIFIER)
    if hasattr(component, 'set_variable_modification_library'):
        component.set_variable_modification_library(_modification_picker_library())
    if hasattr(component, 'set_proforma_copy'):
        component.set_proforma_copy()
    if hasattr(component, 'set_proforma_input'):
        component.set_proforma_input()
    if hasattr(component, 'set_isotope_errors'):
        component.set_isotope_errors(ISOTOPE_ERRORS)
    if hasattr(component, 'set_method_ion_options'):
        component.set_method_ion_options()
    if hasattr(component, 'set_unmatched_masses'):
        component.set_unmatched_masses()


@st.cache_data(max_entries=4, show_spinner=False)
def _sequence_data(sequence, mod_ids):
    fixed_mods = {'fixed_modifications': list(mod_ids)}
    data = getFragmentDataFromSeq(sequence, fixed_mods=fixed_mods)
    # Plain floats: the sequence is part of the component configuration
    return json.loads(json.dumps(data, default=float))


def _search_entry(file_manager, dataset_id, proteoform):
    """The sequence view entry of a proteoform as the search reported it, or None."""
    if proteoform is None:
        return None
    rows = (
        insight_data.proteoform_sequences(file_manager, dataset_id)
        .filter(pl.col('proteoform_index') == int(proteoform))
        .head(1)
        .collect()
        .to_dicts()
    )
    return rows[0] if rows else None


@st.cache_data(max_entries=16, show_spinner=False)
def _edited_entry_data(edit_json, search_json):
    edit = json.loads(edit_json)
    search = json.loads(search_json) if search_json else None
    data = proteoform_edits.entry(
        edit, search,
        lambda sequence: getFragmentDataFromSeq(sequence, fixed_mods={'fixed_modifications': []}),
    )
    # Plain floats: the entry is part of the component configuration
    return json.loads(json.dumps(data, default=float))


def _edited_proteoform(file_manager, dataset_id, state_manager):
    """The selected proteoform as the user edited it in the sequence view
    (applying the edit the view wrote last), or None while it shows the
    search result. Returns (proteoform index, entry, sequence changed)."""
    proteoform = state_manager.get_selection('protein')
    edits = proteoform_edits.load(file_manager, dataset_id)
    edit = state_manager.get_selection(EDIT_IDENTIFIER)
    search = None
    # The edit stays in the selection (every view holds it and would send it
    # back if cleared): its nonce is applied once per session, to the
    # proteoform shown when it was made, so it never reaches another
    # proteoform or dataset
    nonce = edit.get('nonce') if isinstance(edit, dict) else None
    fresh = nonce is not None and nonce != st.session_state.get('_proteoform_edit_nonce')
    if fresh:
        st.session_state['_proteoform_edit_nonce'] = nonce
    if fresh and nonce != edits.get(proteoform_edits.NONCE):
        search = _search_entry(file_manager, dataset_id, proteoform)
        new = proteoform_edits.apply(edits, edit, proteoform, search)
        if new is not None:
            file_manager.store_data(dataset_id, proteoform_edits.RESULT, new)
            edits = new
    if proteoform is None or int(proteoform) not in edits:
        return None
    if search is None:
        search = _search_entry(file_manager, dataset_id, proteoform)
    edit = edits[int(proteoform)]
    search_json = json.dumps(
        {k: search.get(k) for k in ('sequence', 'observed_mass', 'coverage', 'maxCoverage')}
        if search else None, default=float,
    )
    try:
        data = _edited_entry_data(json.dumps(edit, sort_keys=True), search_json)
    except Exception as error:  # pyOpenMS refuses the residues
        st.warning(f'The edited sequence could not be used: {error}')
        return None
    changed = search is None or list(search.get('sequence') or []) != list(edit['sequence'])
    return int(proteoform), data, changed


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


def _set_selected_tag(component, file_manager, dataset_id, state_manager):
    tag = _selected_tag(file_manager, dataset_id, state_manager)
    # Tags FLASHTnT could not place on the proteoform have no positions
    if tag is None or tag['StartPos'] is None or tag['EndPos'] is None:
        component.set_tag(None)
    else:
        component.set_tag(int(tag['StartPos']), int(tag['EndPos']))


def _edited_sequence_view(comp_name, tool, fm, ds, state_manager, proteoform, entry, changed):
    """The identification sequence view of a proteoform the user edited: its
    spectrum matched against the edited sequence and modifications, with the
    button that resets it to the search result."""
    name = f'{comp_name}-{tool}-edited-{proteoform}'
    sources = ['sequence_data', 'protein_dfs', 'scan_table', 'sequence_view', 'settings',
               'out_deconv_mzML']
    component = _component(
        ProteinSequenceView, name, fm, ds, sources,
        lambda: dict(
            # Only this proteoform's spectrum: the cache is rebuilt per edit
            data=insight_data.with_fragmentation(
                insight_data.proteoform_sequences(fm, ds)
                .filter(pl.col('proteoform_index') == proteoform),
                _scan_methods(fm, ds), 'scan_idx',
            ),
            filters={'protein': 'proteoform_index'}, sequence=entry,
            mass_identifier='mass', aa_identifier='aa', fragment_identifier='fragment',
            fragmentation_column='fragmentation',
            settings=_identification_settings(fm, ds),
            entry_title='Proteoform', observed_mass_label='Observed proteoform mass',
            **SEQUENCE_VIEW_OPTIONS, **_spectrum_key('scan_idx'),
        ),
        extra=entry, required=['sequence_data', 'protein_dfs', 'scan_table', 'sequence_view', 'settings'],
    )
    # Earlier edits of this proteoform are not shown again
    keep = f'{name}-{_stamp(fm, ds, sources, entry)}'
    for old in insight_cache_dir(fm.cache_path, ds).glob(f'{name}-*'):
        if old.name != keep:
            shutil.rmtree(old, ignore_errors=True)
    if changed:
        component.set_tag(None)
    else:
        _set_selected_tag(component, fm, ds, state_manager)
    component.set_edit_identifier(EDIT_IDENTIFIER)
    _sequence_view_options(component)
    return component


# A mass shift picked in the Mass shifts component narrows the protein,
# proteoform and PrSM tables to the rows carrying it, besides their protein →
# proteoform filters; with no shift picked they show every row. Needs an
# openms-insight with ALL (list-column filters); older builds leave it out.
_ALL_ROWS = getattr(openms_insight, 'ALL', None)


def _hierarchy_table(comp_name, file_manager, dataset_id):
    """One level of the protein → proteoform → PrSM tables and its height.
    Selecting a proteoform or PrSM selects its (best) PrSM in the other
    cells, as a protein table row does. The tables share one height: that of
    the level with the most rows it can show at once (all proteins, the
    proteoforms of the protein with the most, the PrSMs of the proteoform
    with the most), up to TABLE_HEIGHT."""
    hits = hierarchy_data.prsms(file_manager, dataset_id)
    with_qvalues = hierarchy_data.has_qvalues(hits)

    def with_shifts(rows, key):
        # The mass shifts column and the shift ids the selected mass shift
        # filters on (empty for a row without one); built with the cache only
        groups, _, shifts = mass_shifts.summary(mass_shifts.prsm_shifts(hits))
        if key not in shifts.columns:
            shifts = shifts.join(hits.select('index', key), on='index', how='left')
        return (
            rows.join(mass_shifts.shift_labels(groups, shifts, key), on=key, how='left')
            .join(mass_shifts.shift_ids(shifts, key), on=key, how='left')
            .with_columns(pl.col('shift_ids').fill_null(pl.lit([], dtype=pl.List(pl.Int64))))
        )

    figures = hierarchy_data.key_figures(hits)
    shared = dict(with_qvalues=with_qvalues, table_cls=Table, protein_set=go_terms.RUN_SET, all_rows=_ALL_ROWS)
    if comp_name == 'protein_group_table':
        build = lambda: dict(
            data=with_shifts(hierarchy_data.protein_groups_from(hits), 'protein_group').lazy(),
            **protein_groups.settings(hierarchy_data.level_title('Proteins', figures['proteins']), **shared),
        )
    elif comp_name == 'proteoform_table':
        build = lambda: dict(
            data=with_shifts(hierarchy_data.proteoforms_from(hits), 'proteoform_id').lazy(),
            **proteoforms.settings(hierarchy_data.level_title('Proteoforms', figures['proteoforms']), **shared),
        )
    else:
        build = lambda: dict(
            data=with_shifts(hierarchy_data.prsm_rows_from(hits), 'index').lazy(),
            **proteoform_prsms.settings(hierarchy_data.level_title('PrSMs', figures['prsms']), **shared),
        )
    # All three are as tall as the one with the most rows, so that tables
    # sitting side by side in a layout row line up
    height = _hierarchy.height(max(hierarchy_data.most_rows(hits).values()))
    table = _component(
        Table, comp_name, file_manager, dataset_id, ['protein_dfs', 'scan_table'], build,
    )
    return table, height


def _protein_coverage(comp_name, file_manager, dataset_id, state_manager):
    """The proteoforms of the protein picked in the protein table, drawn
    where they lie on its sequence; a click picks the proteoform."""
    if SequenceCoverage is None:
        raise MissingResult('a newer openms-insight')
    hits = hierarchy_data.prsms(file_manager, dataset_id)
    if not hits.height:
        raise MissingResult('PrSMs')
    track = file_manager.get_display_name(dataset_id) or str(dataset_id)
    component = _component(
        SequenceCoverage, comp_name, file_manager, dataset_id, ['protein_dfs', 'scan_table'],
        lambda: dict(
            data=hierarchy_data.coverage_rows(hits, track).lazy(),
            **sequence_coverage.settings(
                '', count_label='PrSMs', match='proteoform', reference_filter='protein_group',
                title='Sequence coverage', track_order=[track],
            ),
        ),
        extra=track,
    )
    protein = hierarchy_data.protein_of(hits, state_manager.get_selection('protein_group'))
    component.set_reference(*(protein[::-1] if protein else ('', 'Pick a protein in the protein table')))
    return component


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
            lambda: dict(data=insight_data.scan_table(fm, ds), **scan_table.settings()),
        ), {'height': TABLE_HEIGHT}

    if comp_name == 'reporter_table':
        plex, _ = _plex(fm, ds)
        return _component(
            Table, comp_name, fm, ds, ['isobaric_dfs', 'scan_table'],
            lambda: dict(
                data=insight_data.reporter_table(fm, ds),
                **reporter_table.settings(isobaric.label(plex), fm.result_exists(ds, 'protein_dfs')),
            ),
            required=['isobaric_dfs', 'scan_table'],
        ), {'height': TABLE_HEIGHT}

    if comp_name == 'reporter_ions':
        plex, channels = _plex(fm, ds)
        mzs = dict(isobaric.PLEXES.get(plex, ('', []))[1])
        channel_mzs = [mzs.get(c, 0.0) for c in channels]
        data_args, markers = reporter_plot_args(channels, channel_mzs, X_MARKERS)
        return _component(
            LinePlot, comp_name, fm, ds, ['isobaric_dfs', 'scan_table'],
            lambda: dict(
                **data_args(insight_data.reporter_ions(fm, ds, channels, channel_mzs)),
                **reporter_ions.settings(),
            ),
            render=markers,
        ), {}

    if comp_name == 'mass_table':
        return _component(
            Table, comp_name, fm, ds, ['mass_table'],
            lambda: dict(data=insight_data.mass_table(fm, ds), **mass_table.settings()),
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
                **spectrum.spectrum_settings(view, scan_identifier, styling=spectrum_styling(), **_svg('lineplot')),
            ),
        )
        # The fragment ions the sequence view matched on this spectrum
        if FRAGMENT_MARKS:
            component.set_fragment_matches_identifier(FRAGMENT_MATCHES_IDENTIFIER)
        # The raw peaks of the selected mass get its theoretical isotope envelope
        if view in ('annotated', 'ladder') and hasattr(component, 'set_envelope'):
            mass = state_manager.get_selection('mass')
            component.set_envelope(mass, envelope_overlay.selected_mass_envelope(
                fm, ds, state_manager.get_selection(scan_identifier), mass,
            ))
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

    if comp_name == 'survey_spectrum':
        # FLASHTnT spectra follow the spectrum of the selected proteoform
        scan_identifier = 'pscan' if tool == 'flashtnt' else 'scan'
        _follow_selected_spectrum(fm, ds, state_manager, scan_identifier)
        component = _component(
            DeconvolvedSpectrum, f'{comp_name}-{tool}', fm, ds, ['combined_spectrum', 'scan_table'],
            lambda: dict(
                data=insight_data.spectra(fm, ds, 'combined_spectrum'),
                **precursor_spectrum.settings(SURVEY, SURVEY_MASS, spectrum_styling(), **_svg('precursor-spectrum')),
            ),
        )
        # The m/z range isolated for the selected spectrum (none for an MS1
        # spectrum); left out with an openms-insight release without it
        if hasattr(component, 'set_isolation_window'):
            window = linked_spectra.isolation_window_of(
                fm, ds, state_manager.get_selection(scan_identifier),
            )
            component.set_isolation_window(*(window or (None, None, None)))
        return component, {}

    if comp_name == '3D_SN_plot':
        return _component(
            Plot3D, comp_name, fm, ds,
            ['threedim_SN_plot', 'scan_table', 'deconv_spectrum'],
            lambda: dict(data=insight_data.signal_peaks(fm, ds), **signal_noise_3d.settings(**_svg('3d-plot'))),
        ), {'height': 800}

    if comp_name == 'ms1_3d_map':
        # The MS1 heatmap opened in 3D (layouts saved before the two were one)
        if PeakMap is not None:
            return _heatmap('ms1_raw_heatmap', fm, ds, view='3d'), {'height': 600}
        return _ms1_3d_map(comp_name, fm, ds, tool, state_manager), {'height': 600}

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

    if comp_name == 'protein_coverage':
        return _protein_coverage(comp_name, fm, ds, state_manager), {'height': 420}

    if comp_name == 'protein_table':
        return _component(
            Table, comp_name, fm, ds,
            ['protein_dfs', 'scan_table', 'sequence_data', 'sequence_view', 'settings', 'spec2_tsv'],
            lambda: dict(
                data=_proteins_with_characterization(fm, ds).join(
                    prsm_overview.protein_overview(fm, ds), on='index', how='left', maintain_order='left',
                ),
                **prsm_table.settings(
                    with_qvalues=hierarchy_data.has_qvalues(insight_data.proteins(fm, ds).collect()),
                    table_cls=Table, protein_set=go_terms.RUN_SET,
                ),
            ),
            required=['protein_dfs', 'scan_table'],
        ), {'height': TABLE_HEIGHT}

    if comp_name == 'prsm_mass_table':
        return _component(
            Table, comp_name, fm, ds,
            ['sequence_data', 'protein_dfs', 'scan_table', 'sequence_view', 'mass_table', 'settings'],
            lambda: dict(data=prsm_data.prsm_masses(fm, ds), **prsm_mass_table.settings()),
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
                **tag_table.settings(),
            ),
        ), {'height': TABLE_HEIGHT}

    if comp_name == 'sequence_view':
        # The proteoform of an identified run, else the sequence set in the view
        if tool == 'flashtnt' or (tool == 'topdown' and fm.result_exists(ds, 'protein_dfs')):
            # Editing needs an openms-insight whose view writes edits
            editable = hasattr(ProteinSequenceView, 'set_edit_identifier')
            edited = _edited_proteoform(fm, ds, state_manager) if editable else None
            if edited is not None:
                return _edited_sequence_view(comp_name, tool, fm, ds, state_manager, *edited), {}
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
                    **sequence_view.identified_settings(
                        _identification_settings(fm, ds), **_SEQUENCE_VIEW_NAMES,
                        spectrum_key='scan_idx' if FRAGMENT_MARKS else None,
                    ),
                ),
                required=['sequence_data', 'protein_dfs', 'scan_table', 'sequence_view', 'settings'],
            )
            _set_selected_tag(component, fm, ds, state_manager)
            _sequence_view_options(component)
            if hasattr(component, 'set_edit_identifier'):
                component.set_edit_identifier(EDIT_IDENTIFIER)
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
                **sequence_view.sequence_settings(
                    sequence, **_SEQUENCE_VIEW_NAMES, spectrum_key='index' if FRAGMENT_MARKS else None,
                ),
            ),
            extra=sequence, required=['sequence_view'],
        )
        # The "Change sequence" dialog also picks the fixed modifications
        if hasattr(component, 'set_modification_picker'):
            mod_ids = _fixed_modification_ids(_stored_sequence(fm), state_manager)
            component.set_modification_picker(
                _modification_picker_library(), FIXED_MODS_IDENTIFIER, mod_ids
            )
        _sequence_view_options(component)
        return component, {}

    if comp_name == 'structure_view':
        component = _component(
            ProteinStructureView, f'{comp_name}-{tool}', fm, ds,
            ['sequence_data', 'protein_dfs', 'scan_table', 'sequence_view', 'settings'],
            lambda: dict(
                data=insight_data.proteoform_sequences(fm, ds),
                **structure_view.settings(fm.get_results(ds, ['settings'])['settings']),
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


# (name, data, results it is built from); titles and axes are the library's
RUN_QC_PLOTS = [
    (name, builder, sources, *run_qc.PLOTS[name]) for name, builder, sources in (
        ('qc_qscore', qc_data.qscore_distribution, ['mass_table', 'scan_table']),
        ('qc_charge', qc_data.charge_distribution, ['mass_table', 'scan_table']),
        ('qc_mass', qc_data.mass_distribution, ['mass_table', 'scan_table']),
        ('qc_masses_per_spectrum', qc_data.masses_per_spectrum, ['scan_table']),
    )
]
ID_QC_PLOTS = [
    (name, builder, sources, *identification_qc.PLOTS[name]) for name, builder, sources in (
        ('qc_identifications', qc_data.identifications_over_time, ['protein_dfs', 'scan_table']),
        ('qc_coverage', qc_data.sequence_coverage, ['protein_dfs', 'scan_table']),
        ('qc_fragment_error', qc_data.fragment_error_over_time,
         ['sequence_data', 'protein_dfs', 'scan_table', 'sequence_view', 'settings']),
    )
]


def _render_qc(plots, file_manager, dataset_id, tool, state_manager, key, per_row=2, leading=()):
    """A panel of run-level plots, `per_row` per row (one when the panel
    shares its row of the layout, so no plot gets squeezed); plots whose
    results are missing are left out. `leading` builds components drawn
    first, among the panel's plots."""
    components = []
    for build_component in leading:
        try:
            components.append(build_component())
        except MissingResult:
            pass
    for name, builder, sources, title, x_label, y_label, mode, y_range in plots:
        def build(builder=builder, name=name, title=title, x_label=x_label,
                  y_label=y_label, mode=mode, y_range=y_range):
            data = builder(file_manager, dataset_id)
            # A legend for one series only repeats the axis title
            series = data.lazy().select(pl.col('series').n_unique()).collect().item()
            return dict(data=data, **run_qc.plot_settings(
                title, x_label, y_label, mode, y_range, series_count=series, name=name, **_svg(name),
            ))
        try:
            components.append(_component(
                SeriesPlot, f'{name}-{tool}', file_manager, dataset_id, sources, build,
                required=[source for source in sources if source != 'settings'],
            ))
        except MissingResult:
            continue
    if not components:
        st.info('This dataset has no data for the run QC plots.')
        return
    for row_start in range(0, len(components), per_row):
        columns = st.columns(per_row)
        for offset, component in enumerate(components[row_start:row_start + per_row]):
            with columns[offset]:
                component(
                    key=f'{key}_{row_start + offset}', state_manager=state_manager,
                    height=320 if per_row > 1 else 280, help=teaching.help(),
                )


def _render_mass_shifts(file_manager, dataset_id, state_manager, key):
    """The mass shifts FLASHTnT placed on the identified proteoforms: a
    histogram and the recurring shifts with the Unimod modifications they
    match, linked to each other. A picked shift narrows the protein,
    proteoform and PrSM tables to the rows carrying it; clicking it again in
    the table unpicks it."""
    sources = ['protein_dfs', 'scan_table']
    try:
        for name_tag in sources:
            if insight_data.result_path(file_manager, dataset_id, name_tag) is None:
                raise MissingResult(name_tag)
    except MissingResult as e:
        st.info(f'This dataset has no data for this view ({e}).')
        return
    hits = hierarchy_data.prsms(file_manager, dataset_id)
    groups, histogram, prsms = mass_shifts.summary(mass_shifts.prsm_shifts(hits))
    if not groups.height:
        st.info('No identified proteoform of this run carries a mass shift.')
        return
    shifted = prsms['index'].n_unique()
    columns = st.columns(3)
    columns[0].metric('PrSMs with a mass shift', f'{shifted:,} of {hits.height:,}')
    columns[1].metric('Shifts', f'{groups.height:,}', help='Groups of shifts within '
                      f'{mass_shifts.CLUSTER_GAP * 1000:g} mDa of each other.')
    columns[2].metric(
        'Unexplained', f"{(groups['interpretation'] == '').sum():,}",
        help='Shifts with no Unimod modification of that mass, '
             'also one isotope peak off.',
    )
    plot = _component(
        LinePlot, 'mass_shift_histogram', file_manager, dataset_id, sources,
        lambda: dict(
            data=histogram.lazy(), **mass_shift_histogram.settings(),
        ),
    )
    shift_table = _component(
        Table, 'mass_shift_table', file_manager, dataset_id, sources,
        lambda: dict(
            data=groups.lazy(), **mass_shift_table.settings(mass_shifts.MATCH_TOLERANCE),
        ),
    )
    plot(key=f'{key}_histogram', state_manager=state_manager, height=360, help=teaching.help())
    shift_table(key=f'{key}_shifts', state_manager=state_manager,
                height=mass_shift_table.height(groups.height, TABLE_HEIGHT))


INSPECT_LABEL = 'Inspect this spectrum in Match fragments'


def _inspect_button(component, comp_name, dataset_id, on_inspect):
    """A button per row of a spectrum table that hands the row's spectrum
    (and proteoform) to ``on_inspect(dataset_id, scan_idx, proteoform_index)``.
    An openms-insight release without row actions shows no button."""
    from src.render import inspect_scan
    if on_inspect is None or comp_name not in inspect_scan.TABLES or not hasattr(component, 'on_row_action'):
        return
    component.on_row_action(
        'inspect', lambda row: on_inspect(dataset_id, *inspect_scan.from_row(comp_name, row)),
        icon='mdi-puzzle-outline', label=INSPECT_LABEL,
    )


def _render_cell(comp_name, file_manager, dataset_id, tool, state_manager, key, on_missing=None,
                 row=None, on_inspect=None):
    with teaching.drawing(comp_name):
        _draw_cell(comp_name, file_manager, dataset_id, tool, state_manager, key, on_missing, row, on_inspect)


def _draw_cell(comp_name, file_manager, dataset_id, tool, state_manager, key, on_missing=None,
               row=None, on_inspect=None):
    row = row or [comp_name]
    per_row = 2 if len(row) == 1 else 1
    missing = []
    if (
        tool == 'topdown' and comp_name in IDENTIFICATION_ONLY
        and not file_manager.result_exists(dataset_id, 'protein_dfs')
    ):
        missing.append(IDENTIFICATION)
    if comp_name in ISOBARIC_ONLY and not file_manager.result_exists(dataset_id, 'isobaric_dfs'):
        missing.append(ISOBARIC)
    if missing:
        if on_missing is not None:
            from src.render.missing_analysis import ANALYSES
            missing_note(_COMPONENT_LABELS.get(comp_name, comp_name), missing, ANALYSES, on_missing,
                         f'{key}_missing', comp_name)
        elif IDENTIFICATION in missing:
            st.info('This run was deconvolved only. Run it with **Identify proteoforms** on to fill this view.')
        else:
            st.info('This run has no reporter ions. Run it with **Quantification: Isobaric labels** to fill this view.')
        return
    if comp_name == 'internal_fragment_map':
        st.info('The internal fragment map is not available yet.')
        return
    if comp_name == 'run_qc':
        _render_qc(RUN_QC_PLOTS, file_manager, dataset_id, tool, state_manager, key, per_row)
        return
    if comp_name == 'id_qc':
        # The identification score plot leads the panel, unless it has a
        # cell of its own next to it
        leading = [] if 'id_fdr_plot' in row else [
            lambda: _id_score_plot('id_fdr_plot', file_manager, dataset_id)]
        _render_qc(ID_QC_PLOTS, file_manager, dataset_id, tool, state_manager, key, per_row,
                   leading=leading)
        return
    if comp_name == 'mass_shifts':
        _render_mass_shifts(file_manager, dataset_id, state_manager, key)
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
        elif str(e) == 'a newer openms-insight':
            st.info('The sequence coverage map needs a newer openms-insight than the one installed.')
        else:
            st.info(f'This dataset has no data for this view ({e}).')
        return
    _inspect_button(component, comp_name, dataset_id, on_inspect)
    component(key=key, state_manager=state_manager, **{'help': teaching.help(), **render_args})


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
    grid_key='flash_viewer_grid', on_missing=None, on_inspect=None,
):
    """Draw one experiment. ``on_missing(missing, component)``: the button of the
    note a component draws when the run lacks the analysis it needs.
    ``on_inspect(dataset_id, scan_idx, proteoform_index)``: the inspect button
    of each row of the spectrum tables (``inspect_scan.TABLES``); returns True to
    rerun."""
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
                            f'{grid_key}_{row_index}_{col_index}', on_missing,
                            row=row, on_inspect=on_inspect,
                        )
