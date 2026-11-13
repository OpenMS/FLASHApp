"""Match fragments: observed fragment masses against a sequence, without a
FLASHDeconv run.

Two boxes on top hold the masses and the sequence (plain or ProForma). The
masses box mirrors the stored input, which a file (text, tsv, csv or
.msalign) or a scan of a run can also set; the sequence box is linked both
ways with the sequence view (its Change sequence dialog with UniProt, its
residue menus) through the 'proforma' selection. Monoisotopic or
average masses, given as M or MH+. The
input is kept in <workspace>/fragment-matcher/input.json (see
src/fragment_matcher.py, which other pages use to open this one with masses
of their own).
"""

import hashlib
import json
import shutil
from pathlib import Path

import polars as pl
import streamlit as st
from openms_insight import ProteinSequenceView, StateManager, Table

from src import fragment_matcher as fm_input
from src.common.common import page_setup
from src.common.guide import page_intro
from src.render import insight_data, modifications
from src.render.render import (
    FIXED_MODS_IDENTIFIER, FRAGMENT_MATCHES_IDENTIFIER, SEQUENCE_VIEW_OPTIONS,
    _modification_picker_library,
    _sequence_view_options,
)
from src.render.sequence import getFragmentDataFromSeq
from src.workflow.FileManager import FileManager

STATE_KEY = 'fragment_matcher_state'
# The shown proteoform as ProForma, shared by the sequence box and the view
PROFORMA_IDENTIFIER = 'proforma'
EXAMPLE_RUN = 'example_lysate_identified'
EXAMPLE_SCAN = 4268
# DNA-binding protein HU-alpha (E. coli), identified in that scan
EXAMPLE_SEQUENCE = 'MNKTQLIDVIAEKAELSKTQAKAALESTLAAITESLKEGDAVQLVGFGTFKVNHRAERTGRNPQTGKEIKIAAANVPAFVSGKALKDAVK'

params = page_setup()
workspace = Path(st.session_state['workspace'])
file_manager = FileManager(workspace, Path(workspace, 'cache'))

st.title('🧩 Match fragments')
page_intro(
    'Match a list of fragment masses against a sequence, without running a '
    'deconvolution. Type or paste the masses and the sequence (plain or '
    'ProForma) in the two boxes below; the sequence view follows them, and a '
    'change in the view (**Change sequence** with UniProt, a residue\'s menu) '
    'shows in the sequence box. Files and scans of a run: **More ways to give masses**.',
    background=(
        '**Masses.** One mass per line, optionally followed by its intensity, '
        'separated by spaces, tabs, commas or semicolons; lines that are not '
        'numbers are skipped. Or all on one line, separated by semicolons '
        '(`1043.512; 2511.374; 980.1`). An `.msalign` file gives the peaks of its first '
        'spectrum and its precursor mass.\n\n'
        '**Monoisotopic or average.** Choose what the masses are. With average '
        'masses the sequence and its fragments are computed as average masses '
        'too, and the view starts with a tolerance of 0.5 Da; modifications '
        'added in the view shift them by their monoisotopic mass.\n\n'
        '**M or MH+.** Masses given as MH+ (singly protonated) are made '
        'neutral by subtracting one proton.\n\n'
        '**Sequence.** In **Change sequence**, a UniProt accession (e.g. '
        '`P0ACF0`, an isoform such as `P0ACF0-2`, or `sp|P0ACF0|DBHA_ECOLI`) '
        'fills in the sequence of that entry. Edit it there, or write '
        'modifications in ProForma: `K[Acetyl]`, `S[+79.9663]`, an N-terminal '
        '`[Acetyl]-` or a fixed `<[Carbamidomethyl]@C>` in front.'
    ),
)


def state():
    return StateManager(session_key=STATE_KEY)


def save(masses, intensities, precursor, mass_type, given_as, sequence=None, name=None):
    data = fm_input.store_input(
        workspace, fm_input.to_neutral(masses, given_as), intensities,
        fm_input.to_neutral([precursor], given_as)[0] if precursor else None,
        mass_type, sequence, name,
    )
    # A new list starts from its own sequence, if it brings one
    if data.get('sequence'):
        manager = state()
        manager.set_selection('sequence', data['sequence'])
        manager.set_selection(PROFORMA_IDENTIFIER, data['sequence'])
    return data


def mass_lines(data):
    """The masses of an input as the masses box shows them."""
    if not data:
        return ''
    intensities = data.get('intensities')
    return '\n'.join(
        f'{m:.4f}' + (f' {intensities[i]:g}' if intensities else '')
        for i, m in enumerate(data['masses'])
    )


def masses_changed():
    """The masses box (or its mass type, given as, precursor) changed: store it."""
    st.session_state['matcher_mass_error'] = ''
    text = st.session_state.get('matcher_text', '').strip()
    if not text:
        fm_input.clear_input(workspace)
        st.session_state['matcher_text_stamp'] = ''
        return
    try:
        parsed = fm_input.parse_mass_list(text)
    except fm_input.MassListError as error:
        st.session_state['matcher_mass_error'] = str(error)
        return
    data = save(parsed['masses'], parsed['intensities'],
                st.session_state.get('matcher_precursor') or parsed['precursor_mass'],
                st.session_state.get('matcher_mass_type', fm_input.MONOISOTOPIC),
                st.session_state.get('matcher_given_as', fm_input.NEUTRAL), name='Masses')
    # The box keeps what was typed
    st.session_state['matcher_text_stamp'] = fm_input.input_stamp(data)


def sequence_changed():
    """The sequence box changed: the view takes it (plain or ProForma)."""
    text = st.session_state.get('matcher_sequence', '').strip()
    st.session_state['matcher_sequence_seen'] = text
    state().set_selection(PROFORMA_IDENTIFIER, text)


def delete_all():
    """Empty both boxes and the precursor: no masses, no sequence, no modifications."""
    clear_masses()
    st.session_state['matcher_sequence'] = ''
    st.session_state['matcher_sequence_seen'] = ''
    manager = state()
    manager.set_selection(PROFORMA_IDENTIFIER, '')
    manager.set_selection('sequence', '')
    manager.set_selection(FIXED_MODS_IDENTIFIER, [])


def clear_masses():
    """Forget the masses and empty the fields they were given in (part of Delete all)."""
    fm_input.clear_input(workspace)
    for key in ('matcher_precursor', 'matcher_file', 'matcher_mass_error'):
        st.session_state.pop(key, None)


def run_masses(dataset_id, scan_number):
    """Deconvolved masses and precursor mass of a scan of a run, or None."""
    spectra = insight_data.sequence_view_spectra(file_manager, dataset_id)
    scans = insight_data.scan(file_manager, dataset_id, 'scan_table').select('index', 'Scan')
    rows = (
        spectra.join(scans, on='index')
        .filter(pl.col('Scan') == int(scan_number))
        .collect()
    )
    if rows.height == 0:
        return None
    return list(rows['MonoMass'][0]), rows['PrecursorMass'][0]


stored = fm_input.load_input(workspace)
manager = state()

# The masses box shows the stored masses whenever they change elsewhere (a
# file, a run, the example, Delete all); the user's own edit keeps its text
stamp = fm_input.input_stamp(stored) if stored else ''
if st.session_state.get('matcher_text_stamp') != stamp:
    st.session_state['matcher_text'] = mass_lines(stored)
    st.session_state['matcher_precursor'] = float((stored or {}).get('precursor_mass') or 0.0)
    st.session_state['matcher_text_stamp'] = stamp
    if stored:
        st.session_state['matcher_given_as'] = fm_input.NEUTRAL
# The sequence box follows the view (its dialog, residue menus)
from_view = manager.get_selection(PROFORMA_IDENTIFIER)
if from_view is not None and from_view != st.session_state.get('matcher_sequence_seen'):
    st.session_state['matcher_sequence'] = from_view
    st.session_state['matcher_sequence_seen'] = from_view

masses_box, sequence_box = st.columns(2)
with masses_box:
    st.text_area(
        'Fragment masses', height=160, key='matcher_text', on_change=masses_changed,
        placeholder='One per line, optionally with its intensity, e.g.\n1043.512 25000\n2511.374 18000\nor: 1043.512; 2511.374; 980.1',
        help='Applied when you leave the field or press Ctrl+Enter. Mass type, '
             'a file or a scan of a run: below, under More ways to give masses.',
    )
    with st.container(horizontal=True, vertical_alignment='bottom'):
        st.number_input(
            'Precursor mass (optional)', min_value=0.0, step=1.0, format='%.4f',
            key='matcher_precursor', on_change=masses_changed,
            help='Neutral or MH+ as the masses (Given as, below). A precursor in an '
                 '.msalign file or the scan of a run fills it in.',
        )
        st.button('Delete all', icon=':material/delete_sweep:', on_click=delete_all,
                  help='Empty the masses, the precursor and the sequence with its '
                       'modifications.')
    if st.session_state.get('matcher_mass_error'):
        st.error(st.session_state['matcher_mass_error'])
with sequence_box:
    st.text_area(
        'Sequence (plain or ProForma)', height=160, key='matcher_sequence',
        on_change=sequence_changed,
        placeholder='e.g. MNKTQLIDVIAEK...\nor ProForma: [Acetyl]-MNK[Trimethyl]TQLIDVIAEK...\nor <[Carbamidomethyl]@C>PEPTCIDE',
        help='Applied when you leave the field or press Ctrl+Enter, and kept in step '
             'with the sequence view: a change there (Change sequence, UniProt, the '
             'residue menus) shows here. Empty it to clear the sequence.',
    )

with st.expander('More ways to give masses'):
    with st.container(horizontal=True):
        mass_type = st.radio(
            'Masses are', [fm_input.MONOISOTOPIC, fm_input.AVERAGE],
            format_func=str.capitalize, horizontal=True, key='matcher_mass_type',
            index=1 if stored and stored.get('mass_type') == fm_input.AVERAGE else 0,
            on_change=masses_changed,
        )
        given_as = st.radio(
            'Given as', [fm_input.NEUTRAL, fm_input.PROTONATED], horizontal=True,
            key='matcher_given_as', on_change=masses_changed,
            format_func=lambda v: 'Neutral (M)' if v == fm_input.NEUTRAL else 'Protonated (MH+)',
        )
    upload, from_run = st.tabs(['From a file', 'From a run'])
    with upload:
        file = st.file_uploader(
            'Mass list', type=['txt', 'tsv', 'csv', 'msalign'], key='matcher_file',
        )
        if file is not None and st.button('Match the masses of this file', type='primary'):
            try:
                parsed = fm_input.parse_mass_list(file.getvalue())
            except fm_input.MassListError as error:
                st.error(str(error))
            else:
                save(parsed['masses'], parsed['intensities'], parsed['precursor_mass'],
                     mass_type, given_as, name=file.name)
                st.rerun()
    with from_run:
        runs = file_manager.get_results_list(['sequence_view'])
        if not runs:
            st.info('No run in this workspace has deconvolved spectra yet.')
        else:
            names = {r: file_manager.get_display_name(r) for r in runs}
            with st.container(horizontal=True, vertical_alignment='bottom'):
                run = st.selectbox('Run', runs, format_func=lambda r: names[r], key='matcher_run')
                scan_number = st.number_input('Scan number', min_value=0, step=1, key='matcher_scan')
                take = st.button('Match this scan', type='primary')
            if take:
                found = run_masses(run, scan_number)
                if found is None:
                    st.error(f'Scan {scan_number} has no deconvolved masses in this run.')
                else:
                    save(found[0], None, found[1], fm_input.MONOISOTOPIC, fm_input.NEUTRAL,
                         name=f'{names[run]}, scan {scan_number}')
                    st.rerun()
    if EXAMPLE_RUN in file_manager.get_results_list(['sequence_view']):
        if st.button(f'Load an example (scan {EXAMPLE_SCAN} of the lysate example, with HU-alpha)'):
            found = run_masses(EXAMPLE_RUN, EXAMPLE_SCAN)
            if found is not None:
                save(found[0], None, found[1], fm_input.MONOISOTOPIC, fm_input.NEUTRAL,
                     sequence=EXAMPLE_SEQUENCE, name=f'Lysate example, scan {EXAMPLE_SCAN}')
                st.rerun()

# The view is shown without masses too, so the sequence can be set first
shown = stored if stored is not None else fm_input.matcher_input([], mass_type=mass_type)
average = shown.get('mass_type') == fm_input.AVERAGE
if stored is not None:
    summary = f"**{stored.get('name') or 'Masses'}**: {len(stored['masses'])} {'average' if average else 'monoisotopic'} masses"
    if stored.get('precursor_mass'):
        summary += f", precursor {stored['precursor_mass']:.2f} Da"
    st.markdown(summary)

if shown.get('sequence') and manager.get_selection('sequence') is None:
    manager.set_selection('sequence', shown['sequence'])
sequence = manager.get_selection('sequence') or ''
mod_ids = modifications.fixed_ids({'fixed_modifications': []}, manager.get_selection(FIXED_MODS_IDENTIFIER))


@st.cache_data(max_entries=8, show_spinner=False)
def entry_for(sequence, mod_ids, average):
    data = getFragmentDataFromSeq(sequence, fixed_mods={'fixed_modifications': list(mod_ids)},
                                  average=average)
    return json.loads(json.dumps(data, default=float))


try:
    entry = entry_for(sequence, tuple(mod_ids), average)
except Exception as error:  # pyOpenMS refuses the residues
    st.error(f'The sequence could not be used: {error}')
    entry = entry_for('', (), average)

# One cache per mass list and sequence; the earlier ones are not shown again
cache_path = Path(workspace, 'cache', 'insight', 'fragment-matcher')
sequence_stamp = hashlib.sha1(json.dumps([sequence, mod_ids]).encode()).hexdigest()[:10]
cache_id = f'matcher-{fm_input.input_stamp(shown)}-{sequence_stamp}'
if cache_path.exists():
    for old in cache_path.glob('matcher-*'):
        if old.name != cache_id:
            shutil.rmtree(old, ignore_errors=True)
view = ProteinSequenceView(
    cache_id=cache_id,
    cache_path=str(cache_path),
    data=pl.LazyFrame([fm_input.spectrum_row(shown)], schema={
        'index': pl.Int64, 'MonoMass': pl.List(pl.Float64), 'PrecursorMass': pl.Float64}),
    filters={'matcher_spectrum': 'index'}, filter_defaults={'matcher_spectrum': 0},
    sequence=entry, sequence_identifier='sequence',
    settings=fm_input.DEFAULT_SETTINGS[shown.get('mass_type', fm_input.MONOISOTOPIC)],
    title='Sequence View',
    **{k: v for k, v in SEQUENCE_VIEW_OPTIONS.items() if k != 'annotation_title'},
)
if hasattr(view, 'set_modification_picker'):
    view.set_modification_picker(_modification_picker_library(), FIXED_MODS_IDENTIFIER, mod_ids)
_sequence_view_options(view)
if hasattr(view, 'set_fragment_matches_proteoform'):
    # The matched proteoform comes with its fragments, for the comparison below
    view.set_fragment_matches_identifier(FRAGMENT_MATCHES_IDENTIFIER)
    view.set_fragment_matches_proteoform()
# The Change sequence dialog can fetch a sequence from UniProt
if hasattr(view, 'set_uniprot_lookup'):
    view.set_uniprot_lookup()
# The sequence box and the view show the same proteoform
if hasattr(view, 'set_proforma_identifier'):
    view.set_proforma_identifier(PROFORMA_IDENTIFIER)
view(key='fragment_matcher_view', state_manager=manager)


HYPOTHESIS_IDENTIFIER = 'hypothesis'


def _fixed(precision):
    return {'formatter': 'fixed',
            'formatterParams': {'precision': precision, 'placeholders': {'NaN': ''}}}


HYPOTHESIS_COLUMNS = [
    {'field': 'best', 'title': '', 'width': 40, 'hozAlign': 'center', 'headerSort': False,
     'headerTooltip': 'The proteoform with the most matched fragments'},
    {'field': 'proforma', 'title': 'Proteoform', 'hozAlign': 'left', 'width': 260,
     'tooltip': True, 'headerTooltip': 'ProForma'},
    {'field': 'matched', 'title': 'Matched', 'sorter': 'number',
     'headerTooltip': 'Matched fragments, counted as the fragment table above counts them'},
    {'field': 'cleavages', 'title': 'Cleav. %', 'sorter': 'number', **_fixed(0),
     'headerTooltip': 'Bonds with a matched fragment on either side'},
    {'field': 'delta_da', 'title': 'Δ Da', 'sorter': 'number', **_fixed(2),
     'headerTooltip': 'Precursor mass minus theoretical mass, in Da'},
    {'field': 'delta_ppm', 'title': 'Δ ppm', 'sorter': 'number', **_fixed(1),
     'headerTooltip': 'Precursor mass minus theoretical mass, in ppm'},
    {'field': 'theoretical_mass', 'title': 'Mass', 'headerTooltip': 'Theoretical mass of the proteoform', 'sorter': 'number', **_fixed(2)},
    {'field': 'n_terminal', 'title': 'N', 'sorter': 'number', 'headerTooltip': 'N-terminal fragments'},
    {'field': 'c_terminal', 'title': 'C', 'sorter': 'number', 'headerTooltip': 'C-terminal fragments'},
    {'field': 'tolerance', 'title': 'Tol.', 'hozAlign': 'left',
     'headerTooltip': 'Fragment mass tolerance the proteoform was matched with'},
]


def _hypotheses_table(workspace, stored, rows, manager):
    """The proteoforms tried, as an Insight table; a picked row can be shown
    in the view again or removed."""
    best = fm_input.best_hypothesis(rows)
    has_precursor = any(r.get('delta_da') is not None for r in rows)
    frame = pl.DataFrame([
        {
            'id': i,
            'best': '★' if i == best else '',
            'matched': r['matched_fragments'],
            'cleavages': r['cleavages'],
            'delta_da': r['delta_da'],
            'delta_ppm': r['delta_ppm'],
            'theoretical_mass': r['theoretical_mass'],
            'n_terminal': r['n_terminal'],
            'c_terminal': r['c_terminal'],
            'tolerance': r['tolerance'],
            'proforma': r['proforma'],
        }
        for i, r in enumerate(rows)
    ], schema_overrides={'delta_da': pl.Float64, 'delta_ppm': pl.Float64})
    columns = [c for c in HYPOTHESIS_COLUMNS
               if has_precursor or c['field'] not in ('delta_da', 'delta_ppm')]
    digest = hashlib.sha1(
        json.dumps(rows, sort_keys=True, default=str).encode()).hexdigest()[:12]
    cache_path = Path(workspace, 'fragment-matcher', 'insight')
    cache_id = f'hypotheses-{digest}'
    # Tables of earlier lists are not shown again
    if cache_path.exists():
        for old in cache_path.glob('hypotheses-*'):
            if old.name != cache_id:
                shutil.rmtree(old, ignore_errors=True)
    table = Table(
        cache_id=cache_id, cache_path=str(cache_path), data=frame.lazy(),
        column_definitions=columns, index_field='id', title=f'Proteoforms tried ({len(rows)})',
        interactivity={HYPOTHESIS_IDENTIFIER: 'id'}, layout='fitDataFill',
        initial_sort=[{'column': 'matched', 'dir': 'desc'}], pagination=False,
        default_row=-1, default_row_selects=False, float_precision='float64',
    )
    table(key=f'matcher_hypotheses_{digest}', state_manager=manager, height=min(80 + 40 * len(rows), 360))
    picked = manager.get_selection(HYPOTHESIS_IDENTIFIER)
    picked = picked if isinstance(picked, int) and 0 <= picked < len(rows) else None
    with st.container(horizontal=True):
        st.button('Show in the view', icon=':material/visibility:', disabled=picked is None,
                  key='matcher_show_hypothesis', on_click=_show_hypothesis,
                  args=(rows[picked]['proforma'] if picked is not None else '',),
                  help='Show the picked proteoform in the sequence view and the sequence box')
        st.button('Remove', icon=':material/delete:', disabled=picked is None,
                  key='matcher_remove_hypotheses', on_click=_remove_hypothesis,
                  args=(stored, picked), help='Remove the picked proteoform from the table')


def _show_hypothesis(proforma):
    """Show a proteoform of the table in the view (as the sequence box does)."""
    st.session_state['matcher_sequence'] = proforma
    st.session_state['matcher_sequence_seen'] = proforma
    manager = state()
    manager.set_selection(PROFORMA_IDENTIFIER, proforma)
    manager.set_selection(HYPOTHESIS_IDENTIFIER, None)


def _remove_hypothesis(stored, index):
    fm_input.remove_hypotheses(workspace, stored, [index])
    state().set_selection(HYPOTHESIS_IDENTIFIER, None)


# Proteoforms tried against these masses
if hasattr(view, 'set_fragment_matches_proteoform'):
    st.subheader('Compare proteoforms')
    st.caption(
        'Change the sequence, its modifications or the ion types in the view, '
        'and add each proteoform you try: the table compares them against '
        'these masses. ★ marks the one with the most matched fragments; pick a '
        'row to show it in the view again or to remove it.'
    )
    current = fm_input.hypothesis(
        manager.get_selection(FRAGMENT_MATCHES_IDENTIFIER), stored.get('precursor_mass'),
    )
    rows = fm_input.load_hypotheses(workspace, stored)
    with st.container(horizontal=True, vertical_alignment='center'):
        if st.button('Add this proteoform', icon=':material/add:', type='primary',
                     disabled=current is None, key='matcher_add_hypothesis'):
            rows = fm_input.add_hypothesis(workspace, stored, current)
        if current is not None:
            st.caption(
                f"{current['matched_fragments']} matched fragments, "
                f"{current['cleavages']:.0f}% of the bonds"
            )
    if rows:
        _hypotheses_table(workspace, stored, rows, manager)
