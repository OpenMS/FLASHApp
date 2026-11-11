"""Match fragments: observed fragment masses against a sequence, without a
FLASHDeconv run.

The masses come from a pasted list, a file (text, tsv, csv or .msalign) or a
scan of a run in the workspace; the sequence is set in the sequence view's
own Change sequence dialog (plain sequence or ProForma), and modifications in
its residue menus. Monoisotopic or average masses, given as M or MH+. The
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
from openms_insight import ProteinSequenceView, StateManager

from src import fragment_matcher as fm_input
from src.common.common import page_setup
from src.common.guide import page_intro
from src.render import insight_data, modifications
from src.render.render import (
    FIXED_MODS_IDENTIFIER, SEQUENCE_VIEW_OPTIONS, _modification_picker_library,
    _sequence_view_options,
)
from src.render.sequence import getFragmentDataFromSeq
from src.workflow.FileManager import FileManager

STATE_KEY = 'fragment_matcher_state'
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
    'deconvolution: paste the masses, load them from a file, or take them from '
    'a scan of a run. Set the sequence in the view below (**Change sequence**, '
    'plain or ProForma) and add modifications by right-clicking a residue.',
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
        'neutral by subtracting one proton.'
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
    manager = state()
    if data.get('sequence'):
        manager.set_selection('sequence', data['sequence'])
    st.session_state.pop('matcher_seen_sequence', None)
    return data


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

with st.expander('Observed masses', expanded=stored is None):
    with st.container(horizontal=True):
        mass_type = st.radio(
            'Masses are', [fm_input.MONOISOTOPIC, fm_input.AVERAGE],
            format_func=str.capitalize, horizontal=True, key='matcher_mass_type',
            index=1 if stored and stored.get('mass_type') == fm_input.AVERAGE else 0,
        )
        given_as = st.radio(
            'Given as', [fm_input.NEUTRAL, fm_input.PROTONATED], horizontal=True,
            key='matcher_given_as',
            format_func=lambda v: 'Neutral (M)' if v == fm_input.NEUTRAL else 'Protonated (MH+)',
        )
    paste, upload, from_run = st.tabs(['Paste', 'From a file', 'From a run'])
    with paste:
        text = st.text_area(
            'Masses', height=180, key='matcher_text',
            placeholder='e.g.\n1043.512 25000\n2511.374 18000\n...\nor: 1043.512; 2511.374; 980.1',
        )
        precursor = st.number_input(
            'Precursor mass (optional)', min_value=0.0, value=0.0, step=1.0, format='%.4f',
            key='matcher_precursor',
        )
        if st.button('Match these masses', type='primary', disabled=not text.strip()):
            try:
                parsed = fm_input.parse_mass_list(text)
            except fm_input.MassListError as error:
                st.error(str(error))
            else:
                save(parsed['masses'], parsed['intensities'],
                     precursor or parsed['precursor_mass'], mass_type, given_as, name='Pasted masses')
                st.rerun()
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

if stored is None:
    st.info('Give the observed masses above; the sequence view then matches them.')
    st.stop()

average = stored.get('mass_type') == fm_input.AVERAGE
summary = f"**{stored.get('name') or 'Masses'}**: {len(stored['masses'])} {'average' if average else 'monoisotopic'} masses"
if stored.get('precursor_mass'):
    summary += f", precursor {stored['precursor_mass']:.2f} Da"
st.markdown(summary)

manager = state()
if stored.get('sequence') and manager.get_selection('sequence') is None:
    manager.set_selection('sequence', stored['sequence'])
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
cache_id = f'matcher-{fm_input.input_stamp(stored)}-{sequence_stamp}'
if cache_path.exists():
    for old in cache_path.glob('matcher-*'):
        if old.name != cache_id:
            shutil.rmtree(old, ignore_errors=True)
view = ProteinSequenceView(
    cache_id=cache_id,
    cache_path=str(cache_path),
    data=pl.LazyFrame([fm_input.spectrum_row(stored)]),
    filters={'matcher_spectrum': 'index'}, filter_defaults={'matcher_spectrum': 0},
    sequence=entry, sequence_identifier='sequence',
    settings=fm_input.DEFAULT_SETTINGS[stored.get('mass_type', fm_input.MONOISOTOPIC)],
    title='Sequence View',
    **{k: v for k, v in SEQUENCE_VIEW_OPTIONS.items() if k != 'annotation_title'},
)
if hasattr(view, 'set_modification_picker'):
    view.set_modification_picker(_modification_picker_library(), FIXED_MODS_IDENTIFIER, mod_ids)
_sequence_view_options(view)
view(key='fragment_matcher_view', state_manager=manager)
