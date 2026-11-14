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
import re
import shutil
from pathlib import Path

import polars as pl
import streamlit as st
from openms_insight import ProteinSequenceView, StateManager, Table

try:
    from openms_insight.components.protein_sequence_view import FRAGMENTATION_METHODS
except ImportError:
    FRAGMENTATION_METHODS = {'CID/HCD': (), 'ETD/ECD': (), 'EThcD': (), 'UVPD': ()}
# ProForma 2.0 read by pyOpenMS (openms-insight with core.proforma)
try:
    from openms_insight.core import proforma as proforma_reader
    if not proforma_reader.available():
        proforma_reader = None
except ImportError:
    proforma_reader = None

from src import fragment_matcher as fm_input
from src.common.common import page_setup
from src.common.guide import page_intro
from src.render import insight_data, modifications
from src.render.render import (
    FIXED_MODS_IDENTIFIER, FRAGMENT_MATCHES_IDENTIFIER, SEQUENCE_VIEW_OPTIONS,
    _modification_picker_library, _scan_methods, _sequence_view_options,
)
from src.render.sequence import getFragmentDataFromSeq
from src.workflow.FileManager import FileManager

STATE_KEY = 'fragment_matcher_state'
# What the masses box holds
DECONVOLVED = 'Deconvolved masses'
MZ_VALUES = 'm/z values'
# The shown proteoform as ProForma, shared by the sequence box and the view
PROFORMA_IDENTIFIER = 'proforma'
# DNA-binding protein HU-alpha (E. coli), identified in that scan
EXAMPLE_SEQUENCE = 'MNKTQLIDVIAEKAELSKTQAKAALESTLAAITESLKEGDAVQLVGFGTFKVNHRAERTGRNPQTGKEIKIAAANVPAFVSGKALKDAVK'

params = page_setup()
workspace = Path(st.session_state['workspace'])
file_manager = FileManager(workspace, Path(workspace, 'cache'))

st.title('🧩 Match fragments')
page_intro(
    'Match fragment masses against a sequence. Give the sequence (plain or '
    'ProForma, or fetched by UniProt accession) and the fragment masses, or the '
    'm/z peaks of a fragment spectrum, which FLASHDeconv deconvolves, in the two '
    'boxes below. The sequence view follows them, and a change in the view '
    '(**Change sequence**, a residue\'s menu) shows in the sequence box. '
    '**Example** fills in a scan of the lysate example; files and scans of a run: '
    '**More ways to give masses**.',
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
        '**m/z values.** With **m/z values** chosen above the box, give the '
        'peaks of a fragment spectrum, m/z and intensity per line; FLASHDeconv '
        '(as in OpenMS) turns them into neutral monoisotopic masses. Without '
        'intensities it finds no isotope patterns.\n\n'
        '**Sequence.** A UniProt accession (e.g. `P0ACF0`, an isoform such as '
        '`P0ACF0-2`, or `sp|P0ACF0|DBHA_ECOLI`) above the sequence box fills in '
        'the sequence of that entry. Write modifications in ProForma 2.0, read '
        'by OpenMS: `K[Acetyl]`, `S[+79.9663]`, `[UNIMOD:21]`, an N-terminal '
        '`[Acetyl]-`, a fixed `<[Carbamidomethyl]@C>` in front, a range '
        '`(EVTS)[+79.97]` or anywhere in the sequence `[+42.0106]?SEQ`.\n\n'
        '**Fragmentation.** The method the fragments were made with picks the '
        'ion series the view matches.'
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
        request_proforma(data['sequence'])
    return data


def link_text(value):
    """The ProForma in the link: a string, or {'text', 'by', 'n'}."""
    if isinstance(value, dict):
        value = value.get('text')
    return value if isinstance(value, str) else ''


def request_proforma(text):
    """Ask the view for this proteoform: a numbered request, so an answer of
    the view to an earlier one never overwrites the sequence box."""
    n = st.session_state.get('matcher_proforma_n', 0) + 1
    st.session_state['matcher_proforma_n'] = n
    st.session_state['matcher_sequence_seen'] = text
    state().set_selection(PROFORMA_IDENTIFIER, {'text': text, 'by': 'app', 'n': n})


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
    name = 'Masses'
    given_as = st.session_state.get('matcher_given_as', fm_input.NEUTRAL)
    mass_type = st.session_state.get('matcher_mass_type', fm_input.MONOISOTOPIC)
    if st.session_state.get('matcher_values') == MZ_VALUES:
        # FLASHDeconv gives neutral monoisotopic masses
        intensities = parsed['intensities']
        if intensities and any(i is None for i in intensities):
            intensities = None
        deconvolved = fm_input.deconvolve_mz(parsed['masses'], intensities)
        if not deconvolved['masses']:
            st.session_state['matcher_mass_error'] = (
                f"FLASHDeconv found no isotope pattern in the {len(parsed['masses'])} m/z values. "
                'Give each peak with its intensity, and enough peaks per pattern.')
            return
        name = f"Deconvolved from {len(parsed['masses'])} m/z values"
        parsed = {**parsed, **deconvolved}
        given_as, mass_type = fm_input.NEUTRAL, fm_input.MONOISOTOPIC
    data = save(parsed['masses'], parsed['intensities'],
                st.session_state.get('matcher_precursor') or parsed['precursor_mass'],
                mass_type, given_as, name=name)
    # The box keeps what was typed
    st.session_state['matcher_text_stamp'] = fm_input.input_stamp(data)


def sequence_changed():
    """The sequence box changed: the view takes it (plain or ProForma)."""
    request_proforma(st.session_state.get('matcher_sequence', '').strip())


EXAMPLE_FILES = {
    DECONVOLVED: Path('example-data', 'fragment-matcher', 'scan4268_masses.tsv'),
    MZ_VALUES: Path('example-data', 'fragment-matcher', 'scan4268_mz.tsv'),
}
EXAMPLE_PRECURSOR = 9529.1849


def values_switched():
    """Masses and m/z values do not mix: switching empties the masses box."""
    fm_input.clear_input(workspace)
    st.session_state['matcher_text'] = ''
    st.session_state['matcher_text_stamp'] = ''
    st.session_state['matcher_mass_error'] = ''


def fetch_uniprot():
    """The sequence of the UniProt accession goes into the sequence box."""
    st.session_state['matcher_uniprot'] = None
    st.session_state['matcher_uniprot_error'] = ''
    text = (st.session_state.get('matcher_accession') or '').strip()
    if not text:
        return
    try:
        entry = fm_input.uniprot_entry(text)
    except fm_input.UniProtLookupError as error:
        st.session_state['matcher_uniprot_error'] = str(error)
        return
    st.session_state['matcher_uniprot'] = entry
    st.session_state['matcher_sequence'] = entry['sequence']
    sequence_changed()


def load_example():
    """Scan 4268 of the lysate example (HU-alpha) in the boxes: its deconvolved
    masses or its m/z peaks, as the masses box is set, with precursor and sequence."""
    values = st.session_state.get('matcher_values') or DECONVOLVED
    st.session_state['matcher_text'] = EXAMPLE_FILES[values].read_text(encoding='utf-8')
    st.session_state['matcher_precursor'] = EXAMPLE_PRECURSOR
    st.session_state['matcher_given_as'] = fm_input.NEUTRAL
    st.session_state['matcher_mass_type'] = fm_input.MONOISOTOPIC
    masses_changed()
    st.session_state['matcher_sequence'] = EXAMPLE_SEQUENCE
    sequence_changed()


def delete_all():
    """Empty both boxes and the precursor: no masses, no sequence, no modifications."""
    clear_masses()
    st.session_state['matcher_sequence'] = ''
    request_proforma('')
    manager = state()
    manager.set_selection('sequence', '')
    manager.set_selection(FIXED_MODS_IDENTIFIER, [])


def clear_masses():
    """Forget the masses and empty the fields they were given in (part of Delete all)."""
    fm_input.clear_input(workspace)
    for key in ('matcher_precursor', 'matcher_file', 'matcher_mass_error'):
        st.session_state.pop(key, None)


def run_masses(dataset_id, scan_number):
    """Deconvolved masses, precursor mass and fragmentation method (None when
    the run does not say) of a scan of a run, or None."""
    spectra = insight_data.sequence_view_spectra(file_manager, dataset_id)
    scans = insight_data.scan(file_manager, dataset_id, 'scan_table').select('index', 'Scan')
    rows = (
        spectra.join(scans, on='index')
        .filter(pl.col('Scan') == int(scan_number))
        .collect()
    )
    if rows.height == 0:
        return None
    index = int(rows['index'][0])
    methods = _scan_methods(file_manager, dataset_id)
    method = methods[index] if index < len(methods) else None
    if method not in FRAGMENTATION_METHODS:
        method = None
    return list(rows['MonoMass'][0]), rows['PrecursorMass'][0], method


# Input another page handed over (the Viewer's inspect button) is taken; else
# the first visit of a session starts empty, dropping an earlier session's input
handed_over = fm_input.claim_hand_over(workspace)
if not st.session_state.get('matcher_started'):
    if handed_over is None:
        fm_input.clear_input(workspace)
    st.session_state['matcher_started'] = True
# A scan of a run brings its fragmentation method (set before the widget is drawn)
if st.session_state.get('matcher_fragmentation_next'):
    st.session_state['matcher_fragmentation'] = st.session_state.pop('matcher_fragmentation_next')
stored = fm_input.load_input(workspace)
manager = state()
# A handed-over sequence goes to the sequence box and the view, as a new list's does
if handed_over and handed_over.get('sequence'):
    st.session_state['matcher_sequence'] = handed_over['sequence']
    manager.set_selection('sequence', handed_over['sequence'])
    request_proforma(handed_over['sequence'])

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
# (only its answer to the latest request, or a request of its dialog)
link = manager.get_selection(PROFORMA_IDENTIFIER)
if isinstance(link, dict):
    n = link.get('n') or 0
    latest = st.session_state.get('matcher_proforma_n', 0)
    if link.get('by') == 'dialog' and n > latest:
        st.session_state['matcher_proforma_n'] = n
        latest = n
    from_view = link_text(link)
    if (link.get('by') in ('view', 'dialog') and n >= latest
            and from_view != st.session_state.get('matcher_sequence_seen')):
        st.session_state['matcher_sequence'] = from_view
        st.session_state['matcher_sequence_seen'] = from_view

# Headings first, then what fills each box, then the boxes
sequence_heading, masses_heading = st.columns(2)
sequence_heading.markdown(
    '**Sequence** (plain or ProForma)',
    help='Type, paste or fetch it by UniProt accession (the field below). Applied '
             'when you leave the box or press Ctrl+Enter, and kept in step with the '
             'sequence view: a change there (Change sequence, the residue menus) '
             'shows here. Empty the box to clear the sequence.\n\n'
             'ProForma 2.0 (read by OpenMS): `K[Acetyl]`, `S[+79.9663]`, `[UNIMOD:21]`, '
             'terminal `[Acetyl]-SEQ`, fixed `<[Carbamidomethyl]@C>SEQ`, a range '
             '`(EVTS)[+79.97]` or anywhere in the sequence `[+42.0106]?SEQ`; the '
             'matched fragments narrow a range down. Position groups (`#g1`), labile '
             'modifications and mass gaps are not supported yet.',
)
masses_heading.markdown(
    '**Fragment m/z values**' if st.session_state.get('matcher_values') == MZ_VALUES
    else '**Fragment masses**',
    help='Deconvolved masses or the m/z peaks of a fragment spectrum, as chosen '
         'in the field below; one per line, optionally with its intensity. Applied '
         'when you leave the box or press Ctrl+Enter. Switching between the two '
         'empties the box. Mass type, a file or a scan of a run: under More ways '
         'to give masses.',
)
# Above the boxes, in a row of their own, so both boxes start level
accession_choice, values_choice = st.columns(2)
with accession_choice:
    with st.container(horizontal=True, vertical_alignment='bottom'):
        st.text_input(
            'UniProt accession', key='matcher_accession', on_change=fetch_uniprot,
            placeholder='UniProt accession, e.g. P0ACF0', label_visibility='collapsed',
            help='The sequence of a UniProt entry goes into the sequence box: an '
                 'accession, an isoform (P0ACF0-2) or a FASTA identifier '
                 '(sp|P0ACF0|DBHA_ECOLI). Enter or Get sequence fetches it.',
        )
        st.button('Get sequence', on_click=fetch_uniprot,
                  disabled=not (st.session_state.get('matcher_accession') or '').strip())
with values_choice:
    st.segmented_control(
        'Values are', [DECONVOLVED, MZ_VALUES], key='matcher_values', default=DECONVOLVED,
        on_change=values_switched, label_visibility='collapsed',
        help='Deconvolved masses: neutral fragment masses as FLASHDeconv or another '
             'deconvolution gives them. m/z values: the peaks of a fragment spectrum '
             '(m/z and intensity per line); FLASHDeconv turns them into masses.',
    )
sequence_box, masses_box = st.columns(2)
with masses_box:
    st.text_area(
        'Fragment m/z values' if st.session_state.get('matcher_values') == MZ_VALUES
        else 'Fragment masses',
        height=160, key='matcher_text', on_change=masses_changed,
        placeholder='One per line, optionally with its intensity, e.g.\n1043.512 25000\n2511.374 18000\nor: 1043.512; 2511.374; 980.1',
        label_visibility='collapsed',
    )
    if st.session_state.get('matcher_mass_error'):
        st.error(st.session_state['matcher_mass_error'])
with sequence_box:
    st.text_area(
        'Sequence (plain or ProForma)', height=160, key='matcher_sequence',
        on_change=sequence_changed,
        placeholder='e.g. MNKTQLIDVIAEK...\nor ProForma: [Acetyl]-MNK[Trimethyl]TQLIDVIAEK...\nor <[Carbamidomethyl]@C>PEPTCIDE',
        label_visibility='collapsed',
    )
    # Ranges and unlocalized shifts are read by pyOpenMS (openms-insight with
    # core.proforma): without it they are dropped, so say so instead of silently
    if proforma_reader is None and re.search(r'\(|\]\??\^?\d*\?', st.session_state.get('matcher_sequence') or ''):
        st.warning('Modification ranges `(EVTS)[+79.97]` and unlocalized shifts `[+42.01]?SEQ` need an '
                   'openms-insight that reads ProForma with pyOpenMS (with `openms_insight.core.proforma`) '
                   'and pyOpenMS 3.6 or later. This installation shows the sequence without them.')
    # The UniProt entry the sequence came from
    if st.session_state.get('matcher_uniprot_error'):
        st.error(st.session_state['matcher_uniprot_error'])
    elif entry_info := st.session_state.get('matcher_uniprot'):
        st.caption(
            f"[{entry_info['accession']}](https://www.uniprot.org/uniprotkb/{entry_info['accession']}/entry)"
            f" · {entry_info['entry']} · {entry_info['name']}"
            + (f" · *{entry_info['organism']}*" if entry_info['organism'] else '')
            + f" · {len(entry_info['sequence'])} residues"
        )

with st.container(horizontal=True, horizontal_alignment='distribute', vertical_alignment='bottom'):
    st.number_input(
        'Precursor mass (optional)', min_value=0.0, step=1.0, format='%.4f',
        key='matcher_precursor', on_change=masses_changed,
        help='The precursor mass, neutral or MH+ as the masses are given (under '
             'More ways to give masses); with m/z values always the neutral mass. '
             'An .msalign file, a scan of a run or Example fills it in.',
    )
    st.selectbox(
        'Fragmentation', list(FRAGMENTATION_METHODS), key='matcher_fragmentation',
        help='How the fragments were made: the sequence view matches its ion '
             'series (CID/HCD: b, y; ETD/ECD: c, z•; EThcD: b, c, y, z•; UVPD: '
             'a, b, c, x, y, z and their variants). The series can still be '
             'changed in the view\'s settings.',
    )
    st.button('Example', icon=':material/science:', on_click=load_example,
              help='Scan 4268 of the lysate example: its deconvolved masses or its '
                   'm/z peaks (as chosen above the masses), its precursor mass and '
                   'its protein, HU-alpha.')
    st.button('Delete all', icon=':material/delete_sweep:', on_click=delete_all,
              help='Empty the masses, the precursor and the sequence with its '
                   'modifications.')

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
                    # The scan's own fragmentation method, set before the next run draws it
                    if found[2]:
                        st.session_state['matcher_fragmentation_next'] = found[2]
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


def proforma_spans(sequence):
    """Modifications of the sequence box over several residues (a range,
    `(EVTS)[+79.97]`, or anywhere, `[+42.01]?SEQ`), as (start, end, mass,
    label), 0-based. One-residue modifications are the view's own."""
    text = link_text(manager.get_selection(PROFORMA_IDENTIFIER))
    if not text or not sequence or proforma_reader is None:
        return ()
    try:
        read = proforma_reader.read_proforma(text)
    except proforma_reader.ProFormaReadError:
        return ()
    if read['sequence'] != sequence:
        return ()
    return tuple((m['start'], m['end'], m['mass'], m['label'])
                 for m in read['modifications'] if m['start'] != m['end'])


@st.cache_data(max_entries=8, show_spinner=False)
def entry_for(sequence, mod_ids, average, spans=()):
    data = getFragmentDataFromSeq(
        sequence, fixed_mods={'fixed_modifications': list(mod_ids)}, average=average,
        # 1-based, as getFragmentDataFromSeq takes them
        modifications=[(s + 1, e + 1, m) for s, e, m, _ in spans] if spans else None,
    )
    if spans:
        data['theoretical_mass'] += sum(m for _, _, m, _ in spans)
        data['modifications'] = [
            {'start': s, 'end': e, 'mass_diff': m, 'labels': label} for s, e, m, label in spans
        ]
    return json.loads(json.dumps(data, default=float))


# No sequence: an empty entry, built here. pyOpenMS on an empty sequence
# warns and has crashed the server (access violation) now and then
EMPTY_ENTRY = {
    'sequence': [], 'theoretical_mass': 0.0, 'fixed_modifications': [],
    **{f'fragment_masses_{ion}': [] for ion in 'abcxyz'},
}

spans = proforma_spans(sequence)
if not sequence:
    entry = EMPTY_ENTRY
else:
    try:
        entry = entry_for(sequence, tuple(mod_ids), average, spans)
    except Exception as error:  # pyOpenMS refuses the residues
        st.error(f'The sequence could not be used: {error}')
        entry = EMPTY_ENTRY

# One cache per mass list and sequence; the earlier ones are not shown again
cache_path = Path(workspace, 'cache', 'insight', 'fragment-matcher')
sequence_stamp = hashlib.sha1(json.dumps([sequence, mod_ids, spans]).encode()).hexdigest()[:10]
fragmentation = st.session_state.get('matcher_fragmentation') or next(iter(FRAGMENTATION_METHODS))
method_stamp = hashlib.sha1(fragmentation.encode()).hexdigest()[:4]
cache_id = f'matcher-{fm_input.input_stamp(shown)}-{sequence_stamp}-{method_stamp}'
if cache_path.exists():
    for old in cache_path.glob('matcher-*'):
        if old.name != cache_id:
            shutil.rmtree(old, ignore_errors=True)
view = ProteinSequenceView(
    cache_id=cache_id,
    cache_path=str(cache_path),
    data=pl.LazyFrame([{**fm_input.spectrum_row(shown), 'Fragmentation': fragmentation}], schema={
        'index': pl.Int64, 'MonoMass': pl.List(pl.Float64), 'PrecursorMass': pl.Float64,
        'Fragmentation': pl.Utf8}),
    fragmentation_column='Fragmentation',
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
