"""The proteoform sequence of the Viewer's sequence view for a run that was not
identified: entered in a window opened from the Viewer, saved with its fixed
modifications to the workspace's ``'sequence'`` dataset."""

import re

import streamlit as st

from src.render import modifications
from src.render.layout_editor.dialogs import _own_file_manager

DATASET = 'sequence'
# Sites offered for a fixed modification: the termini, then the residues
MOD_SITES = [modifications.N_TERM, modifications.C_TERM] + list('ACDEFGHIKLMNPQRSTUVWY')
_VALID = re.compile('^[ac-ik-wyAC-IK-WYXx]+$')  # amino acid letters and X, no B, J or Z
_TEXT_KEY = 'sequence_text'
_MODS_KEY = 'fixed_modifications'
_LOADED_KEY = 'sequence_input_loaded'


def clean(text):
    """The sequence without white space, in upper case."""
    return ''.join((text or '').split()).upper()


def is_valid(text):
    """True for a sequence of amino acid letters (X allowed, no B, J, Z)."""
    seq = clean(text)
    return bool(seq) and bool(_VALID.match(seq))


def saved(file_manager):
    """(sequence, fixed modification ids) saved for the workspace, or None."""
    if not file_manager.result_exists(DATASET, 'sequence'):
        return None
    stored = file_manager.get_results(DATASET, 'sequence')['sequence']
    mod_ids = stored.get('fixed_modifications')
    if mod_ids is None:
        # Saved before the modification library: Cys/Met names
        mod_ids = modifications.legacy_ids(stored.get('fixed_mod_cysteine'), stored.get('fixed_mod_methionine'))
    return stored['input_sequence'], list(mod_ids)


def save(file_manager, sequence, mod_ids):
    file_manager.store_data(DATASET, 'sequence', {
        'input_sequence': clean(sequence),
        # Kept for readers of the earlier format; the library ids replace them
        'fixed_mod_cysteine': None,
        'fixed_mod_methionine': None,
        'fixed_modifications': list(mod_ids),
    })


def clear(file_manager):
    if file_manager.result_exists(DATASET, 'sequence'):
        file_manager.remove_results(DATASET)


def _mod_label(entry):
    label = f"{entry['id']}   {entry['mass']:+.4f} Da"
    if entry['unimod']:
        label += f"   {entry['unimod']}"
    if entry['description'] and entry['description'] != entry['name']:
        label += f"   {entry['description']}"
    return label


def _load(file_manager):
    """Fill the window's fields from the saved sequence when it opens."""
    if st.session_state.get(_LOADED_KEY):
        return
    current = saved(file_manager)
    st.session_state[_TEXT_KEY] = current[0] if current else ''
    st.session_state[_MODS_KEY] = current[1] if current else []
    st.session_state[_LOADED_KEY] = True


@st.dialog('Proteoform sequence', width='large')
def _dialog(file_manager):
    file_manager = _own_file_manager(file_manager)
    _load(file_manager)
    st.caption('Shown in the sequence view of runs that were not identified.')
    st.text_area('Sequence', key=_TEXT_KEY, height=120, placeholder='e.g. MKVLAAGIVGLLLAQ...')

    c1, c2, c3 = st.columns([1, 5, 1], vertical_alignment='bottom')
    site = c1.selectbox('Site', MOD_SITES, index=MOD_SITES.index('C'), key='fixed_mod_site')
    candidates = modifications.search(site=site)
    choice = c2.selectbox(
        'Fixed modification', candidates, index=None, format_func=_mod_label,
        placeholder='Search by name, Unimod accession or mass, e.g. Carbamidomethyl',
        key=f'fixed_mod_choice_{site}',
    )
    if c3.button('Add', disabled=choice is None, width='stretch'):
        if choice['id'] not in st.session_state[_MODS_KEY]:
            st.session_state[_MODS_KEY] = st.session_state[_MODS_KEY] + [choice['id']]

    sequence = clean(st.session_state.get(_TEXT_KEY))
    for mod_id in list(st.session_state[_MODS_KEY]):
        entry = modifications.by_id(mod_id)
        if entry is None:
            continue
        row = st.columns([5, 2, 1], vertical_alignment='center')
        sites = f", {len(modifications.sites_in(sequence, entry))} sites" if sequence else ''
        row[0].write(f"**{entry['site']}** · {entry['id']}")
        row[1].write(f"{entry['mass']:+.4f} Da{sites}")
        if row[2].button('Remove', key=f'remove_fixed_mod_{mod_id}'):
            st.session_state[_MODS_KEY] = [m for m in st.session_state[_MODS_KEY] if m != mod_id]
            st.rerun(scope='fragment')
    if sequence and is_valid(sequence):
        mass = modifications.modified_mass(sequence.replace('X', ''), st.session_state[_MODS_KEY])
        st.markdown(f'Theoretical mass: **{mass:.4f} Da**')

    with st.container(horizontal=True, horizontal_alignment='right'):
        if st.button('Remove sequence', disabled=saved(file_manager) is None):
            clear(file_manager)
            _close()
        if st.button('Save', type='primary'):
            if not sequence:
                clear(file_manager)
                _close()
            elif is_valid(sequence):
                save(file_manager, sequence, st.session_state[_MODS_KEY])
                _close()
            else:
                st.error('Only amino acid letters (and X) are allowed.')


def _close():
    st.session_state.pop(_LOADED_KEY, None)
    st.rerun()


def button(file_manager, key='sequence_button', label='Sequence'):
    """A button opening the sequence window."""
    if st.button(label, icon=':material/edit_note:', key=key, width='content',
                 help='The proteoform sequence shown for runs that were not identified'):
        st.session_state.pop(_LOADED_KEY, None)
        _dialog(file_manager)
