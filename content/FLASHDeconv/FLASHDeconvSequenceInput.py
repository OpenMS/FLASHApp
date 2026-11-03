import streamlit as st
import re
from src.common.common import page_setup, save_params, v_space
from src.workflow.FileManager import FileManager
from src.render import modifications
from pathlib import Path

# page initialization
params = page_setup()

# Setup cache access
file_manager = FileManager(
    st.session_state["workspace"],
    Path(st.session_state['workspace'], 'cache')
)

# Sites offered for a fixed modification: the termini, then the residues
MOD_SITES = [modifications.N_TERM, modifications.C_TERM] + list('ACDEFGHIKLMNPQRSTUVWY')


def set_sequence(input_sequence, fixed_modifications):
    file_manager.store_data('sequence', 'sequence',
        {
            'input_sequence' : input_sequence,
            # Kept for readers of the earlier format; the library ids replace them
            'fixed_mod_cysteine' : None,
            'fixed_mod_methionine' : None,
            'fixed_modifications' : list(fixed_modifications),
        }
    )

def get_sequence():
    # Check if a sequence has been saved
    if not file_manager.result_exists('sequence', 'sequence'):
        return None
    sequence = file_manager.get_results('sequence', 'sequence')['sequence']
    mod_ids = sequence.get('fixed_modifications')
    if mod_ids is None:
        # Saved before the modification library: Cys/Met names
        mod_ids = modifications.legacy_ids(
            sequence.get('fixed_mod_cysteine'), sequence.get('fixed_mod_methionine')
        )
    return sequence['input_sequence'], list(mod_ids)

def emptySequenceInput():
    if file_manager.result_exists('sequence', 'sequence'):
        file_manager.remove_results('sequence')
    st.session_state['reset_sequence_input'] = True


def validateSequenceInput(input_seq):
    # remove all white spaces
    seq = ''.join(input_seq.split())
    if not seq: return False

    pattern = re.compile("^[ac-ik-wyAC-IK-WYXx]+$")  # only alphabet except for BJXZ
    if not pattern.match(seq):
        return False
    return True


def save_modifications():
    """Store the fixed modifications with the saved sequence, if there is one."""
    cached = get_sequence()
    if cached is not None:
        set_sequence(cached[0], st.session_state['fixed_modifications'])


def mod_label(entry):
    label = f"{entry['id']}   {entry['mass']:+.4f} Da"
    if entry['unimod']:
        label += f"   {entry['unimod']}"
    if entry['description'] and entry['description'] != entry['name']:
        label += f"   {entry['description']}"
    return label


# for resetting the form (cannot be done after form is instantiated)
if 'reset_sequence_input' not in st.session_state:
    st.session_state['reset_sequence_input'] = False

# title and reset buttons
c1, c2 = st.columns([8, 1])
c1.title("Proteoform Sequence Input")
v_space(1, c2)
if c2.button('Reset'):
    emptySequenceInput()

cached_data = get_sequence()
if cached_data is not None:
    seq, mod_ids = cached_data
    if 'sequence_text' not in st.session_state:
        st.session_state['sequence_text'] = seq
    if 'fixed_modifications' not in st.session_state:
        st.session_state['fixed_modifications'] = mod_ids
if 'fixed_modifications' not in st.session_state:
    st.session_state['fixed_modifications'] = []

# clean up the entries of form, if needed
if st.session_state['reset_sequence_input']:
    st.session_state['sequence_text'] = ''
    st.session_state['fixed_modifications'] = []
    st.session_state['reset_sequence_input'] = False

with st.form('sequence_input'):
    # sequence
    st.text_area('Proteoform sequence', key='sequence_text')
    _, c2 = st.columns([8, 1])
    submitted = c2.form_submit_button("Save")
    if submitted:
        if st.session_state['sequence_text'] == '':
            emptySequenceInput()
            st.rerun()
        elif validateSequenceInput(st.session_state['sequence_text']):

            st.success('Proteoform sequence is submitted: ' + st.session_state['sequence_text'])

            # save information for sequence view
            seq = ''.join(st.session_state['sequence_text'].split()).upper()
            set_sequence(seq, st.session_state['fixed_modifications'])

        else:
            st.error('Error: sequence input is not valid')

# Fixed modifications from the Unimod library (pyOpenMS ModificationsDB)
st.subheader('Fixed modifications')
st.caption(
    'Each modification goes on every residue of its kind, or on the terminus of '
    'the proteoform. Search by name, Unimod accession or mass shift.'
)
c1, c2, c3 = st.columns([1, 5, 1], vertical_alignment='bottom')
site = c1.selectbox('Site', MOD_SITES, index=MOD_SITES.index('C'), key='fixed_mod_site')
candidates = modifications.search(site=site)
choice = c2.selectbox(
    f'Modification ({len(candidates)} at {site})',
    candidates,
    index=None,
    format_func=mod_label,
    placeholder='Type to search, e.g. Carbamidomethyl, UniMod:4 or 57.02',
    key=f'fixed_mod_choice_{site}',
)
if c3.button('Add', disabled=choice is None, use_container_width=True):
    if choice['id'] not in st.session_state['fixed_modifications']:
        st.session_state['fixed_modifications'] = st.session_state['fixed_modifications'] + [choice['id']]
        save_modifications()
    st.rerun()

current_sequence = ''.join(st.session_state.get('sequence_text', '').split()).upper()
mod_ids = st.session_state['fixed_modifications']
if mod_ids:
    header = st.columns([1, 4, 2, 2, 1])
    for col, text in zip(header, ['**Site**', '**Modification**', '**Mass shift**', '**Sites in sequence**', '']):
        col.markdown(text)
    for mod_id in mod_ids:
        entry = modifications.by_id(mod_id)
        if entry is None:
            continue
        row = st.columns([1, 4, 2, 2, 1], vertical_alignment='center')
        row[0].write(entry['site'])
        row[1].write(entry['id'] if entry['description'] == entry['name'] else f"{entry['id']} · {entry['description']}")
        row[2].write(f"{entry['mass']:+.4f} Da")
        row[3].write(str(len(modifications.sites_in(current_sequence, entry))) if current_sequence else '-')
        if row[4].button('Remove', key=f'remove_fixed_mod_{mod_id}'):
            st.session_state['fixed_modifications'] = [m for m in mod_ids if m != mod_id]
            save_modifications()
            st.rerun()
    if current_sequence and validateSequenceInput(current_sequence):
        mass = modifications.modified_mass(current_sequence.replace('X', ''), mod_ids)
        st.markdown(f'Theoretical mass with fixed modifications: **{mass:.4f} Da**')
else:
    st.write('No fixed modifications.')

st.info("""
**💡 NOTE**

- This is only needed when the "Sequence View" component will be used in 👀Viewer

- Fixed modifications are saved with the sequence. Variable modifications can be specified within the "Sequence View" component in 👀Viewer.

- Only one protein sequence is allowed
""")

save_params(params)
