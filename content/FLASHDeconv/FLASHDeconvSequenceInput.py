import streamlit as st
import re
from src.common.common import page_setup, save_params, v_space
from src.workflow.FileManager import FileManager
from pathlib import Path


# Setup cache access
file_manager = FileManager(
    st.session_state["workspace"],
    Path(st.session_state['workspace'], 'flashdeconv', 'cache')
)

def set_sequence(input_sequence, fixed_mod_cysteine=None, fixed_mod_methionine=None):
    file_manager.store_data('sequence', 'sequence', 
        {
            'input_sequence' : input_sequence, 
            'fixed_mod_cysteine' : fixed_mod_cysteine, 
            'fixed_mod_methionine' : fixed_mod_methionine
        }
    )

def get_sequence():
    # Check if layout has been set
    if not file_manager.result_exists('sequence', 'sequence'):
        return None
    # fetch layout from cache
    sequence = file_manager.get_results('sequence', 'sequence')['sequence']

    return sequence['input_sequence'], sequence['fixed_mod_cysteine'], sequence['fixed_mod_methionine'] 

def emptySequenceInput():
    if file_manager.result_exists('sequence', 'sequence'):
        file_manager.remove_results('sequence')
    st.session_state['reset_sequence_input'] = True

fixed_mod_cysteine = ['No modification',
                      'Carbamidomethyl (+57)',
                      'Carboxymethyl (+58)',
                      'Xlink:Disulfide (-1 per C)']
                      # 'S-carboxamidoethyl-L-cysteine',
                      # 'S-carboxamidoethly-L-cysteine',
                      # 'S-pyridylethyl-L-cysteine',
                      # 'S-carboxamidomethly-L-cysteine',
                      # 'cyteine mercaptoethanol']
fixed_mod_methionine = ['No modification',
                        'L-methionine sulfoxide (+16)',
                        'L-methionine sulfone (+32)']


def validateSequenceInput(input_seq):
    # remove all white spaces
    seq = ''.join(input_seq.split())
    if not seq: return False

    pattern = re.compile("^[ac-ik-wyAC-IK-WYXx]+$")  # only alphabet except for BJXZ
    if not pattern.match(seq):
        return False
    return True

# page initialization
params = page_setup()

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
    seq, cys_mod, met_mod = cached_data

    if 'sequence_text' not in st.session_state:
        st.session_state['sequence_text'] = seq
    if (
        (cys_mod is not None) 
        and ('selected_fixed_mod_cysteine' not in st.session_state)
    ):
        st.session_state['selected_fixed_mod_cysteine'] = cys_mod
    if (
        (met_mod is not None)
        and ('selected_fixed_mod_methionine' not in st.session_state)
    ):
        st.session_state['selected_fixed_mod_methionine'] = met_mod

# clean up the entries of form, if needed
if st.session_state['reset_sequence_input']:
    st.session_state['sequence_text'] = ''
    st.session_state['selected_fixed_mod_cysteine'] = 'No modification'
    st.session_state['selected_fixed_mod_methionine'] = 'No modification'
    st.session_state['reset_sequence_input'] = False

with st.form('sequence_input'):
    # sequence
    st.text_area('Proteoform sequence', key='sequence_text')

    # fixed modification
    c1, c2 = st.columns(2)
    c1.selectbox('Fixed modification: Cysteine', fixed_mod_cysteine,
                 key='selected_fixed_mod_cysteine', placeholder='No modification')
    c2.selectbox('Fixed modification: Methionine', fixed_mod_methionine,
                 key='selected_fixed_mod_methionine', placeholder='No modification')
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
            cys_mod = None
            met_mod = None
            if 'selected_fixed_mod_cysteine' in st.session_state \
                    and st.session_state['selected_fixed_mod_cysteine'] != 'No modification':
                cys_mod = st.session_state.selected_fixed_mod_cysteine
            if 'selected_fixed_mod_methionine' in st.session_state \
                    and st.session_state['selected_fixed_mod_methionine'] != 'No modification':
                met_mod = st.session_state.selected_fixed_mod_methionine

            set_sequence(seq, cys_mod, met_mod)

        else:
            st.error('Error: sequence input is not valid')

st.info("""
**💡 NOTE** 

- This is only needed when the "Sequence View" component will be used in 👀Viewer
        
- Variable modifications can be specified within the ”Sequence View” component in 👀Viewer.

- Only one protein sequence is allowed
""")

save_params(params)
