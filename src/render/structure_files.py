"""Upload section for the Structure View's structure files (FLASHTnT Viewer)."""
import pandas as pd
import streamlit as st

from src.render.structure import (
    StructureStore, structure_dir, structure_format, structure_key, MAX_BYTES
)


def _selected_protein(identifier):
    tracker = st.session_state.get('state_tracker', {}).get('flashtnt', {}).get(identifier)
    if tracker is None:
        return None
    return tracker.currentState.get('proteinIndex')


def structure_files_section(file_manager, dataset_id, layout, identifier):
    """Expander below the grid to upload a PDB/mmCIF file for a protein, which
    the Structure View then shows instead of the AlphaFold DB model. Needed
    offline and for proteins without a UniProt accession. Shown only when the
    experiment's layout has a Structure View."""
    if not any('structure_view' in row for row in layout):
        return
    store = StructureStore(structure_dir(file_manager))
    proteins = file_manager.get_results(dataset_id, ['protein_dfs'])['protein_dfs']
    proteins = proteins[~proteins['accession'].astype(str).str.startswith('DECOY_')]

    options = {}
    for row in proteins.itertuples():
        key = structure_key(row.accession)
        if key is not None and key not in options:
            options[key] = f'{row.accession} · {row.description}'
    if not options:
        return

    selected = _selected_protein(identifier)
    default_key = None
    if selected is not None:
        match = proteins[proteins['index'] == selected]
        if len(match):
            default_key = structure_key(match['accession'].iloc[0])
    keys = list(options)
    stored = store.list()

    with st.expander('🧬 Structure files', expanded=False):
        st.caption(
            'The Structure View shows the AlphaFold DB model of the selected protein '
            'and keeps it in this workspace. Upload a PDB or mmCIF file to show another '
            'structure, or when this computer is offline.'
        )
        counter_key = f'structure_upload_counter_{identifier}'
        counter = st.session_state.setdefault(counter_key, 0)
        # Follow the protein selected in the viewer, until the user picks another
        select_key = f'structure_upload_protein_{identifier}'
        followed_key = f'structure_upload_followed_{identifier}'
        if default_key in keys and st.session_state.get(followed_key) != default_key:
            st.session_state[followed_key] = default_key
            st.session_state[select_key] = default_key
        if st.session_state.get(select_key) not in keys:
            st.session_state.pop(select_key, None)
        c1, c2 = st.columns(2)
        with c1:
            key = st.selectbox(
                'Protein', keys, format_func=lambda k: options[k], key=select_key,
            )
        with c2:
            uploaded = st.file_uploader(
                'PDB or mmCIF file', type=['pdb', 'ent', 'cif', 'mmcif'],
                key=f'structure_upload_{identifier}_{counter}',
            )
        if uploaded is not None:
            data = uploaded.getvalue()
            fmt = structure_format(uploaded.name, data) if len(data) <= MAX_BYTES else None
            if fmt is None:
                st.error(f'{uploaded.name} is not a PDB or mmCIF file (or larger than 50 MB).')
            else:
                store.save(key, data, fmt, 'upload', uploaded.name)
                # A fresh uploader, so the file is not stored again on the next rerun
                st.session_state[counter_key] = counter + 1
                st.rerun()

        if stored:
            st.dataframe(pd.DataFrame([
                {
                    'Protein': k,
                    'Structure': meta.get('name'),
                    'Source': 'AlphaFold DB' if meta.get('source') == 'alphafold' else 'Uploaded',
                } for k, meta in stored.items()
            ]), hide_index=True, use_container_width=True)
            if key in stored and st.button(
                f'Remove the structure of {key}', key=f'structure_remove_{identifier}'
            ):
                store.remove(key)
                st.rerun()
