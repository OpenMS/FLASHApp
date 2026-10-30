"""Structure files uploaded or removed in the Structure View (FLASHTnT Viewer)."""
import streamlit as st

from src.render.structure import StructureStore, apply_structure_action, structure_key


def upload_errors():
    """{structure key: message} of the last rejected upload per protein."""
    return st.session_state.setdefault('structure_upload_errors', {})


def handle_structure_action(action, additional_data) -> bool:
    """Apply an action sent by the Structure View once. True if the stored
    structures changed or an error is to be shown, so the page reruns."""
    if not isinstance(action, dict) or additional_data is None:
        return False
    if 'structure_dir' not in additional_data:
        return False
    # Streamlit returns the component's last value on every rerun
    handled = st.session_state.setdefault('structure_actions_handled', set())
    nonce = action.get('nonce')
    if nonce is None or nonce in handled:
        return False
    handled.add(nonce)
    changed, error = apply_structure_action(
        StructureStore(additional_data['structure_dir']), action
    )
    key = structure_key(action.get('accession'))
    errors = upload_errors()
    if error is not None:
        errors[key] = error
        return True
    errors.pop(key, None)
    return changed
