"""Structure files uploaded or removed in the Structure View (FLASHTnT Viewer)."""
import streamlit as st

from src.render.structure import StructureStore, apply_structure_action, structure_key


def upload_errors():
    """{structure key: message} of the last rejected upload per protein."""
    return st.session_state.setdefault('structure_upload_errors', {})


def structure_action_handler(directory):
    """Handler for the Structure View's uploads and removals (applied once per
    action by openms-insight). True if the stored structures changed or an
    error is to be shown, so the page reruns."""
    def handle(action) -> bool:
        changed, error = apply_structure_action(StructureStore(directory), action)
        key = structure_key(action.get('accession'))
        errors = upload_errors()
        if error is not None:
            errors[key] = error
            return True
        errors.pop(key, None)
        return changed
    return handle
