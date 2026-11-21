"""What a view keeps of a grid's selections: for 📌 and for a link to it.

The view bar stores an experiment's selection and zoom (``capture_defaults``)
from what its ``StateManager`` holds. Some identifiers there are not picked
by the user but published by a component from what it shows, such as the
fragment ions a sequence view matched (``fragment_matches``): the component
writes them again once it is drawn. Keeping them would only make a link long
(a large proteoform's matches alone pass the length a link may have, and the
link then loses the selection). Streamlit-free.
"""

# Published by the components, never picked: left out of a view's defaults and links
DERIVED = ('fragment_matches',)


class ViewState:
    """A ``StateManager`` whose stored selections leave out the derived ones.

    Everything else (getting, setting and clearing a selection) goes to the
    state itself, so applying a view's defaults or a link works as before."""

    def __init__(self, state, derived=DERIVED):
        self._state = state
        self._derived = set(derived)

    def stored_selections(self):
        stored = getattr(self._state, 'stored_selections', None)
        values = stored() if callable(stored) else self._state.get_all_selections()
        return {k: v for k, v in (values or {}).items() if k not in self._derived}

    def __getattr__(self, name):
        return getattr(self._state, name)
