"""A list of runs to put in order by dragging, with Randomize and inserted runs.

The rows of the Acquisition table are Streamlit widgets and cannot be dragged,
so the order is set in this list (``frontend/index.html``, plain JavaScript
with no build step), drawn in a dialog. Dragging, Randomize and inserting
washes, QC standards or blanks never rerun the script; Python hears from it
only when the user presses Save.
"""

from pathlib import Path

import streamlit as st
import streamlit.components.v1 as st_components

_list = st_components.declare_component('flash_run_order', path=str(Path(__file__).parent / 'frontend'))


def run_order(items: list, key: str, inserts: list = ()):
    """Draw the list; in the run Save was pressed, the new runs in order, else None.

    ``items``: one dict per run, ``{'name', 'method', 'kind'}``; ``kind`` is
    None for a sample run, else the kind of an inserted run (a wash, …), which
    Randomize leaves in place and Insert may reuse or drop.
    ``inserts``: the kinds that can be inserted, ``{'kind', 'label', 'method',
    'base', 'every'}``; new runs of a kind are named ``<base>_<n>``, n the first
    free number, and get its method. ``every`` is the default spacing.

    The result has, per position, the index of a run of ``items`` or a new
    run ``{'name', 'method'}``. Runs may only be dropped when they are inserted
    runs; anything else (a list from before the table changed) gives None.
    """
    handled_key = f'{key}_handled'
    value = _list(items=[{'name': i.get('name') or '', 'method': i.get('method') or '',
                          'kind': i.get('kind')} for i in items],
                  inserts=list(inserts), key=key, default=None)
    if not value or value.get('nonce') == st.session_state.get(handled_key):
        return None
    st.session_state[handled_key] = value.get('nonce')
    order = value.get('order')
    if value.get('action') != 'save' or not isinstance(order, list):
        return None
    kept = [e for e in order if isinstance(e, int)]
    new = [e for e in order if not isinstance(e, int)]
    methods = {i['method'] for i in inserts}
    if (len(set(kept)) != len(kept) or any(not 0 <= i < len(items) for i in kept)
            or any(not items[i].get('kind') for i in set(range(len(items))) - set(kept))
            or any(not isinstance(e, dict) or not e.get('name') or e.get('method') not in methods
                   for e in new)):
        return None
    return [e if isinstance(e, int) else {'name': str(e['name']), 'method': e['method']} for e in order]
