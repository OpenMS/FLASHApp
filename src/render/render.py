import threading

import streamlit as st
from streamlit.runtime.scriptrunner import get_script_run_ctx

from src.render.util import payload_hash
from src.render.StateTracker import StateTracker
from src.render.initialize import initialize_data
from src.render.update import update_data, filter_data
from src.render.components import get_component_function

# @st.fragment()
def render_component(
    components, data, component_key='flash_viewer_grid', on_change=None, 
    additional_data=None, tool=None, state_tracker=None
):
    # Map arguments
    out_components = []
    for row in components:
        out_components.append(list(map(
            lambda component: {
                "componentArgs": component.componentArgs.__dict__
            }, 
            row
        )))
    
    # Get State
    state = state_tracker.getState()

    # Cleared selections now arrive (and are stored) as `None` rather than being
    # dropped, so the frontend can round-trip a deselect. update/filter logic uses
    # the "key not in selection_store" convention, so drop None-valued keys for the
    # data computation while still echoing the full state (incl. nulls) back so the
    # frontend can clear those fields in every component.
    active_state = {k: v for k, v in state.items() if v is not None}

    # Update data with current session state
    data = update_data(data, out_components, active_state, additional_data, tool)

    # Filter data based on selection
    data = filter_data(
        data, out_components, active_state, additional_data, tool
    )

    # Hash the filtered data together with the tracker id: the frontend skips renders
    # whose hash is unchanged, and a new experiment must rebuild every cell even when
    # its data is byte-identical to the previous experiment's (re-run of the same file).
    data['hash'] = payload_hash(data, state['id'])

    # Render component
    data['selection_store'] = state
    new_state = get_component_function()(
        components=out_components,
        key=component_key,
        **data
    )

    # Update state
    if new_state is not None:
        updated = state_tracker.updateState(new_state)

        if updated:
            st.rerun(scope='app')


# One lock per Streamlit session. With runner.fastReruns (the default) Streamlit
# starts a new script thread while the previous one is still executing; two threads
# rendering the grid at the same time share the session's polars/pyarrow objects and
# crashed the server (see docs/white-screen-root-cause.md). Serialising render_grid
# per session removes that overlap: a runner that has been stopped releases the lock
# as soon as its next Streamlit call raises. Keyed by session id (not stored in
# st.session_state) so two concurrent first calls cannot create two different locks.
_GRID_LOCKS: dict[str, threading.Lock] = {}
_GRID_LOCKS_GUARD = threading.Lock()


def _grid_lock() -> threading.Lock:
    ctx = get_script_run_ctx()
    key = ctx.session_id if ctx is not None else '__no_session__'
    with _GRID_LOCKS_GUARD:
        return _GRID_LOCKS.setdefault(key, threading.Lock())


def render_grid(
    selected_data, layout_info_per_exp, file_manager, tool, identifier,
    grid_key='flash_viewer_grid'
):
    with _grid_lock():
        _render_grid_unlocked(
            selected_data, layout_info_per_exp, file_manager, tool, identifier,
            grid_key
        )


def _render_grid_unlocked(
    selected_data, layout_info_per_exp, file_manager, tool, identifier,
    grid_key='flash_viewer_grid'
):
    default_data = {'dataset' : selected_data}
    default_state = StateTracker()
    
    # Set up session state
    for name, default in zip(
        ['plot_data', 'state_tracker'], [default_data, default_state]
    ):
        if name not in st.session_state:
            st.session_state[name] = {}
        if tool not in st.session_state[name]:
            st.session_state[name][tool] = {}
        if identifier not in st.session_state[name][tool]:
            st.session_state[name][tool][identifier] = default

    # Check if dataset has changed
    if st.session_state['plot_data'][tool][identifier]['dataset'] != selected_data:
        st.session_state['plot_data'][tool][identifier] = default_data
        st.session_state['state_tracker'][tool][identifier] = default

    for row_index, row in enumerate(layout_info_per_exp):
        columns = st.columns(len(row))
        for col, (col_index, comp_name) in zip(columns, enumerate(row)):

            
            # Inititalize component data
            if comp_name not in st.session_state.plot_data[tool][identifier]:
                st.session_state.plot_data[tool][identifier][comp_name] = initialize_data(
                    comp_name, selected_data, file_manager, tool
                )

            # Get State
            state_tracker = st.session_state.state_tracker[tool][identifier]

            # Get data
            data_to_send, components, additional_data = (
                st.session_state.plot_data[tool][identifier][comp_name]
            )

            # Create component
            with col:
                render_component(
                    components=components, 
                    data=data_to_send, 
                    component_key=f"{grid_key}_{row_index}_{col_index}",
                    additional_data=additional_data,
                    tool=tool,
                    state_tracker=state_tracker
                )
