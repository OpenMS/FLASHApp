import streamlit as st

from src.render.util import hash_complex
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

    # Hash updated. filtered data
    data['hash'] = hash_complex(data)

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


def render_grid(
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
