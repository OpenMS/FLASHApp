"""A view's defaults and links leave out the selections components publish."""

from src.render.view_state import ViewState


class FakeState:
    def __init__(self, values):
        self.values = dict(values)

    def get_all_selections(self):
        return dict(self.values)

    def get_selection(self, identifier):
        return self.values.get(identifier)

    def set_selection(self, identifier, value):
        self.values[identifier] = value
        return True


def test_derived_selections_are_left_out():
    state = FakeState({'pscan': 4, 'protein': 1, 'fragment_matches': {'fragments': [1] * 500}})
    view = ViewState(state)
    assert view.stored_selections() == {'pscan': 4, 'protein': 1}
    # The rest goes to the state itself
    assert view.get_all_selections()['fragment_matches']
    view.set_selection('pscan', 7)
    assert state.values['pscan'] == 7 and view.get_selection('pscan') == 7


def test_stored_selections_of_a_state_that_keeps_some_apart():
    class Private(FakeState):
        def stored_selections(self):
            return {**self.values, 'row@table': 3}

    view = ViewState(Private({'fragment_matches': [], 'scan': 2}))
    assert view.stored_selections() == {'scan': 2, 'row@table': 3}


def test_capture_leaves_out_the_matches():
    from src.layout_manager.layout import capture_defaults, link_params, read_link
    state = FakeState({'pscan': 4, 'fragment_matches': {'fragments': [{'name': f'b{i}'} for i in range(2000)]}})
    entry = capture_defaults(ViewState(state).stored_selections(), 'run_1')
    assert entry == {'dataset': 'run_1', 'selection': {'pscan': 4}}
    assert read_link(link_params('native', [entry]))['experiments'][0]['selection'] == {'pscan': 4}
