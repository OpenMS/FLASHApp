"""Links and 📌 of the Viewer leave out the fragments the sequence views publish."""

from pathlib import Path


def test_the_run_kind_names_the_matched_fragments():
    source = Path('content/TopDown/TopDownViewer.py').read_text()
    assert 'not_kept=(FRAGMENT_MATCHES_IDENTIFIER,)' in source


def test_a_link_leaves_out_the_matches():
    from src.layout_manager.layout import capture_defaults, link_params, read_link
    from src.layout_manager.ui import _stored, _without
    from src.render.render import FRAGMENT_MATCHES_IDENTIFIER

    class State:
        values = {'pscan': 4, FRAGMENT_MATCHES_IDENTIFIER: {'fragments': [{'name': f'b{i}'} for i in range(2000)]}}

        def get_all_selections(self):
            return dict(self.values)

    state = _without(lambda i: State(), (FRAGMENT_MATCHES_IDENTIFIER,))(0)
    entry = capture_defaults(_stored(state), 'run_1')
    assert entry == {'dataset': 'run_1', 'selection': {'pscan': 4}}
    assert read_link(link_params('native', [entry]))['experiments'][0]['selection'] == {'pscan': 4}
    # The components still see them
    assert FRAGMENT_MATCHES_IDENTIFIER in state.get_all_selections()
