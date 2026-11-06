"""Previous / Next row buttons on the Viewer's tables."""
import inspect

import pytest

openms_insight = pytest.importorskip("openms_insight")

from src.render import render  # noqa: E402


def test_row_step_follows_the_installed_insight():
    has_option = 'row_step' in inspect.signature(openms_insight.Table).parameters
    assert render.ROW_STEP == ({'row_step': True} if has_option else {})


@pytest.mark.parametrize('cell', ['scan_table', 'mass_table', 'protein_table', 'tag_table'])
def test_stepped_tables_pass_the_option(cell):
    source = inspect.getsource(render._build)
    block = source.split(f"if comp_name == '{cell}':", 1)[1].split('if comp_name ==', 1)[0]
    assert '**ROW_STEP' in block
