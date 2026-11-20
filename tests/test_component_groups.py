"""Every Viewer component sits in one of the four component groups."""
from src.render.layout_editor import catalog
from src.workflow.guide import COMPONENT_GROUPS, FEATURES, IDENTIFICATION, MS_DATA, QUANTIFICATION


def test_every_component_is_in_one_of_the_four_groups():
    for components in (catalog.TOPDOWN_COMPONENTS, catalog.TOPDOWN_SEQUENCE_COMPONENTS, catalog.QUANT_COMPONENTS):
        assert {c['group'] for c in components} <= set(COMPONENT_GROUPS)


def test_the_palette_lists_the_groups_in_order():
    order = [COMPONENT_GROUPS.index(c['group']) for c in catalog.TOPDOWN_COMPONENTS]
    assert order == sorted(order)


def test_components_go_with_the_data_they_are_about():
    group = {c['name']: c['group'] for c in catalog.TOPDOWN_COMPONENTS + catalog.QUANT_COMPONENTS}
    assert group['scan_table'] == MS_DATA
    assert group['mass_table'] == FEATURES and group['run_qc'] == FEATURES
    assert group['id_qc'] == IDENTIFICATION and group['go_enrichment'] == IDENTIFICATION
    assert group['reporter_table'] == QUANTIFICATION and group['quant_cv'] == QUANTIFICATION
