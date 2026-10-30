"""Server-side validation of layouts submitted by the Layout Manager editor."""

from src.render.layout_editor.catalog import (
    DECONV_COMPONENTS, DECONV_DEFAULT_LAYOUT, DECONV_SEQUENCE_COMPONENTS,
    TNT_COMPONENTS, TNT_DEFAULT_LAYOUT, validate_layout,
)


def test_default_layouts_are_valid():
    assert validate_layout([DECONV_DEFAULT_LAYOUT], DECONV_COMPONENTS) == ''
    assert validate_layout([DECONV_DEFAULT_LAYOUT + [['sequence_view']]],
                           DECONV_COMPONENTS + DECONV_SEQUENCE_COMPONENTS) == ''
    assert validate_layout([TNT_DEFAULT_LAYOUT], TNT_COMPONENTS) == ''


def test_names_are_unique_per_catalog():
    for catalog in (DECONV_COMPONENTS + DECONV_SEQUENCE_COMPONENTS, TNT_COMPONENTS):
        names = [c['name'] for c in catalog]
        assert len(names) == len(set(names))
        assert all(c.get('requires') in names for c in catalog if c.get('requires'))


def test_rejects_missing_requirement():
    assert 'needs Scan table' in validate_layout([[['mass_table']]], DECONV_COMPONENTS)


def test_rejects_sequence_components_without_sequence():
    assert 'Unknown component' in validate_layout(
        [[['scan_table', 'mass_table'], ['sequence_view']]], DECONV_COMPONENTS)


def test_rejects_shape_problems():
    assert validate_layout([], DECONV_COMPONENTS) == 'Empty layout'
    assert 'empty' in validate_layout([[['scan_table']], []], DECONV_COMPONENTS)
    assert 'at most 3' in validate_layout(
        [[['ms1_raw_heatmap', 'ms2_raw_heatmap', 'scan_table', 'fdr_plot']]], DECONV_COMPONENTS)
    assert 'twice' in validate_layout([[['scan_table'], ['scan_table']]], DECONV_COMPONENTS)
    assert 'At most 5' in validate_layout([[['scan_table']]] * 6, DECONV_COMPONENTS)


def test_same_component_allowed_in_different_experiments():
    assert validate_layout([[['scan_table']], [['scan_table']]], DECONV_COMPONENTS) == ''
