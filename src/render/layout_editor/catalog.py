"""Viewer components offered by the Layout Manager pages, and layout validation.

Each entry's ``name`` is the identifier stored in the saved layout and read by
the Viewer pages, so it must not change. ``requires`` names another component
that has to be in the same experiment (e.g. the mass table follows the scan
table's selection).
"""

DECONV_COMPONENTS = [
    {'name': 'ms1_raw_heatmap', 'label': 'MS1 raw heatmap', 'group': 'Heatmaps'},
    {'name': 'ms2_raw_heatmap', 'label': 'MS2 raw heatmap', 'group': 'Heatmaps'},
    {'name': 'ms1_deconv_heat_map', 'label': 'MS1 deconvolved heatmap', 'group': 'Heatmaps'},
    {'name': 'ms2_deconv_heat_map', 'label': 'MS2 deconvolved heatmap', 'group': 'Heatmaps'},
    {'name': 'scan_table', 'label': 'Scan table', 'group': 'Tables'},
    {'name': 'mass_table', 'label': 'Mass table', 'group': 'Tables', 'requires': 'scan_table'},
    {'name': 'deconv_spectrum', 'label': 'Deconvolved spectrum', 'group': 'Spectra & plots', 'requires': 'scan_table'},
    {'name': 'anno_spectrum', 'label': 'Raw spectrum', 'group': 'Spectra & plots', 'requires': 'scan_table'},
    {'name': '3D_SN_plot', 'label': '3D S/N plot', 'group': 'Spectra & plots', 'requires': 'mass_table'},
    {'name': 'fdr_plot', 'label': 'Score distribution plot', 'group': 'Spectra & plots'},
    {'name': 'isotope_envelope', 'label': 'Isotope envelope fit', 'group': 'Validation', 'requires': 'scan_table'},
    {'name': 'precursor_xic', 'label': 'Precursor XIC by charge', 'group': 'Validation', 'requires': 'scan_table'},
]

# Offered only once a sequence has been submitted on the Sequence Input page.
DECONV_SEQUENCE_COMPONENTS = [
    {'name': 'sequence_view', 'label': 'Sequence view', 'group': 'Sequence', 'requires': 'mass_table'},
    {'name': 'internal_fragment_map', 'label': 'Internal fragment map', 'group': 'Sequence', 'requires': 'mass_table'},
]

TNT_COMPONENTS = [
    {'name': 'protein_table', 'label': 'Protein table', 'group': 'Tables'},
    {'name': 'tag_table', 'label': 'Tag table', 'group': 'Tables', 'requires': 'protein_table'},
    {'name': 'sequence_view', 'label': 'Sequence view', 'group': 'Sequence', 'requires': 'protein_table'},
    {'name': 'structure_view', 'label': 'Structure view (3D)', 'group': 'Sequence', 'requires': 'protein_table'},
    # Internal fragment map is disabled for FLASHTnT for now.
    {'name': 'combined_spectrum', 'label': 'Sequence tag view', 'group': 'Sequence', 'requires': 'tag_table'},
    {'name': 'id_fdr_plot', 'label': 'Score distribution plot', 'group': 'Spectra & plots'},
    {'name': 'fragment_error_plot', 'label': 'Fragment mass errors', 'group': 'Validation', 'requires': 'protein_table'},
    {'name': 'isotope_envelope', 'label': 'Isotope envelope fit', 'group': 'Validation', 'requires': 'protein_table'},
    {'name': 'precursor_xic', 'label': 'Precursor XIC by charge', 'group': 'Validation', 'requires': 'protein_table'},
    {'name': 'ms1_raw_heatmap', 'label': 'MS1 raw heatmap', 'group': 'Heatmaps'},
    {'name': 'ms1_deconv_heat_map', 'label': 'MS1 deconvolved heatmap', 'group': 'Heatmaps'},
]

# Mirrors DEFAULT_LAYOUT in the Viewer pages, which use it when nothing is saved.
DECONV_DEFAULT_LAYOUT = [['ms1_deconv_heat_map'], ['scan_table', 'mass_table'],
                         ['anno_spectrum', 'deconv_spectrum'], ['3D_SN_plot'],
                         ['isotope_envelope', 'precursor_xic']]
TNT_DEFAULT_LAYOUT = [['protein_table'], ['sequence_view'], ['fragment_error_plot', 'isotope_envelope'],
                      ['structure_view'], ['tag_table'], ['combined_spectrum']]

MAX_EXPERIMENTS = 5
MAX_COLUMNS = 3


def validate_layout(layout, components):
    """Return an error message for an invalid layout, or '' if it can be saved.

    ``layout`` is the saved format: experiments -> rows -> component names.
    """
    by_name = {c['name']: c for c in components}
    if not isinstance(layout, list) or not layout:
        return 'Empty layout'
    if len(layout) > MAX_EXPERIMENTS:
        return 'At most %d experiments are supported' % MAX_EXPERIMENTS
    for exp_index, exp in enumerate(layout):
        if not isinstance(exp, list) or not any(exp):
            return 'Experiment %d is empty' % (exp_index + 1)
        names = []
        for row in exp:
            if not isinstance(row, list) or not row:
                return 'Experiment %d contains an empty row' % (exp_index + 1)
            if len(row) > MAX_COLUMNS:
                return 'A row holds at most %d components' % MAX_COLUMNS
            names.extend(row)
        for name in names:
            if name not in by_name:
                return 'Unknown component "%s"' % name
        if len(set(names)) != len(names):
            return 'Experiment %d contains a component twice' % (exp_index + 1)
        for name in names:
            required = by_name[name].get('requires')
            if required and required not in names:
                return '%s needs %s in experiment %d' % (
                    by_name[name]['label'], by_name[required]['label'], exp_index + 1)
    return ''
