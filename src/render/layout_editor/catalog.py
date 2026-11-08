"""Viewer components offered by the Layout Manager pages, and layout validation.

Each entry's ``name`` is the identifier stored in the saved layout and read by
the Viewer pages, so it must not change. ``requires`` names another component
that has to be in the same experiment (e.g. the mass table follows the scan
table's selection), or lists several of which one has to be there
(``requires_label`` then names them briefly in the editor).
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
    {'name': 'survey_spectrum', 'label': 'Precursor spectrum (MS1)', 'group': 'Spectra & plots', 'requires': 'scan_table'},
    {'name': '3D_SN_plot', 'label': '3D S/N plot', 'group': 'Spectra & plots', 'requires': 'mass_table'},
    {'name': 'ms1_3d_map', 'label': '3D MS1 map (peaks by m/z, RT and intensity)', 'group': 'Heatmaps'},
    {'name': 'fdr_plot', 'label': 'Score distribution plot', 'group': 'Spectra & plots'},
    {'name': 'chromatogram', 'label': 'Chromatogram (TIC, base peak, deconvolved)', 'group': 'Spectra & plots'},
    {'name': 'isotope_envelope', 'label': 'Isotope envelope fit', 'group': 'Validation', 'requires': 'scan_table'},
    {'name': 'precursor_xic', 'label': 'Precursor XIC by charge', 'group': 'Validation', 'requires': 'scan_table'},
    {'name': 'run_qc', 'label': 'Run QC (masses, charges, QScores)', 'group': 'Quality control'},
]

# Offered only once a sequence has been set in the sequence view.
DECONV_SEQUENCE_COMPONENTS = [
    {'name': 'sequence_view', 'label': 'Sequence view', 'group': 'Sequence', 'requires': 'mass_table'},
    {'name': 'internal_fragment_map', 'label': 'Internal fragment map', 'group': 'Sequence', 'requires': 'mass_table'},
]

# Components following the selected proteoform: either table selects it
_SELECTED_PROTEOFORM = ['protein_table', 'proteoform_table']

TNT_COMPONENTS = [
    {'name': 'protein_table', 'label': 'Protein table', 'group': 'Tables'},
    {'name': 'tag_table', 'label': 'Tag table', 'group': 'Tables', 'requires': _SELECTED_PROTEOFORM, 'requires_label': 'a protein or proteoform table'},
    {'name': 'protein_group_table', 'label': 'Proteins', 'group': 'Tables'},
    {'name': 'proteoform_table', 'label': 'Proteoforms of the protein', 'group': 'Tables', 'requires': 'protein_group_table'},
    {'name': 'proteoform_prsm_table', 'label': 'PrSMs of the proteoform', 'group': 'Tables', 'requires': 'proteoform_table'},
    {'name': 'prsm_mass_table', 'label': 'PrSM mass table', 'group': 'Tables', 'requires': _SELECTED_PROTEOFORM, 'requires_label': 'a protein or proteoform table'},
    {'name': 'sequence_view', 'label': 'Sequence view', 'group': 'Sequence', 'requires': _SELECTED_PROTEOFORM, 'requires_label': 'a protein or proteoform table'},
    {'name': 'structure_view', 'label': 'Structure view (3D)', 'group': 'Sequence', 'requires': _SELECTED_PROTEOFORM, 'requires_label': 'a protein or proteoform table'},
    # Internal fragment map is disabled for FLASHTnT for now.
    {'name': 'combined_spectrum', 'label': 'Sequence tag view', 'group': 'Sequence', 'requires': 'tag_table'},
    {'name': 'survey_spectrum', 'label': 'Precursor spectrum (MS1)', 'group': 'Spectra & plots', 'requires': _SELECTED_PROTEOFORM, 'requires_label': 'a protein or proteoform table'},
    {'name': 'id_fdr_plot', 'label': 'Score distribution plot', 'group': 'Spectra & plots'},
    {'name': 'chromatogram', 'label': 'Chromatogram (TIC, base peak, deconvolved)', 'group': 'Spectra & plots'},
    {'name': 'fragment_error_plot', 'label': 'Fragment mass errors', 'group': 'Validation', 'requires': _SELECTED_PROTEOFORM, 'requires_label': 'a protein or proteoform table'},
    {'name': 'isotope_envelope', 'label': 'Isotope envelope fit', 'group': 'Validation', 'requires': _SELECTED_PROTEOFORM, 'requires_label': 'a protein or proteoform table'},
    {'name': 'precursor_xic', 'label': 'Precursor XIC by charge', 'group': 'Validation', 'requires': _SELECTED_PROTEOFORM, 'requires_label': 'a protein or proteoform table'},
    {'name': 'run_qc', 'label': 'Run QC (masses, charges, QScores)', 'group': 'Quality control'},
    {'name': 'id_qc', 'label': 'Identification QC (scores, IDs over time, coverage, fragment error)', 'group': 'Quality control'},
    {'name': 'proteoform_heatmap', 'label': 'Proteoform heatmap (PrSMs by RT and mass)', 'group': 'Heatmaps'},
    {'name': 'mass_shifts', 'label': 'Mass shifts (modifications of the proteoforms)', 'group': 'Spectra & plots'},
    {'name': 'go_enrichment', 'label': 'GO enrichment of the selected proteins', 'group': 'Functional analysis',
     'requires': ['protein_table', 'protein_group_table', 'proteoform_table', 'proteoform_prsm_table', 'proteoform_heatmap'],
     'requires_label': 'a protein, proteoform or PrSM table, or the proteoform heatmap'},
    {'name': 'ms1_raw_heatmap', 'label': 'MS1 raw heatmap', 'group': 'Heatmaps'},
    {'name': 'ms1_deconv_heat_map', 'label': 'MS1 deconvolved heatmap', 'group': 'Heatmaps'},
]

# The Viewer's one combined view: every component of both kinds. Selecting a
# proteoform or PrSM also selects its spectrum, so the spectra, the mass table
# and the validation plots follow either the scan table or the identifications.
_SELECTS_SCAN = ['scan_table', 'protein_table', 'proteoform_table', 'proteoform_prsm_table', 'proteoform_heatmap']
_SELECTS_SCAN_LABEL = 'the scan table or an identification table'


# Isobaric labels (TMT, iTRAQ): the reporter ions FLASHDeconv measured in a run
ISOBARIC_COMPONENTS = [
    {'name': 'reporter_table', 'label': 'Reporter ions per spectrum (TMT, iTRAQ)', 'group': 'Quantification'},
    {'name': 'reporter_ions', 'label': 'Reporter ions of the selected spectrum', 'group': 'Quantification'},
]


# What a run needs to have been processed with for a component (``needs``):
# the identification components have nothing to show for a run deconvolved
# only (the sequence view then shows the sequence set in the view), the
# reporter ion components nothing without isobaric labels. The Viewer locks
# the views that need what the shown run lacks and offers to add it.
IDENTIFICATION = 'identification'
ISOBARIC = 'isobaric'
IDENTIFICATION_ONLY = (
    {c['name'] for c in TNT_COMPONENTS} - {c['name'] for c in DECONV_COMPONENTS} - {'sequence_view'}
)
ISOBARIC_ONLY = {c['name'] for c in ISOBARIC_COMPONENTS}


def _combined():
    by_name = {}
    for component in DECONV_COMPONENTS + TNT_COMPONENTS:
        by_name.setdefault(component['name'], dict(component))
    for name in ('mass_table', 'tag_table', 'deconv_spectrum', 'anno_spectrum', 'survey_spectrum', 'isotope_envelope', 'precursor_xic'):
        by_name[name].update(requires=_SELECTS_SCAN, requires_label=_SELECTS_SCAN_LABEL)
    # The proteoform of an identified run, or the sequence set in the view
    by_name['sequence_view'].update(
        requires=_SELECTED_PROTEOFORM + ['mass_table'],
        requires_label='a protein or proteoform table, or the mass table',
    )
    by_name.update({c['name']: dict(c) for c in ISOBARIC_COMPONENTS})
    by_name['reporter_ions'].update(requires=_SELECTS_SCAN + ['reporter_table'], requires_label=_SELECTS_SCAN_LABEL + ' or the reporter ion table')
    for name in IDENTIFICATION_ONLY:
        by_name[name]['needs'] = IDENTIFICATION
    for name in ISOBARIC_ONLY:
        by_name[name]['needs'] = ISOBARIC
    by_name['fdr_plot']['label'] = 'Score distribution plot (masses)'
    by_name['id_fdr_plot']['label'] = 'Score distribution plot (identifications)'
    groups = [
        'Heatmaps', 'Tables', 'Spectra & plots', 'Sequence', 'Quantification', 'Functional analysis',
        'Validation', 'Quality control',
    ]
    return sorted(by_name.values(), key=lambda c: groups.index(c['group']))


TOPDOWN_COMPONENTS = _combined()
# Offered only once a sequence has been set in the sequence view.
TOPDOWN_SEQUENCE_COMPONENTS = [c for c in DECONV_SEQUENCE_COMPONENTS if c['name'] != 'sequence_view']

# Mirrors DEFAULT_LAYOUT in the FLASHDeconv Viewer, which uses it when nothing is saved.
DECONV_DEFAULT_LAYOUT = [['ms1_deconv_heat_map'], ['scan_table', 'mass_table'],
                         ['anno_spectrum', 'deconv_spectrum'], ['3D_SN_plot']]
# FLASHTnT defaults by sample type (src/sample_type.py): a native or targeted
# run has few spectra, each worth inspecting in depth; a lysate run has many
# identifications, browsed protein by protein.
# The raw and the deconvolved spectrum below the sequence view show the
# fragment ions it matched.
TNT_NATIVE_LAYOUT = [['protein_table'], ['sequence_view'], ['anno_spectrum', 'deconv_spectrum'],
                     ['structure_view'], ['tag_table'], ['combined_spectrum']]
TNT_LYSATE_LAYOUT = [
    ['protein_group_table', 'proteoform_table', 'proteoform_prsm_table'],
    ['proteoform_heatmap'],
    ['sequence_view'],
    ['anno_spectrum', 'deconv_spectrum'],
    ['fragment_error_plot', 'precursor_xic'],
    ['id_qc'],
    ['go_enrichment'],
]
TNT_DEFAULT_LAYOUT = TNT_NATIVE_LAYOUT
TNT_DEFAULT_LAYOUTS = {'Native': TNT_NATIVE_LAYOUT, 'Lysate': TNT_LYSATE_LAYOUT}

# A comparison of runs (Viewer entry "Comparison"): one grid, built from these
QUANT_COMPONENTS = [
    {'name': 'quant_summary', 'label': 'Comparison summary', 'group': 'Quantification'},
    {'name': 'quant_table', 'label': 'Comparison table', 'group': 'Quantification'},
    {'name': 'quant_volcano', 'label': 'Volcano plot', 'group': 'Quantification'},
    {'name': 'quant_profile', 'label': 'Abundance per run', 'group': 'Quantification', 'requires': ('quant_table', 'quant_volcano'), 'requires_label': 'the comparison table or the volcano plot'},
    {'name': 'quant_runs', 'label': 'Selected entry per run', 'group': 'Quantification', 'requires': ('quant_table', 'quant_volcano'), 'requires_label': 'the comparison table or the volcano plot'},
    {'name': 'quant_pca', 'label': 'Sample similarity (PCA)', 'group': 'Quantification'},
    {'name': 'quant_heatmap', 'label': 'Clustered heatmap', 'group': 'Quantification'},
    {'name': 'quant_go', 'label': 'GO enrichment of the selected entries', 'group': 'Functional analysis',
     'requires': ('quant_table', 'quant_volcano', 'quant_heatmap'),
     'requires_label': 'the comparison table, the volcano plot or the clustered heatmap'},
    {'name': 'quant_id_counts', 'label': 'Identifications per run', 'group': 'Identification'},
    {'name': 'quant_id_frequency', 'label': 'Found in how many runs', 'group': 'Identification'},
    {'name': 'quant_id_overlap', 'label': 'Shared identifications', 'group': 'Identification'},
    {'name': 'quant_id_table', 'label': 'Identifications across runs', 'group': 'Identification'},
    {'name': 'quant_id_prsms', 'label': 'PrSMs of the selected identification', 'group': 'Identification',
     'requires': 'quant_id_table', 'requires_label': 'the identifications across runs'},
    {'name': 'quant_prsm_sequence', 'label': 'Sequence view of the selected PrSM', 'group': 'Identification',
     'requires': 'quant_id_prsms', 'requires_label': 'the PrSMs of the selected identification'},
    {'name': 'quant_prsm_spectrum', 'label': 'Spectrum of the selected PrSM', 'group': 'Identification',
     'requires': 'quant_id_prsms', 'requires_label': 'the PrSMs of the selected identification'},
    {'name': 'quant_reporter_intensities', 'label': 'Reporter ion intensity per label', 'group': 'Reporter ions'},
    {'name': 'quant_reporter_ratios', 'label': 'Reporter ion ratios', 'group': 'Reporter ions'},
    {'name': 'quant_reporter_labels', 'label': 'Reporter ions found per spectrum', 'group': 'Reporter ions'},
    {'name': 'quant_reporter_spectra', 'label': 'Reporter ions of every spectrum', 'group': 'Reporter ions'},
    {'name': 'quant_reporter_ions', 'label': 'Reporter ions of the selected spectrum', 'group': 'Reporter ions',
     'requires': 'quant_reporter_spectra', 'requires_label': 'the reporter ions of every spectrum'},
    {'name': 'quant_feature_map', 'label': 'Consensus feature map', 'group': 'Consensus features'},
    {'name': 'quant_feature_intensity', 'label': 'Feature intensity per run', 'group': 'Consensus features'},
    {'name': 'quant_feature_table', 'label': 'Consensus features', 'group': 'Consensus features'},
    {'name': 'quant_feature_profile', 'label': 'Selected feature per run', 'group': 'Consensus features',
     'requires': ('quant_feature_table', 'quant_feature_map'),
     'requires_label': 'the consensus features or the feature map'},
    {'name': 'quant_totals', 'label': 'Total intensity per sample', 'group': 'Quality control'},
    {'name': 'quant_correlation', 'label': 'Sample correlation', 'group': 'Quality control'},
    {'name': 'quant_cv', 'label': 'Variation within conditions (CV)', 'group': 'Quality control'},
    {'name': 'quant_match_quality', 'label': 'Match between runs', 'group': 'Quality control'},
    {'name': 'quant_missing_values', 'label': 'Missing values', 'group': 'Quality control'},
    {'name': 'quant_imputed_values', 'label': 'Imputed values', 'group': 'Quality control'},
    {'name': 'quant_normalization', 'label': 'Normalization', 'group': 'Quality control'},
    {'name': 'quant_pvalues', 'label': 'p-value distribution', 'group': 'Quality control'},
]
# The four plots above, which the layouts saved before they were split name
# together (src/render/quant_view.py, split_diagnostics)
QUANT_DIAGNOSTICS = [
    ['quant_missing_values', 'quant_imputed_values'], ['quant_normalization', 'quant_pvalues'],
]
# Identification: what each run identified and what the runs share, then
# every identification with its status run by run, the PrSMs of the selected
# one in every run and the sequence view and spectrum of the selected PrSM
# (from its run). Needs no conditions.
QUANT_IDENTIFICATION_LAYOUT = [
    ['quant_id_counts', 'quant_id_frequency'], ['quant_id_overlap'],
    ['quant_id_table'], ['quant_id_prsms'], ['quant_prsm_sequence'], ['quant_prsm_spectrum'],
]
# Quantification: what changed between the conditions, each entry followed
# run by run, how the samples cluster and what the changed entries do
QUANT_QUANTIFICATION_LAYOUT = [
    ['quant_summary'], ['quant_volcano', 'quant_profile'], ['quant_table'], ['quant_runs', 'quant_pca'],
    ['quant_heatmap'], ['quant_go'],
]
QUANT_DEFAULT_LAYOUT = QUANT_QUANTIFICATION_LAYOUT
# Quality control of a label-free comparison: runs measured one by one, so
# the consensus features (where they are, how each run found them, each
# feature run by run), how much was matched between runs, whether the runs
# were loaded alike and how well replicates agree
QUANT_LABEL_FREE_QC_LAYOUT = [
    ['quant_feature_map', 'quant_feature_intensity'], ['quant_feature_table', 'quant_feature_profile'],
    ['quant_id_counts', 'quant_match_quality'], ['quant_totals', 'quant_normalization'],
    ['quant_correlation', 'quant_cv'], ['quant_missing_values', 'quant_imputed_values'], ['quant_pvalues'],
]
# Quality control of isobaric labels (TMT, iTRAQ): the channels of one pooled
# set, so the reporter ions (intensity and ratio per label, labels found per
# spectrum, every labelled spectrum and its ions), whether the channels were
# loaded alike, how they correlate and cluster, then the comparison steps
QUANT_ISOBARIC_QC_LAYOUT = [
    ['quant_reporter_intensities', 'quant_reporter_ratios'], ['quant_reporter_labels', 'quant_totals'],
    ['quant_reporter_spectra', 'quant_reporter_ions'], ['quant_correlation', 'quant_pca'],
    ['quant_cv', 'quant_normalization'], ['quant_missing_values', 'quant_imputed_values'], ['quant_pvalues'],
]

MAX_EXPERIMENTS = 5
MAX_COLUMNS = 3


def requirements(component):
    """Components one of which has to be in the same experiment."""
    required = component.get('requires')
    if not required:
        return []
    return [required] if isinstance(required, str) else list(required)


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
            required = requirements(by_name[name])
            if required and not any(r in names for r in required):
                return '%s needs %s in experiment %d' % (
                    by_name[name]['label'],
                    ' or '.join(by_name[r]['label'] for r in required), exp_index + 1)
    return ''
