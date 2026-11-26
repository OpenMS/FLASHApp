"""What each component of the Viewer shows and how to use it.

Every component shows an info icon at its top left whose tooltip is its text
here. Names are the catalog's (``layout_editor/catalog.py``).

``help()`` reads the component being drawn from a context variable that
``drawing(name)`` sets, so the builders need not pass their own name down; it
gives the text once per component, to the first plot or table drawn. A
context variable belongs to the thread of one session's script run, so two
sessions never see each other's component.
"""

from contextlib import contextmanager
from contextvars import ContextVar


_current: ContextVar = ContextVar('flash_viewer_component', default=None)

COMPONENTS = {
    # Runs: deconvolution
    'ms1_raw_heatmap': 'Every MS1 peak by retention time (x) and m/z (y), coloured by intensity. Zoom by dragging a box; a click picks the scan nearest to it.',
    'ms2_raw_heatmap': 'Every MS2 peak by retention time (x) and m/z (y), coloured by intensity. Zoom by dragging a box; a click picks the scan nearest to it.',
    'ms1_deconv_heat_map': 'The monoisotopic masses deconvolved from the MS1 spectra, by retention time (x) and mass (y). Proteoforms show as horizontal traces; a click picks the scan.',
    'ms2_deconv_heat_map': 'The fragment masses deconvolved from the MS2 spectra, by retention time (x) and mass (y). A click picks the scan.',
    'ms1_3d_map': 'The MS1 peaks in three dimensions: m/z, retention time and intensity. Drag to rotate, scroll to zoom.',
    'scan_table': 'Every spectrum of the run: scan number, MS level, retention time, precursor and the number of masses found. Picking a scan shows it in the spectra and the mass table.',
    'mass_table': 'The masses deconvolved from the picked scan, with charge range, intensity and quality score (QScore). Picking a mass marks its peaks in the raw spectrum.',
    'deconv_spectrum': 'The deconvolved spectrum of the picked scan: one line per neutral mass. The picked mass is highlighted; fragments matched by a sequence view are marked.',
    'anno_spectrum': 'The raw spectrum of the picked scan. The peaks of the picked mass are coloured by charge; matched fragments are marked and labelled.',
    'survey_spectrum': 'The MS1 spectrum the picked MS2 scan was isolated from, with its isolation window, so the precursor and its neighbours can be judged.',
    '3D_SN_plot': 'The peaks of the picked mass by charge (y), m/z (x) and intensity, signal against noise. A real mass shows a regular charge ladder above the noise.',
    'fdr_plot': 'How many target and decoy masses reach each QScore. Decoys show what wrong masses score; the further the targets lie to the right of them, the better.',
    'chromatogram': 'Total ion current, base peak and deconvolved intensity over retention time. Zoom by dragging; it shows when the sample eluted and whether spray was stable.',
    'isotope_envelope': 'The isotope peaks of the picked mass against the theoretical envelope of a protein of that mass. A good fit means the monoisotopic mass was assigned correctly.',
    'precursor_xic': 'The intensity of the picked precursor over retention time, one trace per charge. Traces that rise and fall together belong to one proteoform.',
    'run_qc': 'Quality of the deconvolution of the run: masses per scan, charge states and the QScore distribution. A run far below others of the same sample is worth a look.',
    'internal_fragment_map': 'Internal fragments of the sequence, matched to the picked spectrum.',
    # Runs: identification
    'sequence_view': 'The sequence of the picked proteoform with its matched fragment ions (N- and C-terminal) and modifications. Many fragments along the sequence make an identification convincing. The pencil edits the sequence and its modifications.',
    'protein_table': 'Every protein identified in the run with its proteoforms, the best score and q-value. Picking a protein narrows the proteoforms and PrSMs to it.',
    'tag_table': 'The sequence tags (short stretches read from consecutive fragment masses) of the picked spectrum that matched the protein. Picking a tag marks its fragments.',
    'protein_group_table': 'The identified proteins of the run. Picking one lists its proteoforms.',
    'proteoform_table': 'The proteoforms of the picked protein: mass, modifications and the PrSMs that support them. Picking one lists its PrSMs.',
    'proteoform_prsm_table': 'The proteoform-spectrum matches (PrSMs) of the picked proteoform: scan, score and q-value. Picking one shows its spectrum and sequence.',
    'prsm_mass_table': 'The fragment masses of the picked PrSM and the ions they matched, with their mass errors.',
    'structure_view': 'A predicted 3D structure of the picked protein (AlphaFold), with the matched fragments mapped onto it. Drag to rotate.',
    'combined_spectrum': 'The picked spectrum with the sequence tags found in it drawn over the peaks, so it shows which stretch of the sequence each tag reads.',
    'id_fdr_plot': 'How many target and decoy PrSMs reach each score. Decoys show what wrong matches score; the q-value cut-off lets through the targets right of where decoys end.',
    'fragment_error_plot': 'The mass error of every matched fragment, in ppm. A tight band around zero means a calibrated instrument; a slope or offset points at a calibration problem.',
    'id_qc': 'Quality of the identification: score distribution, identifications over retention time, sequence coverage and fragment mass errors.',
    'proteoform_heatmap': 'Every PrSM by retention time (x) and proteoform mass (y). Proteoforms of one protein that differ by a modification lie on one vertical line.',
    'mass_shifts': 'How far the masses of the identified proteoforms lie from their unmodified sequence. Shifts shared by many PrSMs are modifications; picking one lists its PrSMs.',
    'go_enrichment': 'Gene Ontology (GO) terms over-represented among the proteins you pick, against all identified proteins of the run. Pick proteins with box or lasso select in the mode bar of the proteoform heatmap (Shift adds), or turn on the check boxes of the protein, proteoform or PrSM table with its check mark button and tick rows. A term with a small adjusted p-value is shared by more picked proteins than chance would give. With nothing picked, it counts the terms of all identified proteins. The annotations come from MyGene.info by UniProt accession; the test is one-sided Fisher’s exact test, Benjamini-Hochberg adjusted per GO category.',
    'reporter_table': 'The reporter ion intensities of every spectrum, one column per label (TMT or iTRAQ channel). Picking a spectrum shows its reporter ions.',
    'reporter_ions': 'The reporter ion region of the picked spectrum, one peak per label. Their heights are the relative amounts of the proteoform in the labelled samples.',
    # Comparisons
    'quant_summary': 'Key figures of the comparison: samples compared, entries tested, and how many are significantly higher or lower in the compared condition. The cut-offs are set with Edit comparison.',
    'quant_table': 'Every tested entry, most significant first: its log2 fold change (1 means twice as abundant), adjusted p-value and whether it counts as higher or lower. Picking one shows it per run; the check boxes pick a set for GO enrichment.',
    'quant_volcano': 'Each entry by its log2 fold change (x) against its significance, -log10 p-value (y). Far left or right and high up differ most reliably. A click picks one entry; All changed, Higher and Lower, or box and lasso select (Shift adds), pick a set for GO enrichment.',
    'quant_profile': 'The abundance of the picked entry in each run, coloured by condition. It shows whether a difference is consistent across replicates or comes from one outlier.',
    'quant_runs': 'The picked entry in each run: its intensity and the PrSM or feature behind it, so it shows in which runs it was measured and in which it was filled in.',
    'quant_pca': 'Principal component analysis places each sample by its abundances, so similar samples lie close together. Samples of one condition should cluster; one far from its condition may be an outlier.',
    'quant_heatmap': 'The abundance of the most variable entries in every sample, scaled per entry (z-score), with entries and samples ordered by similarity.',
    'quant_go': 'Gene Ontology (GO) terms over-represented among the entries you pick, against all tested entries. Pick entries with All changed, Higher or Lower above the volcano plot, box or lasso select in the mode bar of the volcano plot or the clustered heatmap (Shift adds), or the check boxes of the comparison table (its check mark button turns them on). A term with a small adjusted p-value is shared by more picked proteins than chance would give. With nothing picked, it counts the terms of all tested entries. The annotations come from MyGene.info by UniProt accession; the test is one-sided Fisher’s exact test, Benjamini-Hochberg adjusted per GO category.',
    'quant_id_counts': 'How many proteins, proteoforms and PrSMs each run identified. Runs of one sample should reach similar counts.',
    'quant_id_frequency': 'In how many runs each identification was found. Most should be found in every run; many found once point at noise or at runs that differ.',
    'quant_id_overlap': 'The identifications shared by each pair of runs. Replicates should share most of theirs.',
    'quant_id_table': 'Every identification of the compared runs with the runs it was found in and its best score. Picking one lists its PrSMs.',
    'quant_id_prsms': 'The PrSMs of the picked identification in every run. Picking one shows its sequence and spectrum.',
    'quant_prsm_sequence': 'The sequence of the picked PrSM with its matched fragment ions, from its own run.',
    'quant_prsm_spectrum': 'The spectrum of the picked PrSM with its matched fragments marked, from its own run.',
    'quant_reporter_intensities': 'The distribution of reporter ion intensities per label. Labels of equal loading should reach similar heights.',
    'quant_reporter_ratios': 'The ratio of each label to the reference label, per spectrum. Most proteoforms do not change, so the ratios should centre at one.',
    'quant_reporter_labels': 'How many labels were found per spectrum. Spectra missing labels give incomplete ratios.',
    'quant_reporter_spectra': 'The reporter ions of every spectrum of the picked run. Picking one shows its reporter ions.',
    'quant_reporter_ions': 'The reporter ion region of the picked spectrum, one peak per label.',
    'quant_feature_map': 'The consensus features (a proteoform matched across runs) by retention time and mass. Picking one shows it per run.',
    'quant_feature_intensity': 'The intensity distribution of the features per run, detected and matched apart. Matched values should not sit far below the detected ones.',
    'quant_feature_table': 'Every consensus feature with its mass, retention time and intensity in each run. Picking one shows it per run.',
    'quant_feature_profile': 'The picked consensus feature in each run: intensity and retention time, so a shift between runs shows.',
    'quant_totals': 'The summed intensity of each sample. Large differences mean unequal loading, which normalization corrects.',
    'quant_correlation': 'How closely the abundances of each pair of samples agree (Pearson). Replicates should correlate more with each other than with other conditions.',
    'quant_cv': 'The coefficient of variation of each entry within its condition. Low values mean reproducible replicates.',
    'quant_match_quality': 'How well features were matched between runs: the retention time and mass differences of matched features. Tight distributions mean reliable matches.',
    'quant_missing_values': 'How many abundances are missing per sample and entry. Missing values mostly belong to low-abundant entries.',
    'quant_imputed_values': 'The values filled in for missing abundances, against the measured ones. Filled-in values should sit at the low end, where the missing ones come from.',
    'quant_normalization': 'The abundance distribution of each sample before and after normalization. After it, the medians should line up.',
    'quant_pvalues': 'The distribution of the comparison\'s p-values. Flat with a peak near zero means some entries differ; flat only means none do; a peak elsewhere hints at a test that does not fit the data.',
}


@contextmanager
def drawing(name):
    """Draw component ``name``: the first ``help()`` meanwhile returns its text,
    so a component made of several plots shows it once, on the first."""
    if _current.get() is not None:
        # Inside another component (a run's cell drawn in a comparison): its name stays
        yield
        return
    token = _current.set([name, False])
    try:
        yield
    finally:
        _current.reset(token)


def help(name=None):  # pylint: disable=redefined-builtin
    """The explanation of ``name`` (default: the component being drawn, once per component)."""
    if name is not None:
        return COMPONENTS.get(name)
    current = _current.get()
    if current is None or current[1]:
        return None
    current[1] = True
    return COMPONENTS.get(current[0])

