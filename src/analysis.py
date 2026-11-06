"""What a top-down analysis does, in the words of the bench: whether the
proteoforms are identified, for which kind of sample, and how they are
quantified. The workflow's Describe experiment page sets these, its run
follows them and the Viewer picks its default layout from them.

Streamlit-free.
"""

import re
from pathlib import Path

from src import sample_type
from src.quant import isobaric

# Widget keys of the workflow's parameters (Describe experiment)
IDENTIFY_KEY = 'identify'
QUANT_KEY = 'quantification'
PLEX_KEY = 'isobaric-plex'
# Before the quantification choice, comparing runs was a switch of its own
COMPARE_KEY = 'compare-runs'

NO_QUANT = 'none'
LABEL_FREE = 'label-free'
ISOBARIC = 'isobaric'
QUANT_LABELS = {
    NO_QUANT: 'None',
    LABEL_FREE: 'Label-free: compare runs',
    ISOBARIC: 'Isobaric labels (TMT, iTRAQ)',
}
QUANT_HELP = {
    NO_QUANT: 'Deconvolve the runs, and identify the proteoforms when that is on, without comparing samples.',
    LABEL_FREE: (
        'Each sample is its own run. The abundance of a proteoform is the area of its signal, '
        'and the runs are linked by mass and retention time to compare conditions. Works with '
        'or without identification: without it, the masses are compared.'
    ),
    ISOBARIC: (
        'The samples were labelled (TMT or iTRAQ) and pooled into one run, or one set of fractions. '
        'The reporter ion of each label in a fragmentation spectrum gives the abundance in its sample.'
    ),
}

SAMPLE_LABELS = {
    sample_type.NATIVE: 'Native or purified proteins',
    sample_type.LYSATE: 'Cell lysate or complex mixture',
}
SAMPLE_HELP = {
    sample_type.NATIVE: 'A few proteins per sample: native MS, a purified protein or targeted runs.',
    sample_type.LYSATE: 'Many proteins, separated by liquid chromatography before MS/MS.',
}


def quantification(params):
    """NO_QUANT, LABEL_FREE or ISOBARIC from a workflow's parameters."""
    value = params.get(QUANT_KEY)
    if value in QUANT_LABELS:
        return value
    return LABEL_FREE if params.get(COMPARE_KEY, False) else NO_QUANT


def plex(params):
    value = params.get(PLEX_KEY, isobaric.DEFAULT_PLEX)
    return value if value in isobaric.PLEXES else isobaric.DEFAULT_PLEX


def name(identify, quant, plex_type=None):
    """The analysis in a few words, e.g. 'Deconvolution + identification +
    TMT 10-plex quantification'."""
    parts = ['Deconvolution']
    if identify:
        parts.append('identification')
    if quant == LABEL_FREE:
        parts.append('label-free quantification')
    elif quant == ISOBARIC:
        parts.append(f'{isobaric.label(plex_type or isobaric.DEFAULT_PLEX)} quantification')
    return ' + '.join(parts)


def _runs(count):
    return 'one run' if count == 1 else f'{count} runs'


def summary(runs, identify, sample, database, quant, plex_type=None):
    """What a start of the workflow will do, in one line."""
    if not runs:
        return 'Select the MS data to analyse.'
    parts = ['Deconvolution']
    if identify:
        db = f'**{Path(database).name}**' if database else 'no database selected'
        settings = 'native' if sample == sample_type.NATIVE else 'lysate'
        parts.append(f'identification ({settings} settings, {db})')
    if quant == LABEL_FREE:
        parts.append('label-free quantification')
    elif quant == ISOBARIC:
        parts.append(f'{isobaric.label(plex_type or isobaric.DEFAULT_PLEX)} quantification')
    text = ' + '.join(parts) + f' of {_runs(runs)}'
    if quant == ISOBARIC and runs > 1:
        text += ' (fractions of one labelled set)'
    text += '.'
    if quant == LABEL_FREE and runs < 2:
        text += ' Label-free quantification needs at least two runs.'
    return text


# The analysis step each tool or workflow setting of a preset belongs to.
# A tool not listed here (FLASHDeconv) runs in every analysis.
IDENTIFICATION = 'identification'
LABEL_FREE_STEP = 'label-free quantification'
ISOBARIC_STEP = 'isobaric quantification'
PRESET_STEPS = {
    'FLASHTnT': IDENTIFICATION,
    'FLASHAnnotator': IDENTIFICATION,
    'DecoyDatabase': IDENTIFICATION,
    'FLASHQuantConsensus': LABEL_FREE_STEP,
}
GENERAL_STEPS = {
    'few_proteins': IDENTIFICATION,
    PLEX_KEY: ISOBARIC_STEP,
}
# FLASHDeconv parameters that act in one kind of quantification only
DECONV_PREFIX_STEPS = {'iq:': ISOBARIC_STEP}


def active_steps(identify, quant):
    """The optional steps a start of the workflow runs."""
    steps = set()
    if identify:
        steps.add(IDENTIFICATION)
    if quant == LABEL_FREE:
        steps.add(LABEL_FREE_STEP)
    elif quant == ISOBARIC:
        steps.add(ISOBARIC_STEP)
    return steps


def preset_steps(preset):
    """The optional steps whose settings a preset (an entry of presets.json) sets."""
    steps = set()
    for section, values in (preset or {}).items():
        if section == '_description' or not isinstance(values, dict):
            continue
        if section == '_general':
            steps |= {GENERAL_STEPS[k] for k in values if k in GENERAL_STEPS}
        elif section in PRESET_STEPS:
            steps.add(PRESET_STEPS[section])
        elif section == 'FLASHDeconv':
            steps |= {s for k in values for p, s in DECONV_PREFIX_STEPS.items() if k.startswith(p)}
    return steps


def preset_fit(preset, identify, quant):
    """The steps a preset sets that the chosen analysis does not run, sorted.
    Empty when the preset fits the analysis. A preset that sets nothing but
    settings of inactive steps has no effect on a start."""
    return sorted(preset_steps(preset) - active_steps(identify, quant))


def preset_warning(name, preset, identify, quant):
    """A warning for a preset that does not fit the chosen analysis, else None."""
    unused = preset_fit(preset, identify, quant)
    if not unused:
        return None
    steps = ' and '.join(unused)
    switch = []
    if IDENTIFICATION in unused:
        switch.append('switch on **Identify proteoforms**')
    if LABEL_FREE_STEP in unused:
        switch.append(f"choose **Quantification: {QUANT_LABELS[LABEL_FREE]}**")
    if ISOBARIC_STEP in unused:
        switch.append(f"choose **Quantification: {QUANT_LABELS[ISOBARIC]}**")
    return (f'**{name}** sets {steps} settings, but this analysis does not include {steps}: '
            f'they have no effect on this run. To use them, {" and ".join(switch)}.')


def preset_note(key, preset, identify, quant):
    """What to say about a preset for the chosen analysis: None when it fits,
    else (level, text) with level 'warning' or 'info'.

    A sample preset (Native, Lysate) describes the sample, which holds for
    every analysis: deconvolving LC-MS/MS (DDA) runs, or comparing them
    label-free, without identification is a supported analysis. Its
    identification settings then wait for identification to be switched on,
    which is a note, not a warning. Any other setting a preset makes for a
    step the analysis does not run is a warning."""
    unused = preset_fit(preset, identify, quant)
    if not unused:
        return None
    if key in SAMPLE_LABELS and unused == [IDENTIFICATION]:
        return 'info', (
            f'**{SAMPLE_LABELS[key]}** sets how proteoforms are identified. Without identification '
            'the runs are deconvolved (and compared, with quantification) the same way for either sample type; '
            'these settings apply once **Identify proteoforms** is on.'
        )
    return 'warning', preset_warning(SAMPLE_LABELS.get(key, key), preset, identify, quant)


_STAMP = re.compile(r'(\d{8}-\d{6})$')


def newest_result(dataset_ids):
    """The dataset a start made last: runs are ``<file>_<YYYYmmdd-HHMMSS>``,
    comparisons ``comparison_<YYYYmmdd-HHMMSS>``; a comparison wins a tie, as it
    is made after its runs. None when no id carries a time (only examples)."""
    best = None
    for dataset_id in dataset_ids:
        match = _STAMP.search(str(dataset_id))
        if not match:
            continue
        rank = (match.group(1), str(dataset_id).startswith('comparison_'))
        if best is None or rank > best[0]:
            best = (rank, dataset_id)
    return best[1] if best else None
