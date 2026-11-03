"""What a top-down analysis does, in the words of the bench: whether the
proteoforms are identified, for which kind of sample, and how they are
quantified. The workflow's Describe experiment page sets these, its run
follows them and the Viewer picks its default layout from them.

Streamlit-free.
"""

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
    NO_QUANT: 'Find and identify the proteoforms only.',
    LABEL_FREE: (
        'Each sample is its own run. The abundance of a proteoform is the area of its signal, '
        'and the runs are linked by mass and retention time to compare conditions.'
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
