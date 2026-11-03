"""The one layout of the Viewer: deconvolution and identification components
in one combined view, and the default each experiment starts from.

Once a layout is saved (Edit layout), every experiment shows it. Until then
each experiment starts from the default that fits its run: an identified run
from the identification default of its sample type (Native or Lysate), a run
deconvolved only from the deconvolution default. A layout saved in one of the
two former Viewers stays the default of its kind.
"""

from src.render.layout_editor import get_saved_layout
from src.render.layout_editor.catalog import (
    DECONV_DEFAULT_LAYOUT, TNT_DEFAULT_LAYOUTS, TOPDOWN_COMPONENTS, TOPDOWN_SEQUENCE_COMPONENTS,
)

LAYOUT_DATASET = 'topdown_layout'
# A run quantified with isobaric labels: its reporter ions, shown below the
# first component that selects a spectrum
REPORTER_ROW = ('reporter_table', 'reporter_ions')
REPORTER_AFTER = ('proteoform_prsm_table', 'protein_table', 'scan_table')
# The layouts the two former Viewers saved, by "the run is identified"
LEGACY_LAYOUT_DATASET = {True: 'flashtnt_layout', False: 'layout'}


def has_sequence(file_manager):
    return file_manager.result_exists('sequence', 'sequence')


def components(file_manager):
    """What Edit layout offers: every component of both kinds, plus the
    sequence components of the Sequence Input page once a sequence is saved."""
    if has_sequence(file_manager):
        return TOPDOWN_COMPONENTS + TOPDOWN_SEQUENCE_COMPONENTS
    return TOPDOWN_COMPONENTS


def deconvolution_default(file_manager):
    if has_sequence(file_manager):
        return DECONV_DEFAULT_LAYOUT + [['sequence_view']]
    return DECONV_DEFAULT_LAYOUT


def saved(file_manager):
    """(experiments, side_by_side) saved with Edit layout, or None."""
    return get_saved_layout(file_manager, LAYOUT_DATASET)


def has_legacy(file_manager, identified):
    """True when the former Viewer of the run's kind saved a layout."""
    return get_saved_layout(file_manager, LEGACY_LAYOUT_DATASET[identified]) is not None


def with_reporters(rows, after):
    """`rows` with the reporter ion row inserted after the first row holding
    one of `after` (a scan-selecting component), else at the end."""
    rows = [list(row) for row in rows]
    if any(name in row for row in rows for name in REPORTER_ROW):
        return rows
    for i, row in enumerate(rows):
        if any(name in row for name in after):
            return rows[:i + 1] + [list(REPORTER_ROW)] + rows[i + 1:]
    return rows + [list(REPORTER_ROW)]


def default_rows(file_manager, identified, sample_type, index=0, labelled=False):
    """Rows of experiment ``index`` without a saved layout: the former
    Viewer's saved layout of the run's kind, else the built-in default of the
    analysis (deconvolution, native or lysate identification), with the
    reporter ions of a run quantified with isobaric labels."""
    legacy = get_saved_layout(file_manager, LEGACY_LAYOUT_DATASET[identified])
    if legacy:
        experiments = legacy[0]
        rows = experiments[index] if index < len(experiments) else experiments[0]
    elif identified:
        rows = TNT_DEFAULT_LAYOUTS[sample_type]
    else:
        rows = deconvolution_default(file_manager)
    if labelled:
        rows = with_reporters(rows, REPORTER_AFTER)
    return rows


def rows(file_manager, identified, sample_type, index=0, labelled=False):
    """Rows experiment ``index`` shows."""
    layout = saved(file_manager)
    if layout:
        experiments = layout[0]
        return experiments[index] if index < len(experiments) else experiments[0]
    return default_rows(file_manager, identified, sample_type, index, labelled)
