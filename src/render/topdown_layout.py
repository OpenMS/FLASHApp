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


def default_rows(file_manager, identified, sample_type, index=0):
    """Rows of experiment ``index`` without a saved layout: the former
    Viewer's saved layout of the run's kind, else the built-in default."""
    legacy = get_saved_layout(file_manager, LEGACY_LAYOUT_DATASET[identified])
    if legacy:
        experiments = legacy[0]
        return experiments[index] if index < len(experiments) else experiments[0]
    if identified:
        return TNT_DEFAULT_LAYOUTS[sample_type]
    return deconvolution_default(file_manager)


def rows(file_manager, identified, sample_type, index=0):
    """Rows experiment ``index`` shows."""
    layout = saved(file_manager)
    if layout:
        experiments = layout[0]
        return experiments[index] if index < len(experiments) else experiments[0]
    return default_rows(file_manager, identified, sample_type, index)
