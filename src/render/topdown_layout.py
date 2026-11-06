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
    sequence components once a sequence is set in the sequence view."""
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


# ---- Views: the Viewer's view bar (the template's layout manager) ----------
#
# The ready-made views are the three defaults. The user's own views, their
# order and the runs shown are kept in <workspace>/layouts/viewer.json. Until a
# view is picked or saved there, every experiment keeps the default of its run
# (``default_rows``) and the bar marks the view of the first run's kind.

DECONVOLUTION_VIEW = 'deconvolution'
SAMPLE_VIEWS = {'Native': 'native', 'Lysate': 'lysate'}


def views_path(workspace):
    from pathlib import Path
    return Path(workspace, 'layouts', 'viewer.json')


def presets(file_manager):
    """The ready-made views: the deconvolution default, then the identification
    defaults by sample type."""
    return [
        {'name': DECONVOLUTION_VIEW, 'label': 'Deconvolution', 'layout': deconvolution_default(file_manager),
         'help': 'Spectra and their deconvolved masses: heatmap, scan and mass tables, spectra, 3D S/N plot.'},
        {'name': SAMPLE_VIEWS['Native'], 'label': 'Native', 'layout': TNT_DEFAULT_LAYOUTS['Native'],
         'help': 'Identifications of native or purified proteins: each proteoform in depth, with its '
                 'sequence, structure and tags.'},
        {'name': SAMPLE_VIEWS['Lysate'], 'label': 'Lysate', 'layout': TNT_DEFAULT_LAYOUTS['Lysate'],
         'help': 'Identifications of a lysate: proteins, proteoforms and PrSMs, browsed protein by protein.'},
    ]


def default_view(identified, sample_type):
    """The ready-made view that is the default of a run."""
    return SAMPLE_VIEWS.get(sample_type, SAMPLE_VIEWS['Native']) if identified else DECONVOLUTION_VIEW


def adopt_saved_layout(file_manager, workspace):
    """Move a layout saved with the former Edit layout into the views file, once.

    The views file then shows it as the view it matches, or as "My layout".
    """
    from src.layout_manager import save_layout
    path = views_path(workspace)
    if path.exists():
        return
    layout = saved(file_manager)
    if layout:
        save_layout(path, layout[0], layout[1])


def view_rows(file_manager, workspace, identified, sample_type, index=0, labelled=False):
    """Rows experiment ``index`` shows, or None while no view is picked or saved.

    A ready-made view gets the reporter ions of a run quantified with isobaric
    labels, as the defaults do.
    """
    from src.layout_manager import load_layout, load_views
    path = views_path(workspace)
    if load_layout(path) is None:
        return None
    state = load_views(path, components(file_manager), presets(file_manager), deconvolution_default(file_manager))
    current = next(v for v in state['views'] if v['name'] == state['current'])
    rows = [list(r) for r in current['rows']]
    if labelled and current['preset']:
        rows = with_reporters(rows, REPORTER_AFTER)
    return rows


def experiments_shown(workspace):
    """(number of experiments, side by side) of the views file, or None."""
    from src.layout_manager import load_layout
    stored = load_layout(views_path(workspace))
    if stored is None:
        return None
    return len(stored[0]), stored[1]
