"""The Viewer's views: deconvolution and identification components in one
combined view, the views of a comparison of runs, and the default each
experiment starts from.

One view bar lists every view, those of a run and those of a comparison
(Label-free, TMT, iTRAQ; src/render/quant_view.py). A view the shown
experiment cannot fill stays listed, locked, and offers what can fill it.

Every earlier default is a view. The ready-made views are the built-in
defaults: Deconvolution (with the sequence view once a sequence is set), Native
and Lysate, each with the reporter ions of a run with isobaric labels. A layout
saved in one of the two former Viewers becomes "My deconvolution layout" or "My
identification layout" and stays the default of its kind; one saved with the
former Edit layout dialog becomes the current view. Until a view is picked,
each experiment starts from the default that fits its run.
"""

from src.render.layout_editor import get_saved_layout
from src.render.layout_editor.catalog import (
    DECONV_DEFAULT_LAYOUT, TNT_DEFAULT_LAYOUTS, VIEWER_COMPONENTS, VIEWER_SEQUENCE_COMPONENTS,
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
    """What the views editor offers: every component of a run and of a
    comparison, plus the sequence components once a sequence is set in the
    sequence view."""
    if has_sequence(file_manager):
        return VIEWER_COMPONENTS + VIEWER_SEQUENCE_COMPONENTS
    return VIEWER_COMPONENTS


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
    """The ready-made views: the deconvolution default, the identification
    defaults by sample type, then the views of a comparison of runs."""
    from src.render import quant_view
    return [
        {'name': DECONVOLUTION_VIEW, 'label': 'Deconvolution', 'layout': deconvolution_default(file_manager),
         'help': 'Spectra and their deconvolved masses: heatmap, scan and mass tables, spectra, 3D S/N plot.'},
        {'name': SAMPLE_VIEWS['Native'], 'label': 'Native', 'layout': TNT_DEFAULT_LAYOUTS['Native'],
         'help': 'Identifications of native or purified proteins: each proteoform in depth, with its '
                 'sequence, structure and tags.'},
        {'name': SAMPLE_VIEWS['Lysate'], 'label': 'Lysate', 'layout': TNT_DEFAULT_LAYOUTS['Lysate'],
         'help': 'Identifications of a lysate: proteins, proteoforms and PrSMs, browsed protein by protein.'},
    ] + quant_view.presets()


# The layouts the former Viewers saved become views of the user's own, named
# after the Viewer they come from
LEGACY_VIEWS = {False: ('custom-deconvolution', 'My deconvolution layout'),
                True: ('custom-identification', 'My identification layout')}


def default_view(identified, sample_type, file_manager=None):
    """The view that is the default of a run: the layout the former Viewer of
    its kind saved, else the ready-made view of its analysis."""
    legacy = get_saved_layout(file_manager, LEGACY_LAYOUT_DATASET[identified]) if file_manager else None
    if legacy:
        from src.layout_manager import match_preset
        return (match_preset([legacy[0][0]], presets(file_manager), components(file_manager))
                or LEGACY_VIEWS[identified][0])
    return SAMPLE_VIEWS.get(sample_type, SAMPLE_VIEWS['Native']) if identified else DECONVOLUTION_VIEW


def adopt_saved_layout(file_manager, workspace):
    """Move the layouts saved before the view bar into the views file, once.

    Every earlier default stays reachable: the ready-made views are the
    built-in defaults, a layout saved in a former Viewer (deconvolution or
    identification) becomes "My deconvolution layout" or "My identification
    layout", and one saved with Edit layout becomes the current view (the
    ready-made view it matches, or "My layout"). Without an Edit layout one,
    nothing counts as picked, so every run keeps the default of its kind.
    """
    from src.layout_manager import load_views, match_preset
    path = views_path(workspace)
    if path.exists():
        return
    comps = components(file_manager)
    ready = presets(file_manager)
    own = []
    for identified in (False, True):
        legacy = get_saved_layout(file_manager, LEGACY_LAYOUT_DATASET[identified])
        if not legacy or match_preset([legacy[0][0]], ready, comps):
            continue
        name, label = LEGACY_VIEWS[identified]
        own.append({'name': name, 'label': label, 'rows': legacy[0][0], 'in_bar': True})
    import json
    data = {'views': [{'name': p['name'], 'in_bar': True} for p in ready]}
    layout = saved(file_manager)
    if layout:
        # The view the saved layout shows, as the template migrates it
        from src.layout_manager import save_layout
        save_layout(path, layout[0], layout[1])
        state = load_views(path, comps, ready, deconvolution_default(file_manager))
        data = json.loads(path.read_text())
        data['current'] = state['current']
        data['views'] = [{'name': v['name'], 'in_bar': v['in_bar'],
                          **({} if v['preset'] else {'label': v['label'], 'rows': v['rows']})}
                         for v in state['views']]
        own = [v for v in own if all(v['rows'] != w['rows'] for w in state['views'])]
    elif not own:
        return
    data['views'] += own
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=1))


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
