"""Comparisons of runs in the Viewer.

A workflow run that quantified two or more files stores a comparison
(src/quant/session.py) next to its runs. The Viewer lists it in the experiment
selector as "Comparison: ..."; selecting it shows a grid of comparison
components (src/render/quant_cells.py). A view bar of its own picks the grid:
the ready-made views **Quantification**, **Identification** and the QC of the
comparison's kind (**QC: label-free** for runs compared one by one, **QC: TMT,
iTRAQ** for labelled samples), and the user's own (+ and ✏️), in
`<workspace>/layouts/comparison.json`. A comparison of runs that were not
identified shows Identification locked, with + to add identification to its
runs. Until one is picked, a comparison shows Quantification. The layout saved with the
former Edit layout dialog becomes "My comparison layout" and the default of
every comparison. **Edit comparison** opens a window with the conditions of
the runs and the analysis settings.
"""

from pathlib import Path

import streamlit as st

from src.quant import isobaric
from src.quant.session import LEVEL_NAMES, QuantSession, comparisons
from src.render.layout_editor import get_saved_layout
from src.render.layout_editor.catalog import (
    IDENTIFICATION, QUANT_COMPONENTS, QUANT_DEFAULT_LAYOUT, QUANT_DIAGNOSTICS, QUANT_IDENTIFICATION_LAYOUT,
    QUANT_ISOBARIC_QC_LAYOUT, QUANT_LABEL_FREE_QC_LAYOUT, QUANT_QUANTIFICATION_LAYOUT,
)
from src.render.layout_editor.dialogs import _own_file_manager

# One layout serves every comparison. Before the workflows were merged, a
# comparison of identified runs had its own (migrate_layout).
LAYOUT_DATASET = 'quant_layout'
LEGACY_LAYOUT_DATASET = 'flashtnt_quant_layout'


def migrate_layout(file_manager):
    """Keep the layout saved for comparisons of identified runs when no other
    comparison layout was saved; drop it either way."""
    if not file_manager.result_exists(LEGACY_LAYOUT_DATASET, 'layout'):
        return
    if not file_manager.result_exists(LAYOUT_DATASET, 'layout'):
        legacy = file_manager.get_results(LEGACY_LAYOUT_DATASET, 'layout')['layout']
        file_manager.store_data(LAYOUT_DATASET, 'layout', legacy)
    file_manager.remove_results(LEGACY_LAYOUT_DATASET)


def comparison_names(file_manager, tool=None):
    """Display name -> dataset id of the comparisons (of `tool`, default all), newest first."""
    names = {}
    for dataset_id in comparisons(file_manager, tool):
        manifest = file_manager.get_results(dataset_id, ['quant_manifest'])['quant_manifest']
        stamp = dataset_id.rsplit('_', 1)[-1]
        when = f'{stamp[:4]}-{stamp[4:6]}-{stamp[6:8]} {stamp[9:11]}:{stamp[11:13]}' if len(stamp) >= 13 else stamp
        if manifest.get('kind') == 'isobaric':
            what = f"{isobaric.label(manifest.get('plex', ''))}, {len(manifest['runs'])} samples"
        else:
            what = f"{len(manifest['runs'])} runs"
        # A comparison given a name of its own (the examples) keeps it
        custom = file_manager.get_display_name(dataset_id)
        names[custom if custom != dataset_id else f"Comparison: {what} ({when})"] = dataset_id
    return names


def added_comparison(file_manager, workflow_dir):
    """The newest comparison of a start that processed more files for earlier
    runs (Process more files) and compared them with the new ones, else None."""
    import json

    try:
        info = json.loads(Path(workflow_dir, 'last-follow-up.json').read_text(encoding='utf-8'))
    except (OSError, ValueError):
        return None
    earlier = set(info.get('adds_to') or []) if isinstance(info, dict) else set()
    if not earlier:
        return None
    for dataset_id in comparisons(file_manager):
        manifest = file_manager.get_results(dataset_id, ['quant_manifest'])['quant_manifest']
        if earlier <= set(manifest.get('runs') or []):
            return dataset_id
    return None


def comparison_runs(file_manager, dataset_id):
    """The dataset ids of the runs a comparison compares (the fractions of
    an isobaric one, whose ``runs`` are its labelled samples)."""
    manifest = file_manager.get_results(dataset_id, ['quant_manifest'])['quant_manifest']
    if manifest.get('kind') == 'isobaric':
        return list(manifest.get('fractions') or [])
    return list(manifest.get('runs') or [])


def comparison_tool(file_manager, dataset_id):
    """'flashtnt' for a comparison of identified runs, else 'flashdeconv'."""
    manifest = file_manager.get_results(dataset_id, ['quant_manifest'])['quant_manifest']
    return manifest.get('tool') or 'flashdeconv'


# The plots behind the steps of the comparison were one component with a
# switch between them; layouts saved then name it.
LEGACY_DIAGNOSTICS = 'quant_diagnostics'


def split_diagnostics(rows):
    """``rows`` with the former diagnostics component replaced by its four
    plots, two per row, below what else its row held."""
    out = []
    for row in rows:
        if LEGACY_DIAGNOSTICS not in row:
            out.append(row)
            continue
        rest = [name for name in row if name != LEGACY_DIAGNOSTICS]
        if rest:
            out.append(rest)
        out.extend(list(r) for r in QUANT_DIAGNOSTICS)
    return out


def _rename_views(data):
    """Point the saved views at the ready-made views that replaced the
    former Label-free and TMT, iTRAQ views (both now Quantification).
    Returns True when something changed."""
    changed = False
    if data.get('current') in LEGACY_VIEWS:
        data['current'] = LEGACY_VIEWS[data['current']]
        changed = True
    views, seen = [], set()
    for view in data.get('views') or []:
        if isinstance(view, dict) and view.get('name') in LEGACY_VIEWS and 'rows' not in view:
            view = {**view, 'name': LEGACY_VIEWS[view['name']]}
            changed = True
        name = view.get('name') if isinstance(view, dict) else None
        if name is not None and name in seen:
            changed = True
            continue
        seen.add(name)
        views.append(view)
    if changed and 'views' in data:
        data['views'] = views
    return changed


# The order of the ready-made views before Quantification came first
ORDER_VERSION = 2


def _order_marker(path):
    """Beside the views file: its ready-made views have today's order (the
    view bar keeps no keys of its own in the file when it saves it)."""
    return Path(path).with_name(Path(path).stem + f'.order-{ORDER_VERSION}')


def _reorder_presets(data):
    """Put the ready-made views of a views file saved before Quantification
    came first in today's order, in the places they took. Returns True when
    the file changes."""
    if not isinstance(data.get('views'), list):
        return False
    order = [p['name'] for p in presets()]
    views = data['views']
    slots = [i for i, v in enumerate(views) if isinstance(v, dict) and v.get('name') in order]
    ranked = sorted((views[i] for i in slots), key=lambda v: order.index(v['name']))
    if [views[i] for i in slots] == ranked:
        return False
    for i, view in zip(slots, ranked):
        views[i] = view
    return True


def migrate_views(path):
    """Split the former diagnostics component and rename the former
    ready-made views in the saved comparison views, once."""
    import json
    try:
        data = json.loads(Path(path).read_text(encoding='utf-8'))
    except (OSError, ValueError):
        return
    if not isinstance(data, dict):
        return
    changed = _rename_views(data)
    marker = _order_marker(path)
    if not marker.exists():
        changed = _reorder_presets(data) | changed
        marker.touch()
    if LEGACY_DIAGNOSTICS in json.dumps(data):
        for view in data.get('views') or []:
            if isinstance(view, dict) and isinstance(view.get('rows'), list):
                view['rows'] = split_diagnostics(view['rows'])
        if isinstance(data.get('layout'), list):
            data['layout'] = [split_diagnostics(experiment) for experiment in data['layout']]
        changed = True
    if changed:
        Path(path).write_text(json.dumps(data, indent=1))


def layout(file_manager):
    """The layout saved with the former Edit layout dialog, else the default."""
    saved = get_saved_layout(file_manager, LAYOUT_DATASET)
    return split_diagnostics(saved[0][0]) if saved else QUANT_DEFAULT_LAYOUT


IDENTIFICATION_VIEW = 'identification'
QUANTIFICATION_VIEW = 'quantification'
LABEL_FREE_QC_VIEW = 'qc-label-free'
ISOBARIC_QC_VIEW = 'qc-isobaric'
# The ready-made views before the four above, as saved views may name them
LEGACY_VIEWS = {'label-free': QUANTIFICATION_VIEW, 'isobaric': QUANTIFICATION_VIEW}
# The layout saved with the former Edit layout dialog, as a view of one's own
LEGACY_VIEW = ('custom-comparison', 'My comparison layout')


def views_path(workspace):
    return Path(workspace, 'layouts', 'comparison.json')


def presets(isobaric_comparison=None):
    """The ready-made comparison views, in the order of the bar: what changed
    between the conditions, what the runs identified, and the quality control
    of a label-free and of an isobaric (TMT, iTRAQ) comparison.
    ``isobaric_comparison`` True or False keeps only the QC of that kind."""
    views = [
        {'name': QUANTIFICATION_VIEW, 'label': 'Quantification', 'layout': QUANT_QUANTIFICATION_LAYOUT,
         'help': 'What changed between the conditions: volcano plot, comparison table, each entry run by run, '
                 'sample similarity, the clustered heatmap and GO enrichment.'},
        {'name': IDENTIFICATION_VIEW, 'label': 'Identification', 'layout': QUANT_IDENTIFICATION_LAYOUT,
         'help': 'What each run identified and what the runs share, then every proteoform with its status run '
                 'by run; pick one to list its PrSMs in every run, and a PrSM to see its sequence view and '
                 'spectrum. Needs no conditions.'},
        {'name': LABEL_FREE_QC_VIEW, 'label': 'QC: label-free', 'layout': QUANT_LABEL_FREE_QC_LAYOUT,
         'help': 'Quality control of runs measured one by one: the consensus features (feature map, intensity '
                 'of detected and matched values, every feature run by run), match between runs, how the runs '
                 'were loaded, how well replicates agree, missing values and p-values.'},
        {'name': ISOBARIC_QC_VIEW, 'label': 'QC: TMT, iTRAQ', 'layout': QUANT_ISOBARIC_QC_LAYOUT,
         'help': 'Quality control of the labelled samples of one set: the reporter ions (intensity and ratio '
                 'per label, labels found per spectrum, every labelled spectrum and its ions), whether the '
                 'channels were loaded alike, how they correlate and cluster, missing values and p-values.'},
    ]
    if isobaric_comparison is None:
        return views
    other = LABEL_FREE_QC_VIEW if isobaric_comparison else ISOBARIC_QC_VIEW
    return [view for view in views if view['name'] != other]


def comparison_presets(file_manager, dataset_id):
    """The ready-made views of comparison ``dataset_id``: the QC of its kind only."""
    return presets(is_isobaric(file_manager, dataset_id))


def available(file_manager, dataset_id):
    """The analyses comparison ``dataset_id`` has: identification when its
    runs were identified (src/render/missing_analysis.py)."""
    return {IDENTIFICATION} if comparison_tool(file_manager, dataset_id) == 'flashtnt' else set()


def is_isobaric(file_manager, dataset_id):
    """True for a comparison of labelled samples (TMT, iTRAQ)."""
    manifest = file_manager.get_results(dataset_id, ['quant_manifest'])['quant_manifest']
    return manifest.get('kind') == 'isobaric'


def _legacy_rows(file_manager):
    saved = get_saved_layout(file_manager, LAYOUT_DATASET)
    return split_diagnostics(saved[0][0]) if saved else None


def default_view(file_manager, dataset_id=None):
    """The view a comparison shows until one is picked: the layout saved
    before the view bar, else Quantification."""
    del dataset_id
    from src.layout_manager import match_preset
    legacy = _legacy_rows(file_manager)
    if legacy:
        return match_preset([legacy], presets(), QUANT_COMPONENTS) or LEGACY_VIEW[0]
    return QUANTIFICATION_VIEW


def adopt_saved_layout(file_manager, workspace):
    """Make the layout saved with the former Edit layout dialog a view of the
    user's own, once (unless it is one of the ready-made views)."""
    from src.layout_manager import match_preset
    path = views_path(workspace)
    legacy = _legacy_rows(file_manager)
    if path.exists() or not legacy or match_preset([legacy], presets(), QUANT_COMPONENTS):
        return
    import json
    name, label = LEGACY_VIEW
    _order_marker(path).parent.mkdir(parents=True, exist_ok=True)
    _order_marker(path).touch()
    data = {'views': [{'name': p['name'], 'in_bar': True} for p in presets()]
            + [{'name': name, 'label': label, 'rows': legacy, 'in_bar': True}]}
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=1))


def fit_current(path, shown):
    """Keep the current view of the views file one the comparison offers.

    The views file serves every comparison, but each offers only the QC of
    its kind (``shown``, its ready-made views). A picked QC of the other kind
    is parked and Quantification shown; it comes back for a comparison of
    its kind, unless another view was picked meanwhile."""
    import json
    try:
        data = json.loads(Path(path).read_text(encoding='utf-8'))
    except (OSError, ValueError):
        return
    if not isinstance(data, dict):
        return
    names = [p['name'] for p in shown]
    hidden = {p['name'] for p in presets()} - set(names)
    parked_path = Path(path).with_name(Path(path).stem + '.parked.json')
    try:
        parked = json.loads(parked_path.read_text(encoding='utf-8'))
    except (OSError, ValueError):
        parked = None
    current = data.get('current')
    changed = False
    if isinstance(parked, dict):
        if current != parked.get('replaced_by'):
            parked_path.unlink(missing_ok=True)  # another view was picked meanwhile
        elif parked.get('view') in names:
            data['current'] = parked['view']
            parked_path.unlink(missing_ok=True)
            changed = True
    if data.get('current') in hidden:
        parked_path.write_text(json.dumps({'view': data['current'], 'replaced_by': names[0]}))
        data['current'] = names[0]
        changed = True
    if changed:
        rows = next(p['layout'] for p in presets() if p['name'] == data['current']) \
            if data['current'] in {p['name'] for p in presets()} else None
        if rows is not None:
            data['layout'] = [rows]
        Path(path).write_text(json.dumps(data, indent=1))


def view_rows(file_manager, workspace, dataset_id):
    """The rows a comparison shows: the picked view, else the default of its kind."""
    from src.layout_manager import load_layout, load_views
    path = views_path(workspace)
    migrate_views(path)
    fit_current(path, comparison_presets(file_manager, dataset_id))
    state = load_views(path, QUANT_COMPONENTS, comparison_presets(file_manager, dataset_id),
                       QUANT_QUANTIFICATION_LAYOUT)
    name = state['current'] if load_layout(path) is not None else default_view(file_manager, dataset_id)
    view = next((v for v in state['views'] if v['name'] == name), state['views'][0])
    return [list(r) for r in view['rows']]


VIEW_KEY = 'comparison_view'


def _options(workspace, file_manager, dataset_id):
    from src.render import missing_analysis

    def options(missing, view):
        del view
        return missing_analysis.comparison_options(
            workspace, file_manager, dataset_id, comparison_runs(file_manager, dataset_id), missing)
    return options


def view_selector(file_manager, workspace, dataset_id):
    """The view bar of a comparison: the views it can fill, and Identification
    locked (+) for runs that were not identified."""
    from src.layout_manager import view_bar
    from src.render.missing_analysis import ANALYSES
    migrate_views(views_path(workspace))
    fit_current(views_path(workspace), comparison_presets(file_manager, dataset_id))
    view_bar(
        views_path(workspace), QUANT_COMPONENTS, comparison_presets(file_manager, dataset_id),
        QUANT_QUANTIFICATION_LAYOUT,
        key=VIEW_KEY, default_view=default_view(file_manager, dataset_id),
        available=available(file_manager, dataset_id), analyses=ANALYSES,
        unavailable_options=_options(workspace, file_manager, dataset_id),
        dataset_label='comparison', target='comparisons of runs',
        help_text=(
            "Drag a view to order the bar, untick it to keep it under More. These views are used for "
            "every comparison of runs; your own take up to 3 components per row, ready-made views can "
            "be duplicated."
        ),
    )


@st.dialog('Edit comparison', width='large')
def _comparison_dialog(file_manager, dataset_id):
    from src.render import quant_cells  # needs openms_insight
    session = QuantSession(_own_file_manager(file_manager), dataset_id)
    st.caption(
        'Which runs (or labelled samples) belong to which condition, which level to compare and how. '
        'The defaults suit most experiments; **Apply** recomputes the comparison.'
    )
    if quant_cells.settings_form(session):
        st.rerun()


def buttons(file_manager, dataset_id):
    """Edit comparison, right of the experiment selector."""
    if st.button('Edit comparison', icon=':material/tune:', key='edit_comparison', width='content',
                 help='Conditions of the runs, level and analysis settings'):
        _comparison_dialog(file_manager, dataset_id)


def render(file_manager, dataset_id, rows=None):
    """The comparison grid of the Viewer."""
    # openms_insight only where the comparison is drawn, so the layout helpers
    # above import without it (the unit tests do not install it)
    from openms_insight import StateManager, defer_reruns
    from src.layout_manager import missing_analyses, missing_note, request_analysis
    from src.render import quant_cells
    from src.render.missing_analysis import ANALYSES
    has = available(file_manager, dataset_id)
    labels = {c['name']: c['label'] for c in QUANT_COMPONENTS}
    tool = comparison_tool(file_manager, dataset_id)
    session = QuantSession(file_manager, dataset_id)
    level = session.level
    if len(session.levels) == 1 and tool == 'flashtnt':
        st.caption('No PrSM identified a compared feature, so only features can be compared.')
    if session.isobaric:
        fractions = len(session.manifest.get('fractions', []))
        compared = (f"{len(session.runs)} {isobaric.label(session.manifest.get('plex', ''))} samples"
                    + (f' over {fractions} runs' if fractions > 1 else ''))
    else:
        compared = f'{len(session.runs)} runs'
    st.caption(f'Comparing **{LEVEL_NAMES[level].lower()}** of {compared}. '
               'Switch the level and set the conditions with **Edit comparison**.')
    # Selections are per comparison and level: an index means another entry on another level
    state_manager = StateManager(session_key=f'flash_quant_state_{tool}_{dataset_id}_{level}')
    with defer_reruns():
        for row_index, row in enumerate(rows if rows is not None else layout(file_manager)):
            columns = st.columns(len(row))
            for col_index, (column, comp_name) in enumerate(zip(columns, row)):
                with column:
                    missing = missing_analyses([comp_name], QUANT_COMPONENTS, has)
                    if missing:
                        missing_note(labels.get(comp_name, comp_name), missing, ANALYSES,
                                     lambda m, n: request_analysis(VIEW_KEY, m, component=n),
                                     f'quant_missing_{row_index}_{col_index}', comp_name,
                                     dataset_label='comparison')
                        continue
                    quant_cells.render_cell(
                        comp_name, session, state_manager, f'quant_{tool}_{level}_{row_index}_{col_index}_{comp_name}',
                    )
