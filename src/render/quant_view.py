"""Comparisons of runs in the Viewer.

A workflow run that quantified two or more files stores a comparison
(src/quant/session.py) next to its runs. The Viewer lists it in the experiment
selector as "Comparison: ..."; selecting it shows a grid of comparison
components (src/render/quant_cells.py). A view bar of its own picks the grid:
the ready-made views **Identification**, **Quantification**, **QC: label-free**
and **QC: TMT, iTRAQ**, and the user's own (+ and ✏️), in
`<workspace>/layouts/comparison.json`. Until one is picked, a comparison shows
Quantification. The layout saved with the
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
    QUANT_COMPONENTS, QUANT_DEFAULT_LAYOUT, QUANT_DIAGNOSTICS, QUANT_IDENTIFICATION_LAYOUT,
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
    """The dataset ids of the runs a comparison compares."""
    manifest = file_manager.get_results(dataset_id, ['quant_manifest'])['quant_manifest']
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


def presets():
    """The ready-made comparison views: what the runs identified, what
    changed between the conditions, and the quality control of a label-free
    and of an isobaric (TMT, iTRAQ) comparison."""
    return [
        {'name': IDENTIFICATION_VIEW, 'label': 'Identification', 'layout': QUANT_IDENTIFICATION_LAYOUT,
         'help': 'Every proteoform with its status run by run; pick one to list its PrSMs in every run, and a '
                 'PrSM to see its sequence view and spectrum. Then what each run identified and what the runs '
                 'share. Needs no conditions.'},
        {'name': QUANTIFICATION_VIEW, 'label': 'Quantification', 'layout': QUANT_QUANTIFICATION_LAYOUT,
         'help': 'What changed between the conditions: volcano plot, comparison table, each entry run by run, '
                 'sample similarity, the clustered heatmap and GO enrichment.'},
        {'name': LABEL_FREE_QC_VIEW, 'label': 'QC: label-free', 'layout': QUANT_LABEL_FREE_QC_LAYOUT,
         'help': 'Quality control of runs measured one by one: the consensus features (feature map, intensity '
                 'of detected and matched values, every feature run by run), match between runs, how the runs '
                 'were loaded, how well replicates agree, missing values and p-values.'},
        {'name': ISOBARIC_QC_VIEW, 'label': 'QC: TMT, iTRAQ', 'layout': QUANT_ISOBARIC_QC_LAYOUT,
         'help': 'Quality control of the labelled samples of one set: the reporter ions (intensity and ratio '
                 'per label, labels found per spectrum, every labelled spectrum and its ions), whether the '
                 'channels were loaded alike, how they correlate and cluster, missing values and p-values.'},
    ]


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
    data = {'views': [{'name': p['name'], 'in_bar': True} for p in presets()]
            + [{'name': name, 'label': label, 'rows': legacy, 'in_bar': True}]}
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=1))


def view_rows(file_manager, workspace, dataset_id):
    """The rows a comparison shows: the picked view, else the default of its kind."""
    from src.layout_manager import load_layout, load_views
    path = views_path(workspace)
    migrate_views(path)
    state = load_views(path, QUANT_COMPONENTS, presets(), QUANT_IDENTIFICATION_LAYOUT)
    name = state['current'] if load_layout(path) is not None else default_view(file_manager, dataset_id)
    view = next((v for v in state['views'] if v['name'] == name), state['views'][0])
    return [list(r) for r in view['rows']]


def view_selector(file_manager, workspace, dataset_id):
    """The view bar of a comparison."""
    from src.layout_manager import view_bar
    migrate_views(views_path(workspace))
    view_bar(
        views_path(workspace), QUANT_COMPONENTS, presets(), QUANT_IDENTIFICATION_LAYOUT,
        key='comparison_view', default_view=default_view(file_manager, dataset_id),
        target='comparisons of runs',
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
    from src.render import quant_cells
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
                    quant_cells.render_cell(
                        comp_name, session, state_manager, f'quant_{tool}_{level}_{row_index}_{col_index}_{comp_name}',
                    )
