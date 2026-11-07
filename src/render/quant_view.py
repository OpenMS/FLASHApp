"""Comparisons of runs in the Viewer.

A workflow run that quantified two or more files stores a comparison
(src/quant/session.py) next to its runs. The Viewer lists it in the experiment
selector as "Comparison: ..."; selecting it shows a grid of Quantification
components (src/render/quant_cells.py). A view bar of its own picks the grid:
the ready-made views **Label-free** and **TMT, iTRAQ** (isobaric labels), and
the user's own (+ and ✏️), in `<workspace>/layouts/comparison.json`. Until one
is picked, a comparison shows the view of its kind. The layout saved with the
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
    QUANT_COMPONENTS, QUANT_DEFAULT_LAYOUT, QUANT_ISOBARIC_LAYOUT, QUANT_LABEL_FREE_LAYOUT, QUANT_NAMES,
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


def layout(file_manager):
    """The layout saved with the former Edit layout dialog, else the default."""
    saved = get_saved_layout(file_manager, LAYOUT_DATASET)
    return saved[0][0] if saved else QUANT_DEFAULT_LAYOUT


LABEL_FREE_VIEW = 'label-free'
ISOBARIC_VIEW = 'isobaric'
# The layout saved with the former Edit layout dialog, as a view of one's own
LEGACY_VIEW = ('custom-comparison', 'My comparison layout')


def views_path(workspace):
    """Where comparison views were kept before the Viewer listed every view in
    one bar (and where the former layout is adopted, then merged)."""
    return Path(workspace, 'layouts', 'comparison.json')


def presets():
    """The ready-made comparison views, by how the samples were quantified."""
    return [
        {'name': LABEL_FREE_VIEW, 'label': 'Label-free', 'layout': QUANT_LABEL_FREE_LAYOUT,
         'help': 'Runs compared feature by feature: volcano plot, comparison table, each entry run by run, '
                 'sample similarity and GO enrichment.'},
        {'name': ISOBARIC_VIEW, 'label': 'TMT, iTRAQ', 'layout': QUANT_ISOBARIC_LAYOUT,
         'help': 'The labelled samples of one set (TMT, iTRAQ): volcano plot, comparison table, how the '
                 'samples cluster and the diagnostics of each step, then GO enrichment.'},
    ]


def is_isobaric(file_manager, dataset_id):
    """True for a comparison of labelled samples (TMT, iTRAQ)."""
    manifest = file_manager.get_results(dataset_id, ['quant_manifest'])['quant_manifest']
    return manifest.get('kind') == 'isobaric'


def _legacy_rows(file_manager):
    saved = get_saved_layout(file_manager, LAYOUT_DATASET)
    return saved[0][0] if saved else None


def default_view(file_manager, dataset_id=None):
    """The view a comparison shows until one is picked: the layout saved
    before the view bar, else the ready-made view of its kind."""
    from src.layout_manager import match_preset
    legacy = _legacy_rows(file_manager)
    if legacy:
        return match_preset([legacy], presets(), QUANT_COMPONENTS) or LEGACY_VIEW[0]
    if dataset_id is not None and is_isobaric(file_manager, dataset_id):
        return ISOBARIC_VIEW
    return LABEL_FREE_VIEW


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


def merge_views(workspace):
    """Move the comparison views of the former bar into the Viewer's views, once.

    The user's own comparison views are added (renamed where a name is taken),
    ready-made ones keep whether they were in the bar. The file is kept as
    ``comparison.json.merged``.
    """
    import json
    from src.render import topdown_layout
    source = views_path(workspace)
    if not source.exists():
        return
    try:
        old = json.loads(source.read_text())
    except (OSError, ValueError):
        old = {}
    target = topdown_layout.views_path(workspace)
    try:
        data = json.loads(target.read_text()) if target.exists() else {}
    except (OSError, ValueError):
        data = {}
    if not isinstance(data, dict):
        data = {}
    views = [v for v in data.get('views') or [] if isinstance(v, dict) and isinstance(v.get('name'), str)]
    if not views:
        # The views of a run first, as the bar lists them
        views = [{'name': name, 'in_bar': True}
                 for name in [topdown_layout.DECONVOLUTION_VIEW, *topdown_layout.SAMPLE_VIEWS.values()]]
    names = {v['name'] for v in views}
    ready = {p['name'] for p in presets()}
    for entry in (old.get('views') if isinstance(old, dict) else None) or []:
        if not isinstance(entry, dict) or not isinstance(entry.get('name'), str):
            continue
        entry = dict(entry)
        if entry['name'] in ready:
            if entry['name'] not in names:
                views.append({'name': entry['name'], 'in_bar': bool(entry.get('in_bar', True))})
                names.add(entry['name'])
            continue
        if not entry['name'].startswith('custom-') or not isinstance(entry.get('rows'), list):
            continue
        name, n = entry['name'], 1
        while name in names:
            n += 1
            name = f"custom-comparison-{n}" if n > 1 else 'custom-comparison'
        entry['name'] = name
        views.append(entry)
        names.add(name)
    if views:
        data['views'] = views
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(data, indent=1))
    source.replace(source.with_name(source.name + '.merged'))


def view_rows(file_manager, workspace, dataset_id):
    """The rows a comparison shows: the picked view, else the default of its kind."""
    from src.layout_manager import load_layout, load_views
    from src.render import topdown_layout
    path = topdown_layout.views_path(workspace)
    state = load_views(path, topdown_layout.components(file_manager), topdown_layout.presets(file_manager),
                       topdown_layout.deconvolution_default(file_manager))
    name = state['current'] if load_layout(path) is not None else default_view(file_manager, dataset_id)
    view = next((v for v in state['views'] if v['name'] == name),
                next(v for v in state['views'] if v['name'] == LABEL_FREE_VIEW))
    return [list(r) for r in view['rows']]


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


def render(file_manager, dataset_id, rows=None, on_missing=None):
    """The comparison grid of the Viewer. A component of a run in ``rows``
    draws a note instead, whose button calls ``on_missing(missing, component)``."""
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
                    if comp_name not in QUANT_NAMES:
                        from src.layout_manager import missing_note
                        from src.render.layout_editor.catalog import RUN, TOPDOWN_COMPONENTS
                        from src.render.missing_analysis import ANALYSES
                        labels = {c['name']: c['label'] for c in TOPDOWN_COMPONENTS}
                        missing_note(labels.get(comp_name, comp_name), [RUN], ANALYSES, on_missing,
                                     f'quant_{tool}_{level}_{row_index}_{col_index}_missing', comp_name,
                                     dataset_label='comparison')
                        continue
                    quant_cells.render_cell(
                        comp_name, session, state_manager, f'quant_{tool}_{level}_{row_index}_{col_index}',
                    )
