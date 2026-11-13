"""Live overview of what FLASHIda decided in a run (FLASHIda Viewer page).

Reads the run's logs (src/parse/ida_log.py) and shows them with openms-insight:
a map of the survey masses with the precursors FLASHIda chose, the precursors,
and the scans it commanded for the selected precursor. The components are
rebuilt whenever the logs grow; the selection is kept.
"""

import hashlib
import inspect
import shutil
import time
from pathlib import Path

import polars as pl
import streamlit as st

from openms_insight import Heatmap, StateManager, Table, defer_reruns
from openms_insight.core.cache import CacheMissError

from src.parse import ida_log
from src.render.cache_dirs import build_lock, short_name

# Code of the frame builders and of this module: part of every cache id
CODE_VERSION = hashlib.md5(
    Path(__file__).read_bytes() + Path(ida_log.__file__).read_bytes()
).hexdigest()[:8]

COLORS = {'Deconvolved': '#9aa5b1', 'Chosen': '#f08c00', 'Identified': '#2b8a3e'}

# What the map can be colored by besides the intensity (toolbar "Color"), and
# its layers (toolbar "Layers"). Needs an openms-insight whose Heatmap takes
# `color_columns`; an older one keeps the map colored by decision.
COLOR_COLUMNS = {
    'QScore': 'qscore',
    '# deconvolved MS2 masses': 'ms2_masses',
    '# fragment ions': 'fragment_ions',
}
HAS_COLOR_OPTIONS = 'color_columns' in inspect.signature(Heatmap.__init__).parameters


def _fixed(precision=4):
    return {'formatter': 'fixed', 'formatterParams': {'precision': precision, 'minLength': 4}}


def _tick():
    return {'formatter': 'tickCross', 'hozAlign': 'center'}


PRECURSOR_COLUMNS = [
    {'title': '#', 'field': 'precursor', 'sorter': 'number',
     'headerTooltip': 'Precursors in the order FLASHIda chose them.'},
    {'title': 'RT', 'field': 'RT', 'sorter': 'number', **_fixed(2),
     'headerTooltip': 'Retention time (min) of the survey scan it was chosen from.'},
    {'title': 'Mass', 'field': 'Mass', 'sorter': 'number', **_fixed(),
     'headerTooltip': 'Monoisotopic mass.'},
    {'title': 'z', 'field': 'z', 'sorter': 'number'},
    {'title': 'QScore', 'field': 'QScore', 'sorter': 'number', **_fixed(3)},
    {'title': 'MS2/MS3', 'field': 'MS2/MS3',
     'headerTooltip': 'MS2 scans commanded for it / MS3 scans for its fragments.'},
    {'title': 'Scans for', 'field': 'Decisions',
     'headerTooltip': 'What its scans were for: recording, exploration, quantification, ...'},
    {'title': 'Best CE', 'field': 'Best energy', 'sorter': 'number',
     'headerTooltip': 'Energy of the winning variant of an exploration.'},
    {'title': 'Quant', 'field': 'Quant verdict',
     'headerTooltip': 'Quantification verdict.'},
    {'title': 'Identified', 'field': 'Identified', **_tick()},
    {'title': 'Score', 'field': 'Score', 'sorter': 'number', **_fixed(0)},
    {'title': 'Cov. %', 'field': 'Coverage', 'sorter': 'number', **_fixed(1),
     'headerTooltip': 'Sequence coverage of the fragments pooled over its scans.'},
    {'title': 'Proteoform', 'field': 'Proteoform', 'hozAlign': 'left'},
]

SCAN_COLUMNS = [
    {'title': 'Scan', 'field': 'Scan',
     'headerTooltip': "FLASHIda's id of the commanded scan."},
    {'title': 'Status', 'field': 'Status'},
    {'title': 'Decision', 'field': 'Decision'},
    {'title': 'MS', 'field': 'MS', 'sorter': 'number'},
    {'title': 'RT', 'field': 'RT', 'sorter': 'number', **_fixed(2),
     'headerTooltip': 'Retention time (min) of the acquired scan.'},
    {'title': 'Fragment', 'field': 'Fragment',
     'headerTooltip': 'Fragment ion isolated for an MS3 scan.'},
    {'title': 'Mass', 'field': 'Mass', 'sorter': 'number', **_fixed()},
    {'title': 'z', 'field': 'z', 'sorter': 'number'},
    {'title': 'm/z', 'field': 'm/z', 'sorter': 'number', **_fixed()},
    {'title': 'Width', 'field': 'Width', 'sorter': 'number', **_fixed(2),
     'headerTooltip': 'Isolation width (Th).'},
    {'title': 'Activation', 'field': 'Activation'},
    {'title': 'Energy', 'field': 'Energy', 'sorter': 'number'},
    {'title': 'Variant', 'field': 'Variant',
     'headerTooltip': 'Variant of an exploration (e.g. a collision energy step).'},
    {'title': 'Expl. score', 'field': 'Exploration score', 'sorter': 'number', **_fixed(2)},
    {'title': 'Winner', 'field': 'Winner', **_tick(),
     'headerTooltip': 'The variant the exploration chose.'},
    {'title': 'Masses', 'field': 'Masses', 'sorter': 'number',
     'headerTooltip': 'Masses deconvolved from the scan.'},
    {'title': 'Tags', 'field': 'Tags', 'sorter': 'number'},
    {'title': 'Quant', 'field': 'Quant verdict'},
    {'title': 'Score', 'field': 'Score', 'sorter': 'number', **_fixed(0)},
    {'title': 'Proteoform', 'field': 'Proteoform', 'hozAlign': 'left'},
    {'title': 'QScore', 'field': 'QScore', 'sorter': 'number', **_fixed(3)},
    {'title': 'Priority', 'field': 'Priority', 'sorter': 'number'},
    {'title': 'Parent', 'field': 'Parent',
     'headerTooltip': 'The scan this one was decided from.'},
]


def _no_collapse(columns):
    """Columns that never fold into a row's details: a table too narrow for
    them scrolls instead."""
    return [{**c, 'responsive': 0} for c in columns]


def _component(cls, name, cache_path, stamp, build, **styling):
    """Insight component `name` from its cache, or built when the logs changed.
    `styling` (colors, hidden columns) is not stored with the cache, so it is
    passed either way."""
    cache_id = f'{name}-{stamp}-{CODE_VERSION}'
    try:
        return cls(cache_id=cache_id, cache_path=str(cache_path), **styling)
    except CacheMissError:
        pass
    with build_lock(cache_path, cache_id):
        # Another tab may have built it while this one waited
        try:
            return cls(cache_id=cache_id, cache_path=str(cache_path), **styling)
        except CacheMissError:
            pass
        if cache_path.exists():
            for old in cache_path.glob(f'{name}-*'):
                if old.name != cache_id:
                    shutil.rmtree(old, ignore_errors=True)
        return cls(cache_id=cache_id, cache_path=str(cache_path), **build(), **styling)


# Streamlit greys out every element of a fragment while it reruns, until the
# element is sent again. A live view reruns every few seconds and takes a
# second or two to read the logs, so its plots and tables were grey half the
# time; inside this container they keep showing the last data until the new
# data replaces it.
NO_STALE_CSS = (
    '<style>.st-key-{key} [data-stale="true"]'
    '{{opacity:1 !important;transition:none !important}}</style>'
)


def live_container(key: str):
    """Container for a live view whose elements are not greyed out while it
    refreshes."""
    st.html(NO_STALE_CSS.format(key=key))
    return st.container(key=key)


_DECISIONS = {}


def run_decisions(run: Path, stamp: str) -> dict:
    """What FLASHIda decided in `run`, read again only when its logs changed
    (`stamp`): a refresh while FLASHIda writes nothing then costs no parsing."""
    cached = _DECISIONS.get(str(run))
    if cached is not None and cached[0] == stamp:
        return cached[1]
    result = ida_log.decisions(ida_log.read_run(run))
    # One run per Viewer is shown at a time: keep a few, not every run opened
    if len(_DECISIONS) >= 4 and str(run) not in _DECISIONS:
        _DECISIONS.pop(next(iter(_DECISIONS)))
    _DECISIONS[str(run)] = (stamp, result)
    return result


def _frame(df):
    if 'MS3' in df.columns and 'Waiting' in df.columns:
        df = df.assign(**{'MS2/MS3': [f'{a} / {b}' for a, b in zip(df['MS2'], df['MS3'])]})
    return pl.from_pandas(df).lazy()


def _count(n):
    return f'{n:,}'


def summary_metrics(summary, last_change):
    # One tile per count: a combined "MS2 / MS3" value no longer fit its
    # column once both reached four digits
    cols = st.columns(7)
    cols[0].metric('Survey scans', _count(summary['surveys']))
    cols[1].metric('Precursors', _count(summary['precursors']), help='Precursors FLASHIda chose.')
    cols[2].metric('MS2 scans', _count(summary['ms2']), help='MS2 scans FLASHIda commanded.')
    cols[3].metric('MS3 scans', _count(summary['ms3']), help='MS3 scans FLASHIda commanded.')
    missed = summary.get('not_acquired', 0)
    cols[4].metric(
        'Waiting', _count(summary['waiting']),
        help='Commanded scans the instrument has not returned yet.'
        + (f' {missed:,} earlier commanded scans never came back.' if missed else ''))
    cols[5].metric('Identified', _count(summary['identified']), help='Precursors with a proteoform.')
    cols[6].metric('RT (min)', '-' if summary['rt'] is None else f"{summary['rt']:.2f}")
    if last_change is not None:
        age = time.time() - last_change
        if age < 60:
            st.caption(f'🟢 Last log entry {age:.0f} s ago.')
        else:
            st.caption(f"⚪ No new log entries since {time.strftime('%H:%M:%S', time.localtime(last_change))}.")


def show_run(run: Path, cache_root: Path, key: str):
    """Overview of a run's decisions; rebuilt when its logs grow."""
    stamp = ida_log.stamp(run)
    result = run_decisions(run, stamp)
    summary_metrics(result['summary'], ida_log.last_change(run))
    if not len(result['points']) and not len(result['precursors']):
        st.info('FLASHIda has not logged any scans for this run yet.')
        return

    # Run names are long raw file names: a short folder keeps the cache
    # files within Windows' path length limit
    cache_path = Path(cache_root, short_name(run.name))
    state_manager = StateManager(session_key=f'flashida_live_state_{key}')
    run_key = f'flashida_live_run_{key}'
    if st.session_state.get(run_key) != str(run):
        state_manager.clear()
        st.session_state[run_key] = str(run)

    points = result['points']
    common = dict(
        x_column='rt', y_column='mass', intensity_column='intensity',
        interactivity={'precursor': 'precursor'}, highlight={'precursor': 'precursor'},
        hover_column='label', zoom_identifier=f'flashida_live_zoom_{key}',
        title='Survey masses and chosen precursors', x_label='Retention time (min)',
        y_label='Monoisotopic mass (Da)', min_points=50000, downsample_zoomed=False,
        svg_filename='FLASHIda-decisions', scroll_zoom=False,
    )
    if HAS_COLOR_OPTIONS:
        # Colored by a value, chosen in the toolbar; survey masses, chosen and
        # identified precursors are layers to show or hide
        heatmap = _component(
            Heatmap, 'decision_map_values', cache_path, stamp,
            lambda: dict(data=_frame(points), **common),
            # Render-time options, not stored with the cache
            color_columns=COLOR_COLUMNS, layer_column='Layer',
        )
    else:
        heatmap = _component(
            Heatmap, 'decision_map', cache_path, stamp,
            lambda: dict(data=_frame(points), category_column='Decision', **common),
            # Render-time styling, not stored with the cache
            category_colors=COLORS,
        )
    precursors = _component(
        Table, 'precursors', cache_path, stamp, lambda: dict(
            data=_frame(result['precursors']), title='Chosen precursors',
            pagination_identifier=f'flashida_live_precursors_page_{key}',
            index_field='precursor', go_to_fields=['precursor', 'Mass'],
            interactivity={'precursor': 'precursor'}, column_definitions=_no_collapse(PRECURSOR_COLUMNS),
            initial_sort=[{'column': 'precursor', 'dir': 'desc'}],
            reclick='keep', scroll_position='top',
            toggle_filters=[{'label': 'Identified', 'column': 'Identified', 'default': False}],
        ),
        hidden_columns=['Identified'],
    )
    scans = _component(
        Table, 'scans', cache_path, stamp, lambda: dict(
            data=_frame(result['scans']), title='Scans commanded for the precursor',
            pagination_identifier=f'flashida_live_scans_page_{key}',
            index_field='scan', filters={'precursor': 'precursor'},
            column_definitions=_no_collapse(SCAN_COLUMNS), reclick='keep', scroll_position='top',
            default_row_selects=False,
        ),
        hidden_columns=['m/z', 'Width', 'Exploration score', 'Tags', 'Quant verdict', 'QScore',
                        'Priority', 'Parent'],
    )
    with defer_reruns():
        heatmap(key=f'{key}_map', state_manager=state_manager, height=420)
        precursors(key=f'{key}_precursors', state_manager=state_manager, height=360)
        scans(key=f'{key}_scans', state_manager=state_manager, height=320)
