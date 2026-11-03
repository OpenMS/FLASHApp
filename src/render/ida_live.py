"""Live overview of what FLASHIda decided in a run (FLASHIda Viewer page).

Reads the run's logs (src/parse/ida_log.py) and shows them with openms-insight:
a map of the survey masses with the precursors FLASHIda chose, the precursors,
and the scans it commanded for the selected precursor. The components are
rebuilt whenever the logs grow; the selection is kept.
"""

import hashlib
import shutil
import time
from pathlib import Path

import polars as pl
import streamlit as st

from openms_insight import Heatmap, StateManager, Table, defer_reruns
from openms_insight.core.cache import CacheMissError

from src.parse import ida_log

# Code of the frame builders and of this module: part of every cache id
CODE_VERSION = hashlib.md5(
    Path(__file__).read_bytes() + Path(ida_log.__file__).read_bytes()
).hexdigest()[:8]

COLORS = {'Deconvolved': '#9aa5b1', 'Chosen': '#f08c00', 'Identified': '#2b8a3e'}


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
    if cache_path.exists():
        for old in cache_path.glob(f'{name}-*'):
            if old.name != cache_id:
                shutil.rmtree(old, ignore_errors=True)
    return cls(cache_id=cache_id, cache_path=str(cache_path), **build(), **styling)


def _frame(df):
    if 'MS3' in df.columns and 'Waiting' in df.columns:
        df = df.assign(**{'MS2/MS3': [f'{a} / {b}' for a, b in zip(df['MS2'], df['MS3'])]})
    return pl.from_pandas(df).lazy()


def summary_metrics(summary, last_change):
    cols = st.columns(6)
    cols[0].metric('Survey scans', summary['surveys'])
    cols[1].metric('Precursors chosen', summary['precursors'])
    cols[2].metric('MS2 / MS3 commanded', f"{summary['ms2']} / {summary['ms3']}")
    cols[3].metric('Waiting', summary['waiting'], help='Commanded scans not acquired yet.')
    cols[4].metric('Identified', summary['identified'], help='Precursors with a proteoform.')
    cols[5].metric('RT', '-' if summary['rt'] is None else f"{summary['rt']:.2f} min")
    if last_change is not None:
        age = time.time() - last_change
        if age < 60:
            st.caption(f'🟢 Last log entry {age:.0f} s ago.')
        else:
            st.caption(f"⚪ No new log entries since {time.strftime('%H:%M:%S', time.localtime(last_change))}.")


def show_run(run: Path, cache_root: Path, key: str):
    """Overview of a run's decisions; rebuilt when its logs grow."""
    stamp = ida_log.stamp(run)
    result = ida_log.decisions(ida_log.read_run(run))
    summary_metrics(result['summary'], ida_log.last_change(run))
    if not len(result['points']) and not len(result['precursors']):
        st.info('FLASHIda has not logged any scans for this run yet.')
        return

    cache_path = Path(cache_root, run.name)
    state_manager = StateManager(session_key=f'flashida_live_state_{key}')
    run_key = f'flashida_live_run_{key}'
    if st.session_state.get(run_key) != str(run):
        state_manager.clear()
        st.session_state[run_key] = str(run)

    points = result['points']
    heatmap = _component(
        Heatmap, 'decision_map', cache_path, stamp, lambda: dict(
            data=_frame(points), x_column='rt', y_column='mass', intensity_column='intensity',
            category_column='Decision',
            interactivity={'precursor': 'precursor'}, highlight={'precursor': 'precursor'},
            hover_column='label', zoom_identifier=f'flashida_live_zoom_{key}',
            title='Survey masses and chosen precursors', x_label='Retention time (min)',
            y_label='Monoisotopic mass (Da)', min_points=50000, downsample_zoomed=False,
            svg_filename='FLASHIda-decisions', scroll_zoom=False,
        ),
        # Render-time styling, not stored with the cache
        category_colors=COLORS,
    )
    precursors = _component(
        Table, 'precursors', cache_path, stamp, lambda: dict(
            data=_frame(result['precursors']), title='Chosen precursors',
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
