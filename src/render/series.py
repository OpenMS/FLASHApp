"""Draw a small table as an openms-insight ``SeriesPlot`` (lines, bars or boxes).

An openms-insight without bars, boxes or views (0.7.5 and older) gets the same
plot drawn with Plotly, so the pages work on either.

The plots of a comparison's quality control, of its comparison steps and of
the GO enrichment are computed in the session from tables already in the
cache, so they are cached here by their content and settings: a rerun reuses
the plot, a new table or setting builds it again.
"""

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl

# Rows kept per box: evenly spaced values of the sorted ones, which draw the same box
BOX_POINTS = 201


def thin(frame: pl.DataFrame, by: list[str], value: str, points: int = BOX_POINTS) -> pl.DataFrame:
    """At most ``points`` evenly spaced values of ``value`` per group of ``by``."""
    frame = frame.filter(pl.col(value).is_not_null() & pl.col(value).is_finite())
    parts = []
    for _, part in frame.group_by(by, maintain_order=True):
        if part.height > points:
            rows = np.unique(np.round(np.linspace(0, part.height - 1, points)).astype(int))
            part = part.sort(value)[rows.tolist()]
        parts.append(part)
    return pl.concat(parts) if parts else frame


def _fingerprint(frame: pl.DataFrame, settings: dict) -> str:
    digest = hashlib.sha1(repr(sorted(settings.items())).encode())
    digest.update(repr((frame.columns, frame.shape)).encode())
    if frame.height:
        digest.update(str(int(frame.hash_rows(seed=0).sum())).encode())
    return digest.hexdigest()[:16]


def supports(plot_class, plot_kwargs: dict) -> bool:
    """Whether the installed ``SeriesPlot`` takes every argument and the mode.

    Bars, boxes, horizontal orientation and views came after openms-insight
    0.7.5; an install without them draws the plot with Plotly instead.
    """
    import inspect

    params = inspect.signature(plot_class.__init__).parameters
    if plot_kwargs.get('mode', 'lines+markers') not in LINE_MODES and 'orientation' not in params:
        return False
    return all(name in params for name in plot_kwargs)


LINE_MODES = ('lines', 'markers', 'lines+markers')


def show_series(data, name: str, cache_path, key: str, height=None, state_manager=None, **plot_kwargs):
    """Draw ``data`` (polars or pandas) as a ``SeriesPlot``; ``plot_kwargs`` go to it."""
    from openms_insight import SeriesPlot

    frame = pl.from_pandas(data) if isinstance(data, pd.DataFrame) else data
    if not supports(SeriesPlot, plot_kwargs):
        return _show_plotly(frame, key, height, state_manager, **plot_kwargs)
    cache_id = f'{name}-{_fingerprint(frame, plot_kwargs)}'
    plot = SeriesPlot(
        cache_id=cache_id, cache_path=str(Path(cache_path, 'insight', 'series')), data=frame.lazy(),
        **plot_kwargs,
    )
    plot(key=key, state_manager=state_manager, height=height or 400)
    return plot


def figure(frame: pl.DataFrame, x_column='x', y_column='y', series_column=None, series_names=None,
           series_colors=None, mode='lines+markers', orientation='v', hover_column=None, title=None,
           x_label=None, y_label=None, show_legend=True, y_range_mode='nonnegative', **_):
    """The Plotly figure of a ``SeriesPlot`` table (Streamlit-free)."""
    import plotly.graph_objects as go

    horizontal = orientation == 'h'
    fig = go.Figure()
    groups = frame.partition_by(series_column, maintain_order=True, as_dict=False) if series_column else [frame]
    for part in groups:
        label = str(part[series_column][0]) if series_column else None
        shown = (series_names or {}).get(label, label)
        color = (series_colors or {}).get(label)
        xs, ys = part[x_column].to_list(), part[y_column].to_list()
        axes = dict(x=ys, y=xs) if horizontal else dict(x=xs, y=ys)
        common = dict(name=shown, showlegend=bool(show_legend and series_column))
        if hover_column and hover_column in part.columns:
            common.update(hovertext=part[hover_column].to_list())
        if mode == 'bars':
            fig.add_trace(go.Bar(**axes, **common, orientation=orientation, marker_color=color))
        elif mode == 'box':
            fig.add_trace(go.Box(**axes, **common, orientation=orientation, marker_color=color))
        else:
            fig.add_trace(go.Scatter(**axes, **common, mode=mode, line_color=color))
    names = (y_label or y_column, x_label or x_column) if horizontal else (x_label or x_column, y_label or y_column)
    fig.update_layout(title=title, barmode='group', boxmode='group', showlegend=bool(show_legend and series_column),
                      xaxis_title=names[0], yaxis_title=names[1], margin=dict(t=50 if title else 20))
    value_axis = fig.update_xaxes if horizontal else fig.update_yaxes
    value_axis(rangemode=y_range_mode)
    if horizontal:
        fig.update_yaxes(autorange='reversed', type='category')
    return fig


def _show_plotly(frame, key, height, state_manager, view_column=None, view_names=None, view_identifier=None,
                 view_default=None, top_choices=None, top_default=None, bar_height=None, title=None, **plot_kwargs):
    """Draw a ``SeriesPlot`` table with Plotly, for an openms-insight without
    the plot's mode or views; a view and the number of rows are picked above it."""
    import streamlit as st

    from src.common.common import follow_app_theme

    if view_column and view_column in frame.columns:
        values = list(dict.fromkeys(frame[view_column].to_list()))
        order = [v for v in (view_names or {}) if v in values] + [v for v in values if v not in (view_names or {})]
        current = state_manager.get_selection(view_identifier) if state_manager and view_identifier else None
        start = current if current in order else (view_default if view_default in order else order[0])
        left, right = st.columns([3, 1]) if top_choices else (st.container(), None)
        view = left.segmented_control('View', order, default=start, key=f'{key}_view', label_visibility='collapsed',
                                      format_func=lambda v: (view_names or {}).get(v, v)) or start
        if state_manager and view_identifier and view != current:
            state_manager.set_selection(view_identifier, view)
        frame = frame.filter(pl.col(view_column) == view)
        title = (title or '').replace('{view}', (view_names or {}).get(view, str(view))) or None
    else:
        right = st.container() if top_choices else None
    if top_choices:
        default = top_default if top_default in top_choices else top_choices[len(top_choices) // 2]
        top = right.selectbox('Top', top_choices, index=top_choices.index(default), key=f'{key}_top',
                              format_func=lambda n: f'Top {n}', label_visibility='collapsed')
        keep = list(dict.fromkeys(frame[plot_kwargs.get('x_column', 'x')].to_list()))[:top]
        frame = frame.filter(pl.col(plot_kwargs.get('x_column', 'x')).is_in(keep))
    fig = figure(frame, title=title, **plot_kwargs)
    rows = frame[plot_kwargs.get('x_column', 'x')].n_unique() if frame.height else 0
    height = max(height or 400, (bar_height or 0) * rows + 120) if bar_height else (height or 400)
    fig.update_layout(height=height)
    st.plotly_chart(follow_app_theme(fig), key=f'{key}_plotly', use_container_width=True)
    return fig
