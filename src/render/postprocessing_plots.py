"""Diagnostic plots that help choose a setting of the run comparison.

Taken from OpenDDA (OpenMS/quantms-web, src/common/postprocessing_plots.py),
with proteins renamed to proteoforms. Every function takes plain pandas
tables and returns a Plotly figure, so the page only decides where to show
them (via ``show_fig``).
"""
from typing import Callable

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


def _log2(values: pd.DataFrame) -> pd.DataFrame:
    """log2 intensities with zeros and non-positive values treated as missing."""
    values = values.apply(pd.to_numeric, errors="coerce")
    return np.log2(values.where(values > 0))


def filter_threshold_curve(
    run_filter: Callable[[float], int],
    thresholds: list[float],
    current: float,
    total: int,
    x_label: str,
) -> go.Figure:
    """Proteins kept by a filter across its threshold range.

    ``run_filter(threshold)`` returns the number of proteins the filter keeps.
    A plateau means the threshold is not sensitive there; a steep drop shows
    where the filter starts removing many proteins.
    """
    kept = [run_filter(t) for t in thresholds]
    fig = go.Figure(go.Scatter(x=thresholds, y=kept, mode="lines+markers", name="Proteoforms kept"))
    fig.add_hline(y=total, line_dash="dot", line_color="grey", annotation_text=f"All proteoforms ({total})")
    fig.add_vline(x=current, line_dash="dash", line_color="#E74C3C", annotation_text="Current")
    fig.update_layout(
        title="Proteoforms kept at each threshold",
        xaxis_title=x_label,
        yaxis_title="Proteoforms kept",
        yaxis_rangemode="tozero",
        showlegend=False,
    )
    return fig


def missingness_vs_intensity(df: pd.DataFrame, sample_cols: list[str]) -> go.Figure:
    """Mean log2 intensity of proteoforms with and without missing values.

    If proteoforms with missing values sit clearly lower than complete ones,
    values are missing because they fall below the detection limit (MNAR);
    overlapping distributions point to random dropout (MAR).
    """
    log_values = _log2(df[sample_cols])
    has_missing = log_values.isna().any(axis=1)
    plot_df = pd.DataFrame(
        {
            "Mean log2 intensity": log_values.mean(axis=1),
            "Proteoforms": np.where(has_missing, "With missing values", "Complete"),
        }
    ).dropna()
    fig = px.histogram(
        plot_df,
        x="Mean log2 intensity",
        color="Proteoforms",
        barmode="overlay",
        histnorm="probability density",
        opacity=0.6,
        nbins=40,
        color_discrete_map={"Complete": "#3498DB", "With missing values": "#E67E22"},
    )
    fig.update_layout(
        title=f"Where missing values occur ({int(has_missing.sum())} of {len(df)} proteoforms have gaps)",
        yaxis_title="Density",
    )
    return fig


def imputed_value_preview(
    before: pd.DataFrame, after: pd.DataFrame, sample_cols: list[str]
) -> go.Figure:
    """log2 distribution of observed values versus the values imputation adds."""
    before_log = _log2(before[sample_cols])
    after_log = _log2(after[sample_cols])
    missing = before_log.isna()
    observed = before_log.stack().rename("log2 intensity").to_frame()
    observed["Values"] = "Observed"
    imputed = after_log.where(missing).stack().rename("log2 intensity").to_frame()
    imputed["Values"] = "Imputed"
    plot_df = pd.concat([observed, imputed], ignore_index=True)
    fig = px.histogram(
        plot_df,
        x="log2 intensity",
        color="Values",
        barmode="overlay",
        opacity=0.6,
        nbins=60,
        color_discrete_map={"Observed": "#3498DB", "Imputed": "#E74C3C"},
    )
    fig.update_layout(title="Observed values and the values this setting would fill in", yaxis_title="Count")
    return fig


def sample_distributions(
    before: pd.DataFrame,
    after: pd.DataFrame,
    sample_cols: list[str],
    sample_group_map: dict,
    after_is_log: bool,
) -> go.Figure:
    """Per-sample box plots of log2 intensities before and after normalization.

    Well-normalized samples have aligned medians; a sample that stays offset
    after normalization is worth checking. Untransformed tables are shown on
    a log2 scale so both panels are readable.
    """
    frames = []
    for label, table, is_log in [("Before", before, False), ("After", after, after_is_log)]:
        values = table[sample_cols].apply(pd.to_numeric, errors="coerce")
        if not is_log:
            values = _log2(values)
        long_df = values.melt(var_name="Sample", value_name="Intensity").dropna()
        long_df["Stage"] = label
        frames.append(long_df)
    plot_df = pd.concat(frames, ignore_index=True)
    plot_df["Group"] = plot_df["Sample"].map(sample_group_map).fillna("Unassigned")
    fig = px.box(
        plot_df,
        x="Sample",
        y="Intensity",
        color="Group",
        facet_row="Stage",
        category_orders={"Stage": ["Before", "After"], "Sample": sample_cols},
        color_discrete_sequence=px.colors.qualitative.Set2,
        points=False,
    )
    fig.for_each_annotation(lambda a: a.update(text=a.text.replace("Stage=", "")))
    fig.update_yaxes(title_text="log intensity", matches=None)
    fig.update_layout(title="Sample intensity distributions before and after", height=650)
    return fig


def pvalue_histogram(statistics_df: pd.DataFrame) -> go.Figure:
    """Histogram of raw p-values.

    A flat histogram with a peak near 0 is what a sound test looks like. A
    U-shape or a peak near 1 suggests the input scale or the chosen test does
    not fit the data.
    """
    pvalues = pd.to_numeric(statistics_df["p-value"], errors="coerce").dropna()
    fig = px.histogram(x=pvalues, nbins=20, range_x=[0, 1])
    fig.update_traces(xbins=dict(start=0, end=1, size=0.05), marker_color="#3498DB")
    expected = len(pvalues) / 20
    fig.add_hline(y=expected, line_dash="dash", line_color="grey", annotation_text="Expected if nothing changes")
    fig.update_layout(title="p-value distribution", xaxis_title="p-value", yaxis_title="Proteoforms")
    return fig
