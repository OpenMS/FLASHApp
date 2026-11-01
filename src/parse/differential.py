"""Group comparison in one go: filtering, imputation, normalization and statistics.

Taken from OpenDDA (OpenMS/quantms-web, src/common/differential.py) so that
FLASHQuant's run comparison runs the same chain of openms_insight functions,
with the same settings and defaults, on proteoforms instead of proteins. The
differences: the defaults live in DEFAULT_SETTINGS instead of
default-parameters.json, the messages speak of proteoforms, and there are no
step pages to continue from (``done`` is kept for parity).
"""

import pandas as pd
import polars as pl
from openms_insight.analysis.filter import (
    filter_low_abundance,
    filter_low_repeatability,
    filter_low_variance,
)
from openms_insight.analysis.imputation import impute_mar, impute_smallest_value
from openms_insight.analysis.normalization import normalize_samples, scale_data, transform_data
from openms_insight.analysis.statistics import adjust_fdr_lazy, calculate_statistical_tests

STAT_COLUMNS = ["log2FC", "p-value", "p-adj", "stat"]

# Session-state keys of the chain, in order (same as results_helpers.POSTPROCESSING_STEPS).
STEPS = ["filtered_df", "imputed_df", "normalized_df", "statistics_df"]

SETTING_KEYS = [
    "postproc-filter-method",
    "postproc-filter-abundance-percentile",
    "postproc-filter-max-missing",
    "postproc-filter-variance-percentile",
    "postproc-impute-class",
    "postproc-impute-mar-strategy",
    "postproc-impute-mnar-scope",
    "postproc-transform",
    "postproc-normalization",
    "postproc-scaling",
    "postproc-stat-test",
    "postproc-fdr",
]


# OpenDDA's shipped defaults (its default-parameters.json)
DEFAULT_SETTINGS = {
    "postproc-filter-method": "Low Repeatability",
    "postproc-filter-abundance-percentile": 10.0,
    "postproc-filter-max-missing": 50.0,
    "postproc-filter-variance-percentile": 10.0,
    "postproc-impute-class": "MNAR (Missing Not At Random)",
    "postproc-impute-mar-strategy": "median",
    "postproc-impute-mnar-scope": "row",
    "postproc-transform": "log2",
    "postproc-normalization": "median",
    "postproc-scaling": "None",
    "postproc-stat-test": "limma_like",
    "postproc-fdr": "BH",
}

# Choices offered for each setting (as on OpenDDA's Fine-tune pages)
FILTER_METHODS = ["Low Repeatability", "Low Abundance", "Low Variance"]
IMPUTE_CLASSES = ["MNAR (Missing Not At Random)", "MAR (Missing At Random)"]
MAR_STRATEGIES = ["median", "mean"]
MNAR_SCOPES = ["row", "global"]
TRANSFORMS = ["None", "log2", "log10", "square_root", "cube_root"]
NORMALIZATIONS = ["None", "sum", "median", "pqn", "quantile"]
SCALINGS = ["None", "mean_centering", "auto_scaling", "pareto_scaling", "range_scaling"]
FDR_METHODS = ["BH", "Bonferroni", "None"]


def test_options(n_groups: int) -> list[str]:
    return ["limma_like", "welch", "paired"] if n_groups == 2 else ["limma_like", "anova"]


def default_settings() -> dict:
    """The shipped downstream-analysis defaults."""
    return dict(DEFAULT_SETTINGS)


def resolve_settings(params: dict, defaults: dict | None = None) -> dict:
    """Saved settings where present, shipped defaults for the rest."""
    defaults = default_settings() if defaults is None else defaults
    return {k: params.get(k, defaults[k]) for k in SETTING_KEYS}


def strip_stat_columns(df: pd.DataFrame) -> pd.DataFrame:
    return df.drop(columns=[c for c in STAT_COLUMNS if c in df.columns], errors="ignore")


def valid_test(test: str, n_groups: int) -> str:
    """The saved test if it fits the number of groups, else limma_like."""
    return test if test in test_options(n_groups) else "limma_like"


def describe_settings(settings: dict) -> list[str]:
    """One plain-language line per step, for showing what was run."""
    method = settings["postproc-filter-method"]
    if method == "Low Repeatability":
        filt = (
            f"kept proteoforms measured in at least {100 - float(settings['postproc-filter-max-missing']):g} % "
            "of the samples of at least one group"
        )
    elif method == "Low Abundance":
        filt = f"removed the {float(settings['postproc-filter-abundance-percentile']):g} % least abundant proteoforms"
    else:
        filt = f"removed the {float(settings['postproc-filter-variance-percentile']):g} % least variable proteoforms"
    if settings["postproc-impute-class"].startswith("MNAR"):
        imp = (
            "filled missing values with the proteoform's lowest measured value (assumes they were too faint to detect)"
            if settings["postproc-impute-mnar-scope"] == "row"
            else "filled missing values with the lowest value of the whole dataset"
        )
    else:
        imp = f"filled missing values with the group {settings['postproc-impute-mar-strategy']}"
    norm = f"{settings['postproc-transform']} transform, {settings['postproc-normalization']} normalization"
    if settings["postproc-scaling"] != "None":
        norm += f", {settings['postproc-scaling']} scaling"
    stat = f"{settings['postproc-stat-test']} test, {settings['postproc-fdr']} correction for multiple testing"
    return [f"**Filtering:** {filt}", f"**Missing values:** {imp}", f"**Normalization:** {norm}", f"**Statistics:** {stat}"]


def run_differential(
    pivot_df: pd.DataFrame,
    sample_group_map: dict,
    id_col: str,
    settings: dict,
    done: dict | None = None,
) -> dict:
    """Run the steps of the chain that are not in ``done`` and return all outputs.

    ``done`` maps step keys (see ``STEPS``) to tables the user already produced;
    the chain continues from the latest of them. Returns every step's table plus
    ``normalized_transform`` and the test actually used (``stat_test``).
    Raises ValueError when the samples do not form at least two groups.
    """
    done = {k: v for k, v in (done or {}).items() if k in STEPS and v is not None}
    # Samples without a group stay in the tables but take no part in the comparison
    sample_cols = [
        c for c in pivot_df.columns if sample_group_map.get(c, "Unassigned") != "Unassigned"
    ]
    metadata = pl.DataFrame(
        [{"sample_id": s, "group": sample_group_map[s]} for s in sample_cols],
        schema={"sample_id": pl.String, "group": pl.String},
    )
    groups = sorted({sample_group_map[s] for s in sample_cols})
    if len(groups) < 2:
        raise ValueError(
            "Comparing conditions needs runs of at least two conditions. Give each run a "
            "condition (e.g. control, treated)."
        )

    out = dict(done)
    out.setdefault("normalized_transform", settings["postproc-transform"])

    if "filtered_df" not in out:
        data = pl.from_pandas(strip_stat_columns(pivot_df)).lazy()
        method = settings["postproc-filter-method"]
        if method == "Low Abundance":
            data = filter_low_abundance(
                quantification_data=data, metadata=metadata, group_column="group",
                threshold_percentile=float(settings["postproc-filter-abundance-percentile"]),
            )
        elif method == "Low Variance":
            data = filter_low_variance(
                quantification_data=data, metadata=metadata, group_column="group",
                threshold_percentile=float(settings["postproc-filter-variance-percentile"]),
            )
        else:
            data = filter_low_repeatability(
                quantification_data=data, metadata=metadata, group_column="group",
                max_missing_ratio=float(settings["postproc-filter-max-missing"]) / 100.0,
            )
        out["filtered_df"] = strip_stat_columns(data.collect().to_pandas())
        if out["filtered_df"].empty:
            raise ValueError(
                "No proteoform passed the filter. Relax it (e.g. allow more missing values)."
            )

    if "imputed_df" not in out:
        data = pl.from_pandas(out["filtered_df"]).lazy()
        if settings["postproc-impute-class"].startswith("MAR"):
            data = impute_mar(
                quantification_data=data, metadata=metadata, group_column="group",
                strategy=settings["postproc-impute-mar-strategy"],
            )
        else:
            data = impute_smallest_value(
                quantification_data=data, metadata=metadata, scope=settings["postproc-impute-mnar-scope"],
            )
        out["imputed_df"] = strip_stat_columns(data.collect().to_pandas())

    if "normalized_df" not in out:
        normalization = settings["postproc-normalization"]
        if normalization == "reference_feature":
            # Needs a reference feature, which this page does not offer
            normalization = "median"
        data = transform_data(
            quantification_data=pl.from_pandas(out["imputed_df"]).lazy(),
            metadata=metadata, strategy=settings["postproc-transform"],
        )
        data = normalize_samples(
            quantification_data=data, metadata=metadata, strategy=normalization, id_col=id_col,
            reference_feature=None,
        )
        data = scale_data(quantification_data=data, metadata=metadata, strategy=settings["postproc-scaling"])
        out["normalized_df"] = strip_stat_columns(data.collect().to_pandas())
        out["normalized_transform"] = settings["postproc-transform"]

    out["stat_test"] = valid_test(settings["postproc-stat-test"], len(groups))
    if "statistics_df" not in out:
        data = calculate_statistical_tests(
            quantification_data=pl.from_pandas(out["normalized_df"]).lazy(),
            metadata=metadata, method=out["stat_test"],
        )
        data = adjust_fdr_lazy(quantification_data=data, strategy=settings["postproc-fdr"])
        out["statistics_df"] = data.collect().to_pandas()
    return out


def count_significant(statistics_df: pd.DataFrame, fc: float = 1.0, padj: float = 0.05) -> tuple[int, int]:
    """Features up and down at |log2FC| >= fc and adjusted p <= padj."""
    df = statistics_df.dropna(subset=["log2FC", "p-adj"])
    sig = df["p-adj"] <= padj
    return int((sig & (df["log2FC"] >= fc)).sum()), int((sig & (df["log2FC"] <= -fc)).sum())


def comparison_groups(sample_group_map: dict, sample_cols: list) -> tuple[str, str] | None:
    """(reference, compared) groups behind the sign of log2FC.

    openms_insight subtracts the alphabetically first group's mean from the
    last one's, so a positive log2FC means higher in the last group. Samples
    without a group are not part of the comparison.
    """
    groups = sorted({
        sample_group_map[s] for s in sample_cols
        if sample_group_map.get(s, "Unassigned") != "Unassigned"
    })
    if len(groups) < 2:
        return None
    return groups[0], groups[-1]
