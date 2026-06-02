import pandas as pd
import polars as pl


def build_proteoform_scan_map(protein_df, scan_table_df):
    """Map each proteoform index to its scan and the deconv row index.

    protein_df: DataFrame with 'index' (proteoform index) and 'Scan'.
    scan_table_df: DataFrame with 'index' (deconv row index) and 'Scan'.

    Returns {proteoform_index: {'scan': int, 'deconv_index': int}}.
    Proteoforms whose Scan is NaN or absent from scan_table are omitted.
    """
    scan_to_index = (
        scan_table_df.drop_duplicates(subset="Scan", keep="first")
        .set_index("Scan")["index"]
    )
    result = {}
    for proteoform_index, scan in zip(protein_df["index"], protein_df["Scan"]):
        if pd.isna(scan):
            continue
        scan = int(scan)
        if scan in scan_to_index.index:
            result[int(proteoform_index)] = {
                "scan": scan,
                "deconv_index": int(scan_to_index.loc[scan]),
            }
    return result


def build_proteoform_scan_frame(protein_df, scan_table_df):
    """Polars frame surfacing the proteoform->scan resolution as COLUMNS.

    ADDITIVE helper for the OpenMS-Insight TnT viewer (Stage C). It reproduces
    ``build_proteoform_scan_map`` (the legacy PyArrow pushdown in
    ``src/render/update.py``) as a value-filterable frame so OpenMS-Insight
    components can value-filter (``filters={'proteinIndex': 'proteoform_index'}``)
    instead of doing an ``iloc`` / per-scan pushdown by hand.

    Args:
        protein_df: pandas/polars frame with columns ``index`` (proteoform index)
            and ``Scan``.
        scan_table_df: pandas/polars frame with columns ``index`` (deconv row
            index) and ``Scan``.

    Returns:
        Polars DataFrame with columns ``proteoform_index`` (int64), ``scan``
        (int64) and ``deconv_index`` (int64). Proteoforms whose ``Scan`` is null
        or absent from ``scan_table_df`` are omitted (same as the map builder).
    """
    scan_map = build_proteoform_scan_map(
        _to_pandas(protein_df), _to_pandas(scan_table_df)
    )
    if not scan_map:
        return pl.DataFrame(
            schema={
                "proteoform_index": pl.Int64,
                "scan": pl.Int64,
                "deconv_index": pl.Int64,
            }
        )
    rows = [
        {"proteoform_index": int(pid), "scan": int(v["scan"]),
         "deconv_index": int(v["deconv_index"])}
        for pid, v in scan_map.items()
    ]
    return pl.DataFrame(rows).sort("proteoform_index")


def _to_pandas(df):
    """Accept a polars or pandas frame; return pandas (the map builder uses
    pandas indexing semantics)."""
    if isinstance(df, pl.DataFrame):
        return df.to_pandas()
    if isinstance(df, pl.LazyFrame):
        return df.collect().to_pandas()
    return df
