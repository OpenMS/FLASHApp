import pandas as pd


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
