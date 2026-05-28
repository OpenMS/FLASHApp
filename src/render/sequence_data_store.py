import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq

# One row per proteoform. Explicit schema so the always-empty
# fixed_modifications and the empty/variable modifications get consistent types.
SCHEMA = pa.schema([
    ("proteoform_index", pa.int64()),
    ("sequence", pa.list_(pa.string())),
    ("theoretical_mass", pa.float64()),
    ("fixed_modifications", pa.list_(pa.string())),
    ("coverage", pa.list_(pa.float64())),
    ("maxCoverage", pa.float64()),
    ("fragment_masses_a", pa.list_(pa.list_(pa.float64()))),
    ("fragment_masses_b", pa.list_(pa.list_(pa.float64()))),
    ("fragment_masses_c", pa.list_(pa.list_(pa.float64()))),
    ("fragment_masses_x", pa.list_(pa.list_(pa.float64()))),
    ("fragment_masses_y", pa.list_(pa.list_(pa.float64()))),
    ("fragment_masses_z", pa.list_(pa.list_(pa.float64()))),
    ("proteoform_start", pa.int64()),
    ("proteoform_end", pa.int64()),
    ("computed_mass", pa.float64()),
    ("modifications", pa.list_(pa.struct([
        ("start", pa.int64()), ("end", pa.int64()),
        ("mass_diff", pa.float64()), ("labels", pa.string()),
    ]))),
])

ROW_GROUP_SIZE = 64
ENTRY_KEYS = [f.name for f in SCHEMA if f.name != "proteoform_index"]


def _py(x):
    """Recursively convert numpy scalars to builtins so pa.Table.from_pylist
    serializes cleanly (coverage/maxCoverage are np.float64)."""
    import numpy as np
    if isinstance(x, np.generic):
        return x.item()
    if isinstance(x, list):
        return [_py(v) for v in x]
    if isinstance(x, dict):
        return {k: _py(v) for k, v in x.items()}
    return x


def build_table(sequence_data):
    """{proteoform_index: entry} -> pyarrow Table, one row per proteoform,
    sorted by proteoform_index (so row groups carry contiguous index ranges
    and pushdown can skip)."""
    rows = []
    for pid in sorted(sequence_data):
        entry = sequence_data[pid]
        row = {"proteoform_index": int(pid)}
        for k in ENTRY_KEYS:
            row[k] = _py(entry[k])
        rows.append(row)
    return pa.Table.from_pylist(rows, schema=SCHEMA)


def _as_dataset(dataset_or_path):
    if isinstance(dataset_or_path, ds.Dataset):
        return dataset_or_path
    return ds.dataset(str(dataset_or_path), format="parquet")


def load_entry(dataset_or_path, proteoform_index):
    """Pushdown-read one proteoform's row; return its entry dict (native Python
    containers via to_pylist) with proteoform_index removed, or None if absent."""
    dataset = _as_dataset(dataset_or_path)
    table = dataset.to_table(filter=ds.field("proteoform_index") == int(proteoform_index))
    rows = table.to_pylist()
    if not rows:
        return None
    entry = rows[0]
    entry.pop("proteoform_index", None)
    return entry


def reconstruct_all(dataset_or_path):
    """Read every row -> {proteoform_index: entry}. For migration verification
    and the golden adapter only; never the hot render path."""
    if isinstance(dataset_or_path, ds.Dataset):
        table = dataset_or_path.to_table()
    else:
        table = pq.read_table(str(dataset_or_path))
    out = {}
    for row in table.to_pylist():
        pid = row.pop("proteoform_index")
        out[pid] = row
    return out
