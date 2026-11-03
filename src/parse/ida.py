"""Results of FLASHIda acquisitions and of the FLASHIda simulator, stored for
the FLASHIda viewers (src/render/render.py)."""

import pandas as pd

from src.render.compression import compute_compression_levels, downsample_heatmap

# Columns of the MS2 table, from FLASHDeconv's spec2.tsv (one row per mass;
# the precursor columns are the same for every mass of a scan)
MS2_COLUMNS = {
    'ScanNum': 'Scan', 'RetentionTime': 'RT', 'PrecursorScanNum': 'PrecursorScan',
    'PrecursorMz': 'PrecursorMz', 'PrecursorCharge': 'PrecursorCharge',
    'PrecursorMonoisotopicMass': 'PrecursorMass', 'PrecursorIntensity': 'PrecursorIntensity',
    'PrecursorQscore': 'PrecursorQscore', 'MassCountInSpec': '#Masses',
}


def identified_scans(id_tsv) -> set:
    """Scan numbers with a PrSM in a TopPIC *_prsm_single.tsv. Its parameter
    header has no fixed length, so the table starts at the line naming the
    "Scan(s)" column."""
    with open(id_tsv, encoding='utf-8', errors='replace') as f:
        lines = f.readlines()
    start = next((i for i, line in enumerate(lines) if 'Scan(s)' in line.split('\t')), None)
    if start is None:
        return set()
    ids = pd.read_csv(id_tsv, sep='\t', skiprows=start, dtype=str)
    scans = set()
    for value in ids['Scan(s)'].dropna():
        scans.update(int(s) for s in str(value).split() if s.strip().isdigit())
    return scans


def ms2_table(spec2: pd.DataFrame, scan_table: pd.DataFrame, identified=None) -> pd.DataFrame:
    """One row per MS2 scan with its precursor, `index` (the spectrum index
    every viewer selection uses) and `precursor_index` (the precursor
    spectrum's index; -1 when it is not in the file)."""
    present = [c for c in MS2_COLUMNS if c in spec2.columns]
    table = (
        spec2.drop_duplicates(subset='ScanNum', keep='first')[present]
        .rename(columns=MS2_COLUMNS)
    )
    index_of = dict(zip(scan_table['Scan'], scan_table['index']))
    table.insert(0, 'index', table['Scan'].map(index_of))
    table = table.dropna(subset=['index'])
    table['index'] = table['index'].astype('int64')
    if 'PrecursorScan' in table.columns:
        table['precursor_index'] = (
            table['PrecursorScan'].map(index_of).fillna(-1).astype('int64')
        )
    else:
        table['precursor_index'] = -1
    if identified is not None:
        table['Identified'] = table['Scan'].isin(identified)
    return table.reset_index(drop=True)


def parseIda(file_manager, dataset_id, spec2_tsv, id_tsv=None):
    """Store the MS2 table of a FLASHIda acquisition. Needs the dataset's
    scan table (parseDeconv)."""
    spec2 = pd.read_csv(spec2_tsv, sep='\t')
    scan_table = file_manager.get_results(dataset_id, ['scan_table'])['scan_table']
    identified = None if id_tsv is None else identified_scans(id_tsv)
    file_manager.store_data(
        dataset_id, 'ida_ms2_table', ms2_table(spec2, scan_table, identified)
    )


def simulation_heatmap(simulation: pd.DataFrame) -> pd.DataFrame:
    """The simulated precursors as heatmap points; `row` links them to the
    simulation table."""
    heatmap = simulation.loc[:, ['monoMasses', 'rt', 'precursorIntensity']].rename(
        columns={'monoMasses': 'mass', 'precursorIntensity': 'intensity'}
    )
    heatmap.insert(0, 'row', range(len(heatmap)))
    return heatmap.dropna(subset=['mass', 'rt', 'intensity'])


def parseIdaSimulator(file_manager, dataset_id, out_simulation, logger=None):
    """Store a FLASHIda simulator result: the table, and the heatmap of the
    precursors it chose with its downsampled levels."""
    simulation = pd.read_csv(out_simulation, sep='\t')
    simulation = simulation.loc[:, [c for c in simulation.columns if not c.startswith('Unnamed')]]
    simulation.insert(0, 'row', range(len(simulation)))
    file_manager.store_data(dataset_id, 'simulation_dfs', simulation)

    heatmap = simulation_heatmap(simulation.drop(columns='row'))
    file_manager.store_data(dataset_id, 'simulation_heatmap', heatmap)
    for size in reversed(compute_compression_levels(20000, len(heatmap), logger=logger)):
        heatmap = downsample_heatmap(heatmap, max_datapoints=size)
        file_manager.store_data(dataset_id, f'simulation_heatmap_{size}', heatmap)
