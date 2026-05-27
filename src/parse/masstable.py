import atexit
import json
import os
import threading
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED

import pandas as pd
import polars as pl
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from pathlib import Path
from pyopenms import MSExperiment, MzMLFile, SpectrumLookup, Constants


# Fixed prefix of the deconv_dfs schema. Scoremap key columns are inserted between
# prefix and suffix at runtime, after pre-scanning the first deconvolved spectrum
# to discover which scoremap keys this FLASHDeconv run produced.
_DECONV_DFS_SCHEMA_FIXED = [
    ('rt', pa.float64()),
    ('ms_level', pa.int64()),
    ('mz_array', pa.list_(pa.float64())),
    ('intensity_array', pa.list_(pa.float32())),
    ('MinCharges', pa.list_(pa.int64())),
    ('MaxCharges', pa.list_(pa.int64())),
    ('MinIsotopes', pa.list_(pa.int64())),
    ('MaxIsotopes', pa.list_(pa.int64())),
    ('PrecursorScan', pa.float64()),
    ('PrecursorMass', pa.float64()),
]
_DECONV_DFS_SCHEMA_SUFFIX = [
    ('SignalPeaks',   pa.list_(pa.list_(pa.list_(pa.float64())))),
    ('NoisyPeaks',    pa.list_(pa.list_(pa.list_(pa.float64())))),
    ('CombinedPeaks', pa.list_(pa.list_(pa.list_(pa.float64())))),
    ('MSLevel', pa.int64()),
    ('Scan', pa.int64()),
]

_ANNO_DFS_SCHEMA = pa.schema([
    ('rt', pa.float64()),
    ('ms_level', pa.int64()),
    ('mz_array', pa.list_(pa.float64())),
    ('intensity_array', pa.list_(pa.float32())),
    ('MSLevel', pa.int64()),
])

# Batch size for per-spectrum dict accumulation before flushing to the parquet writer.
# Tunable; 500 keeps peak transient batch memory in the ~5-15 MB range on Eclip.
BATCH_SIZE = 500

# Cached pyopenms physical constants — accessed in the inner peak loop millions
# of times per Eclip run. Attribute lookup on the C-extension module is ~50-100 ns
# each; local-binding avoids the lookup overhead.
_PROTON_MASS_U = Constants.PROTON_MASS_U
_C13C12_MASSDIFF_U = Constants.C13C12_MASSDIFF_U

# Bounded-prefetch + reorder-before-flush constants for parallel parseFLASHDeconvOutput.
PREFETCH_PER_WORKER = 4
MAX_PENDING = 2000


def _default_worker_count():
    """Return the worker-pool size. Order of precedence:
    1. FLASH_POSTPROC_WORKERS env var (clamped to ≥1)
    2. settings.json max_threads (if present and > 0), capped at cpu_count-1 and 8
    3. min(cpu_count - 1, 8)
    """
    env_val = os.environ.get('FLASH_POSTPROC_WORKERS')
    if env_val:
        return max(1, int(env_val))
    cpu = (os.cpu_count() or 2) - 1
    try:
        settings_path = Path(__file__).resolve().parents[2] / 'settings.json'
        if settings_path.exists():
            raw = json.loads(settings_path.read_text())
            mt = raw.get('max_threads', 0)
            # max_threads may be a dict ({"local": N, "online": N}) or a plain int.
            if isinstance(mt, dict):
                n = mt.get('local', 0)
            else:
                n = mt
            if n and n > 0:
                return min(n, cpu, 8)
    except Exception:
        pass
    return min(cpu, 8)


_pool = None
_pool_workers = 0
_pool_lock = threading.Lock()


def _get_pool():
    """Return (executor, n_workers) for the module-level ProcessPoolExecutor
    singleton. Created lazily on first call; shut down via atexit hook.
    """
    global _pool, _pool_workers
    if _pool is None:
        with _pool_lock:
            if _pool is None:
                _pool_workers = _default_worker_count()
                _pool = ProcessPoolExecutor(max_workers=_pool_workers)
                atexit.register(_pool.shutdown, wait=True)
    return _pool, _pool_workers


def _parse_deconv_meta(mstr):
    """Parse the DeconvMassInfo meta string into a structured dict.

    Format example (semicolon-separated key=value pairs, with the 'peaks' value
    spilling across subsequent pair positions):
        tol=0;massoffset=0.000000;chargemass=1.007276;scorenames=cos,snr,qscore,qvalue;peaks=2:2,0:4,0.998392:3.97592:0.828587:1;
    Returns a dict with keys: tol, massoffset, chargemass, precursorscan,
    precursormass, peaks (list of [(minCharge,maxCharge),(minIso,maxIso)] tuples),
    scores (dict of key -> list[float]).
    """
    input_pairs = mstr.split(';')

    parsed_dict = {}
    scoremap = {}
    for pair in input_pairs:
        if len(pair) == 0:
            continue
        if '=' in pair:
            key, value = pair.split('=')
            if key == 'peaks':
                peaks_values = []
                peak_values = value.split(',')
                peaks_values.append([tuple(map(int, p.split(':'))) for p in peak_values])
                parsed_dict[key] = peaks_values
            elif ',' in value:  # scores (comma-separated list of floats with trailing comma)
                scoremap[key] = [float(x) for x in value[0:len(value) - 1].split(',')]
            else:
                parsed_dict[key] = float(value)
        else:
            peak_values = pair.split(',')
            parsed_dict['peaks'].append([tuple(map(int, p.split(':'))) for p in peak_values])

    return {
        'tol': parsed_dict['tol'],
        'massoffset': parsed_dict['massoffset'],
        'chargemass': parsed_dict['chargemass'],
        'precursorscan': parsed_dict['precursorscan'],
        'precursormass': parsed_dict['precursormass'],
        'peaks': parsed_dict['peaks'],
        'scores': scoremap,
    }


def _parse_anno_meta(mstr):
    """Parse the DeconvMassPeakIndices meta string into a list of (mass, signal_indices).

    Format example (semicolon-separated peak items, each `mass:idx,idx,idx`):
        1234.567:1,2,3;5678.9:4,5;
    Returns list[(float, list[int])].
    """
    peak_items = mstr.split(';')
    parsed_peaks = []
    for item in peak_items:
        if len(item) == 0:
            continue
        peak_values = item.split(':')
        peak_mass = float(peak_values[0])
        peak_infos = list(map(int, peak_values[1].split(',')))
        parsed_peaks.append((peak_mass, peak_infos))
    return parsed_peaks


def _extract_scan(native_id, source_files, fallback_idx):
    """Extract the integer scan number from a spectrum's native ID.

    Uses pyopenms SpectrumLookup; falls back to fallback_idx (1-based spec ordinal)
    if extraction returns <= 0. type_accession is read from the first source file
    when available, defaulting to MS:1000768 otherwise.
    """
    type_accession = 'MS:1000768'
    if source_files:
        type_accession = source_files[0].getNativeIDTypeAccession()
        if not type_accession:
            type_accession = 'MS:1000768'
    scan_number = SpectrumLookup().extractScanNumber(native_id, type_accession)
    return scan_number if scan_number > 0 else fallback_idx


def _validate_scoremap_keys(scoremap, expected, spec_idx):
    """Raise ValueError if the per-spectrum scoremap keys diverge from the
    first-spectrum reference set. Fail-fast replacement of the old code's
    silent lazy-append behavior."""
    got = set(scoremap.keys())
    if got != expected:
        raise ValueError(
            f"Spectrum {spec_idx}: scoremap keys diverged from first-spectrum set. "
            f"Expected {expected}, got {got}"
        )


def _compute_peak_cells(payload):
    """Worker function: compute (sig_cells, noisy_cells) for one spectrum from a
    picklable payload. Pure Python + numpy + Constants. No MSExperiment access.

    Mirrors the per-spectrum nested loop in parseFLASHDeconvOutput but uses
    np.searchsorted to replicate pyopenms findNearest(mz). findNearest returns
    the index of the peak with minimum absolute distance to mz; we replicate
    this by comparing the two candidates around the searchsorted insertion point
    and picking the closer one.
    """
    if os.environ.get('FLASH_POSTPROC_FORCE_CRASH'):
        raise RuntimeError("simulated worker crash (FLASH_POSTPROC_FORCE_CRASH set)")
    anno_masses   = payload['anno_masses']
    parsed_peaks  = payload['parsed_peaks']
    deconv_masses = payload['deconv_masses']
    aspec_mz      = payload['aspec_mz_sorted']
    aspec_int     = payload['aspec_int_sorted']
    aspec_size    = len(aspec_mz)

    sig_cells, noisy_cells = [], []
    if aspec_size == 0:
        # Defensive: an empty annotated spectrum has no peaks to match; emit
        # empty cells for every mass without iterating the z-loop.
        for _ in range(len(anno_masses)):
            sig_cells.append([])
            noisy_cells.append([])
        return sig_cells, noisy_cells

    for mass_idx in range(len(anno_masses)):
        peakinfo = parsed_peaks[mass_idx]
        anno_mass, sig_indices = anno_masses[mass_idx]
        sig_set = set(sig_indices)
        deconv_mass = float(deconv_masses[mass_idx])

        sig_buf, noisy_buf = [], []
        for z in range(peakinfo[0][0], peakinfo[0][1] + 1):
            minmz = (deconv_mass - 3.0) / z + _PROTON_MASS_U
            maxmz = (deconv_mass + 3.0 + peakinfo[1][1] * _C13C12_MASSDIFF_U) / z + _PROTON_MASS_U
            # Replicate pyopenms findNearest: index of minimum |mz - minmz|.
            # searchsorted gives the insertion point; the nearest is whichever
            # of the two surrounding candidates is closer.
            pos = int(np.searchsorted(aspec_mz, minmz, side='left'))
            if pos == 0:
                start = 0
            elif pos == aspec_size:
                start = aspec_size - 1
            else:
                lo, hi = pos - 1, pos
                start = lo if abs(float(aspec_mz[lo]) - minmz) <= abs(float(aspec_mz[hi]) - minmz) else hi
            for pi in range(start, aspec_size):
                pmz = float(aspec_mz[pi])
                if pmz > maxmz:
                    break
                if z == round(deconv_mass / pmz):
                    rec = (float(pi), pmz, float(aspec_int[pi]),
                           float(round(anno_mass / pmz)))
                    (sig_buf if pi in sig_set else noisy_buf).append(rec)
        sig_cells.append(sig_buf)
        noisy_cells.append(noisy_buf)
    return sig_cells, noisy_cells


def parseFLASHDeconvOutput(annotated, deconvolved, file_manager, dataset_id,
                           deconv_tag='deconv_dfs', anno_tag='anno_dfs', logger=None):
    """Stream-write deconv_dfs and anno_dfs parquets via file_manager; return the
    deconv tolerance (single float).

    Parallel implementation: worker pool computes the inner SignalPeaks/NoisyPeaks
    loop per spectrum; parent assembles EASY fields inline and reorders worker
    completions by spec_idx before flushing to keep parquet rows bit-identical
    to the sequential version.
    """
    annotated_exp = MSExperiment()
    deconvolved_exp = MSExperiment()
    MzMLFile().load(str(Path(annotated)), annotated_exp)
    MzMLFile().load(str(Path(deconvolved)), deconvolved_exp)

    # Pre-scan first deconvolved spectrum to discover the scoremap key set.
    first_parsed = _parse_deconv_meta(deconvolved_exp[0].getMetaValue('DeconvMassInfo'))
    scoremap_keys_ordered = list(first_parsed['scores'].keys())
    scoremap_keys = set(scoremap_keys_ordered)
    deconv_schema = pa.schema(
        _DECONV_DFS_SCHEMA_FIXED
        + [(k, pa.list_(pa.float64())) for k in scoremap_keys_ordered]
        + _DECONV_DFS_SCHEMA_SUFFIX
    )

    source_files = annotated_exp.getSourceFiles()
    executor, n_workers = _get_pool()
    prefetch = n_workers * PREFETCH_PER_WORKER
    spec_iter = enumerate(zip(deconvolved_exp, annotated_exp))

    in_flight = {}              # Future -> spec_idx
    parent_state = {}           # spec_idx -> easy_fields dict
    pending_complete = {}       # spec_idx -> completed deconv row dict
    next_flush_idx = [0]        # single-element list for closure write access
    ready_deconv = []
    anno_batch = []
    tol_by_idx = {}

    with file_manager.parquet_sink(dataset_id, deconv_tag) as deconv_path, \
         file_manager.parquet_sink(dataset_id, anno_tag) as anno_path, \
         pq.ParquetWriter(deconv_path, deconv_schema) as deconv_writer, \
         pq.ParquetWriter(anno_path, _ANNO_DFS_SCHEMA) as anno_writer:

        def submit_one():
            try:
                spec_idx, (spec, aspec) = next(spec_iter)
            except StopIteration:
                return False
            spec.sortByPosition()
            aspec.sortByPosition()
            parsed = _parse_deconv_meta(spec.getMetaValue('DeconvMassInfo'))
            anno_masses = _parse_anno_meta(aspec.getMetaValue('DeconvMassPeakIndices'))
            _validate_scoremap_keys(parsed['scores'], expected=scoremap_keys, spec_idx=spec_idx)
            tol_by_idx[spec_idx] = parsed['tol']
            spec_mz, spec_int = spec.get_peaks()
            aspec_mz, aspec_int = aspec.get_peaks()

            easy = {
                'rt': spec.getRT(),
                'ms_level': spec.getMSLevel(),
                'mz_array': spec_mz,
                'intensity_array': spec_int,
                'MinCharges': [p[0][0] for p in parsed['peaks']],
                'MaxCharges': [p[0][1] for p in parsed['peaks']],
                'MinIsotopes': [p[1][0] for p in parsed['peaks']],
                'MaxIsotopes': [p[1][1] for p in parsed['peaks']],
                'PrecursorScan': parsed['precursorscan'],
                'PrecursorMass': parsed['precursormass'],
                'MSLevel': aspec.getMSLevel(),
                'Scan': _extract_scan(aspec.getNativeID(), source_files,
                                      fallback_idx=spec_idx + 1),
            }
            for k in scoremap_keys_ordered:
                easy[k] = parsed['scores'][k]
            parent_state[spec_idx] = easy

            anno_batch.append({
                'rt': aspec.getRT(),
                'ms_level': aspec.getMSLevel(),
                'mz_array': aspec_mz,
                'intensity_array': aspec_int,
                'MSLevel': aspec.getMSLevel(),
            })

            payload = {
                'anno_masses': anno_masses,
                'parsed_peaks': parsed['peaks'],
                'deconv_masses': spec_mz[:len(anno_masses)],
                'aspec_mz_sorted': aspec_mz,
                'aspec_int_sorted': aspec_int,
            }
            fut = executor.submit(_compute_peak_cells, payload)
            in_flight[fut] = spec_idx
            return True

        def harvest(fut):
            spec_idx = in_flight.pop(fut)
            sig_cells, noisy_cells = fut.result()
            easy = parent_state.pop(spec_idx)
            easy['SignalPeaks'] = sig_cells
            easy['NoisyPeaks'] = noisy_cells
            easy['CombinedPeaks'] = noisy_cells   # same data, separate column for wire compat
            pending_complete[spec_idx] = easy

            while next_flush_idx[0] in pending_complete:
                ready_deconv.append(pending_complete.pop(next_flush_idx[0]))
                next_flush_idx[0] += 1

            if len(ready_deconv) >= BATCH_SIZE:
                deconv_writer.write_table(pa.Table.from_pylist(ready_deconv, schema=deconv_schema))
                ready_deconv.clear()
            if len(anno_batch) >= BATCH_SIZE:
                anno_writer.write_table(pa.Table.from_pylist(anno_batch, schema=_ANNO_DFS_SCHEMA))
                anno_batch.clear()

        # Prime the prefetch window
        for _ in range(prefetch):
            if not submit_one():
                break

        # Sliding-window drain with backpressure
        while in_flight:
            done, _ = wait(list(in_flight), return_when=FIRST_COMPLETED)
            for fut in done:
                harvest(fut)
                if len(pending_complete) < MAX_PENDING:
                    submit_one()

        # Final flush
        if ready_deconv:
            deconv_writer.write_table(pa.Table.from_pylist(ready_deconv, schema=deconv_schema))
        if anno_batch:
            anno_writer.write_table(pa.Table.from_pylist(anno_batch, schema=_ANNO_DFS_SCHEMA))

    return tol_by_idx[max(tol_by_idx)] if tol_by_idx else 0.0


def parseFLASHTaggerOutput(tags, proteins):
    # db = get_sequences(FastaFile.read(db), ProteinSequence)
    return pd.read_csv(tags, sep='\t'), pd.read_csv(proteins, sep='\t')


def getSpectraTableDF(deconv_df: pd.DataFrame):
    out_df = deconv_df[['Scan', 'MSLevel', 'rt', 'PrecursorMass']].copy()
    out_df['#Masses'] = [len(ele) for ele in deconv_df['MinCharges']]
    out_df.reset_index(inplace=True)
    return out_df


def getMSSignalDF(anno_df: pd.DataFrame):
    # Convert to polars for efficient processing
    # pl_df = pl.from_pandas(anno_df.reset_index(names=['scan_idx']))
    pl_df = anno_df.with_row_count(name="scan_idx")

    # Create exploded dataframe with polars efficient operations
    exploded_df = (
        pl_df
        .with_columns([
            # Convert numpy arrays to polars lists
            pl.col("mz_array").map_elements(lambda x: x.tolist() if hasattr(x, 'tolist') else list(x), return_dtype=pl.List(pl.Float32)),
            pl.col("intensity_array").map_elements(lambda x: x.tolist() if hasattr(x, 'tolist') else list(x), return_dtype=pl.List(pl.Float32))
        ])
        .with_row_index("original_idx")
        .select([
            pl.col("scan_idx"),
            pl.col("MSLevel"),
            pl.col("rt"),
            pl.col("mz_array").alias("mass"),
            pl.col("intensity_array").alias("intensity")
        ])
        # Explode arrays to individual rows
        .explode(["mass", "intensity"])
        # Add mass index within each scan
        .with_columns([
            pl.int_range(pl.len()).over("scan_idx").alias("mass_idx")
        ])
        # Filter out NaN and zero intensities
        .filter(
            pl.col("intensity").is_not_null() &
            (pl.col("intensity") > 0)
        )
        # Sort by intensity
        .sort("intensity")
        # Select final columns in correct order
        .select([
            pl.col("mass"),
            pl.col("rt"),
            pl.col("intensity"),
            pl.col("scan_idx"),
            pl.col("mass_idx"),
            pl.col("MSLevel")
        ])
    )

    # Return polars LazyFrame for efficient pipeline processing
    return exploded_df.lazy()
