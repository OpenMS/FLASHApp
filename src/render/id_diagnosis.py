"""Why a run has few identifications: its MS2 spectra step by step to a PrSM.

Every MS2 spectrum of an identified run passes the same steps on its way to
a proteoform: FLASHDeconv finds masses in it (and the precursor's mass in the
MS1 spectrum before it), FLASHTnT finds sequence tags among those masses,
matches the spectrum to a proteoform of the database, and the target/decoy
FDR keeps or drops the match. ``steps`` counts the spectra left after each
step; ``diagnose`` names, for each step that loses more than is usual, the
likely cause and the setting to look at, the most telling step first.

Everything comes from results the workflow already stored (scan table,
proteoform table, tag table, annotation settings). Streamlit-free, so tests
call it.
"""

import polars as pl

from src.render import insight_data

Q_THRESHOLD = 0.01
DECOY_PREFIX = 'DECOY'

# Shares (of the spectra entering a step) below which a step is flagged
MIN_DECONVOLVED = 0.5
MIN_PRECURSOR = 0.5
MIN_TAGGED = 0.5
MIN_MATCHED = 0.05
MIN_KEPT = 0.5
# A median fragment calibration error beyond this share of the tolerance is an offset
OFFSET_SHARE = 0.5
MAX_CAUSES = 3

COLUMNS = ['step', 'spectra', 'share', 'cause']
STEP_LABELS = {
    'ms2': 'MS2 spectra',
    'deconvolved': 'with deconvolved masses',
    'precursor': 'with a precursor mass',
    'tagged': 'with sequence tags',
    'matched': 'matched to a proteoform',
    'kept': f'PrSMs at q ≤ {Q_THRESHOLD:g}',
}


def _scans(frame, column='Scan'):
    return set(frame[column].drop_nulls().to_list()) if frame.height else set()


def steps(file_manager, dataset_id):
    """Counts of an identified run's MS2 spectra after every step, with what
    the diagnosis needs besides (MS1 spectra, q-values estimated, tolerance,
    median fragment calibration error). None without a scan table."""
    if insight_data.result_path(file_manager, dataset_id, 'scan_table') is None:
        return None
    scans = insight_data.scan_table(file_manager, dataset_id).select(
        'Scan', 'MSLevel', 'PrecursorMass', '#Masses').collect()
    ms2 = scans.filter(pl.col('MSLevel') == 2)
    ms2_scans = _scans(ms2)
    counts = {
        'ms1': scans.filter(pl.col('MSLevel') == 1).height,
        'ms2': ms2.height,
        'deconvolved': ms2.filter(pl.col('#Masses') > 0).height,
        'precursor': ms2.filter(pl.col('PrecursorMass') > 0).height,
        'tagged': None, 'matched': None, 'kept': None,
        'q_estimated': None, 'calibration_ppm': None, 'tolerance_ppm': None,
    }
    if insight_data.result_path(file_manager, dataset_id, 'tag_dfs') is not None:
        tags = insight_data.scan(file_manager, dataset_id, 'tag_dfs').select('Scan').collect()
        counts['tagged'] = len(_scans(tags) & ms2_scans)
    if insight_data.result_path(file_manager, dataset_id, 'protein_dfs') is not None:
        lf = insight_data.scan(file_manager, dataset_id, 'protein_dfs')
        names = lf.collect_schema().names()
        hits = lf.select(
            'Scan', 'accession', 'ProteoformLevelQvalue',
            *(['RawCalibrationPPM'] if 'RawCalibrationPPM' in names else []),
        ).collect()
        hits = hits.filter(~pl.col('accession').fill_null('').str.starts_with(DECOY_PREFIX))
        q = pl.col('ProteoformLevelQvalue')
        counts['matched'] = len(_scans(hits) & ms2_scans)
        counts['q_estimated'] = bool(hits.filter(q >= 0).height) if hits.height else None
        counts['kept'] = len(_scans(hits.filter((q >= 0) & (q <= Q_THRESHOLD))) & ms2_scans)
        if 'RawCalibrationPPM' in hits.columns:
            # FLASHAnnotator writes -1 when it did not calibrate a match
            calibration = hits.filter(pl.col('RawCalibrationPPM') != -1)['RawCalibrationPPM']
            if calibration.len():
                counts['calibration_ppm'] = float(calibration.median())
    if insight_data.result_path(file_manager, dataset_id, 'settings') is not None:
        settings = file_manager.get_results(dataset_id, ['settings'])['settings'] or {}
        if settings.get('tolerance'):
            counts['tolerance_ppm'] = float(settings['tolerance'])
    return counts


def _share(part, whole):
    return None if part is None or not whole else part / whole


def _percent(share):
    return f'{100 * share:.0f} %'


def _causes(counts):
    """{step: (weight, text)} for the steps that lose more than is usual."""
    causes = {}
    ms2 = counts['ms2']
    if not ms2:
        causes['ms2'] = (1.0, 'The run has no MS2 spectra; proteoforms are identified from '
                              'fragmentation spectra. An intact-mass run is deconvolved only.')
        return causes
    deconvolved = _share(counts['deconvolved'], ms2)
    if deconvolved is not None and deconvolved < MIN_DECONVOLVED:
        causes['deconvolved'] = (0.9, f'FLASHDeconv found masses in only {_percent(deconvolved)} '
                                      'of the MS2 spectra: check that the spectra are centroided and '
                                      'that the charge and mass ranges of MS2 cover the fragments.')
    precursor = _share(counts['precursor'], ms2)
    if precursor is not None and precursor < MIN_PRECURSOR:
        if not counts['ms1']:
            text = ('The run has no MS1 spectra, so no MS2 spectrum has a precursor mass; '
                    'the proteoform mass comes from the fragments alone.')
        else:
            text = (f'Only {_percent(precursor)} of the MS2 spectra have a precursor mass: '
                    'FLASHDeconv found no mass in their isolation window of the MS1 spectrum '
                    'before them (MS1 charge and mass ranges, or a low precursor signal).')
        causes['precursor'] = (0.5, text)
    tagged = _share(counts['tagged'], counts['deconvolved'])
    if tagged is not None and tagged < MIN_TAGGED:
        causes['tagged'] = (0.8, f'Only {_percent(tagged)} of the spectra with masses gave a sequence '
                                 'tag: a shorter minimum tag length or a wider tolerance finds more.')
    matched = _share(counts['matched'], counts['tagged'] or counts['deconvolved'])
    if matched is not None and matched < MIN_MATCHED:
        causes['matched'] = (0.7, f'Only {_percent(matched)} of the spectra with tags matched a '
                                  'proteoform: check that the database holds the sample\'s proteins '
                                  '(organism, isoforms) and that its modifications are searched.')
    if counts['q_estimated'] is False:
        causes['kept'] = (0.6, 'FLASHTnT estimated no q-values (no decoys were searched), so no '
                               'PrSM passes an FDR threshold.')
    else:
        kept = _share(counts['kept'], counts['matched'])
        if kept is not None and kept < MIN_KEPT:
            causes['kept'] = (0.6, f'Only {_percent(kept)} of the matched spectra pass q ≤ '
                                   f'{Q_THRESHOLD:g}: decoys score like the targets, as with a database '
                                   'much larger than the sample or many unknown modifications allowed.')
    error, tolerance = counts['calibration_ppm'], counts['tolerance_ppm']
    if error is not None and tolerance and abs(error) > OFFSET_SHARE * tolerance:
        text = (f'The median fragment error of the matches is {error:+.1f} ppm in a ±{tolerance:g} ppm '
                'tolerance: the instrument\'s calibration is off, so a wider tolerance finds more.')
        step = 'matched' if 'matched' not in causes else 'kept' if 'kept' not in causes else None
        if step is not None:
            causes[step] = (0.75, text)
    return causes


def diagnose(counts):
    """One row per step (COLUMNS): its spectra, their share of the MS2
    spectra and, where the step loses more than is usual, the likely cause.
    Only the MAX_CAUSES most telling causes are named."""
    causes = _causes(counts)
    named = sorted(causes, key=lambda step: -causes[step][0])[:MAX_CAUSES]
    rows = []
    for step, label in STEP_LABELS.items():
        count = counts.get(step)
        if count is None:
            continue
        share = _share(count, counts['ms2'])
        rows.append({
            'step': label,
            'spectra': int(count),
            'share': None if share is None else round(100 * share, 1),
            'cause': causes[step][1] if step in named else '',
        })
    return pl.DataFrame(rows, schema={'step': pl.String, 'spectra': pl.Int64,
                                      'share': pl.Float64, 'cause': pl.String})


def headline(counts):
    """One sentence on the run: its PrSM rate and how many steps lose spectra."""
    if not counts['ms2']:
        return 'This run has no MS2 spectra to identify.'
    kept = counts['kept'] if counts['kept'] is not None else counts['matched']
    rate = _share(kept, counts['ms2'])
    flagged = len(_causes(counts))
    text = (f'{kept:,} of {counts["ms2"]:,} MS2 spectra ({_percent(rate)}) '
            f'became a PrSM at q ≤ {Q_THRESHOLD:g}.') if rate is not None else ''
    if flagged:
        return f'{text} {flagged} step{"s" if flagged > 1 else ""} lose more spectra than usual.'.strip()
    return f'{text} No step loses more spectra than usual.'.strip()
