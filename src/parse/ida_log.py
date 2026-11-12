"""The logs of a FLASHIda run, read while the run is going on.

FLASHIda writes every run into a folder of its own,
``<runtime.log_dir>/<raw name>_<yyyy-MM-dd-HH-mm-ss>/``, with five streams:
``ida.log`` and four tab-separated tables. The viewer reads the tables:

- ``scan_commands.tsv``: one row per scan FLASHIda commanded (survey, recording,
  exploration, quantification, ...), with the precursor it chose and why
  (qscore, charge, isolation window, activation and energy);
- ``scan_results.tsv``: one row per scan the instrument returned, with its
  retention time (minutes), the masses deconvolved from it and the outcome of
  explorations and quantification;
- ``identification.tsv``: the proteoform a scan identified;
- ``pooled_identification.tsv``: the proteoform of each precursor with the
  fragments pooled over all of its scans.

The files grow while FLASHIda runs, so a last line without a line break is not
read yet. Streamlit-free.
"""

import csv
import hashlib
import io
from pathlib import Path

import numpy as np
import pandas as pd

STREAMS = (
    'scan_commands.tsv', 'scan_results.tsv', 'identification.tsv',
    'pooled_identification.tsv',
)


def is_run(folder) -> bool:
    return Path(folder, 'scan_commands.tsv').is_file()


def find_runs(log_dir) -> list:
    """Run folders in a log folder (or the folder itself if it is one),
    newest first."""
    log_dir = Path(log_dir)
    if not log_dir.is_dir():
        return []
    if is_run(log_dir):
        return [log_dir]
    runs = [p for p in log_dir.iterdir() if p.is_dir() and is_run(p)]
    return sorted(runs, key=lambda p: Path(p, 'scan_commands.tsv').stat().st_mtime, reverse=True)


def stamp(run) -> str:
    """Changes whenever one of the run's tables does."""
    parts = []
    for name in STREAMS:
        try:
            stat = Path(run, name).stat()
            parts.append(f'{name}:{stat.st_size}:{stat.st_mtime_ns}')
        except OSError:
            parts.append(f'{name}:-')
    return hashlib.md5('|'.join(parts).encode()).hexdigest()[:12]


def last_change(run):
    """Time of the last write to the run's tables (seconds since the epoch)."""
    times = []
    for name in STREAMS:
        try:
            times.append(Path(run, name).stat().st_mtime)
        except OSError:
            pass
    return max(times) if times else None


def read_table(path) -> pd.DataFrame:
    """A log table as strings. A last line still being written is left out."""
    try:
        text = Path(path).read_text(encoding='utf-8', errors='replace')
    except OSError:
        return pd.DataFrame()
    if not text.endswith('\n'):
        text = text[:text.rfind('\n') + 1]
    if not text.strip():
        return pd.DataFrame()
    rows = list(csv.reader(io.StringIO(text), delimiter='\t'))
    header, rows = rows[0], [r for r in rows[1:] if len(r) == len(rows[0])]
    return pd.DataFrame(rows, columns=header, dtype=str)


def read_run(run) -> dict:
    return {name.split('.')[0]: read_table(Path(run, name)) for name in STREAMS}


def _num(df, column, default=np.nan):
    if column not in df.columns:
        return pd.Series(default, index=df.index, dtype='float64')
    return pd.to_numeric(df[column], errors='coerce')


def _col(df, column, default=''):
    if column not in df.columns:
        return pd.Series(default, index=df.index, dtype=object)
    return df[column]


def _first(df, column):
    """First value of a column whose MS3 rows list one value per stage."""
    return pd.to_numeric(_col(df, column).str.split(';').str[0], errors='coerce')


def _last(df, column):
    return pd.to_numeric(_col(df, column).str.split(';').str[-1], errors='coerce')


def survey_masses(results: pd.DataFrame) -> pd.DataFrame:
    """Every mass deconvolved from a survey scan: rt (minutes), mass,
    intensity (summed over its charges) and qscore (NaN when not logged)."""
    surveys = results[_col(results, 'ms_level') == '1']
    rt, masses, intensities, qscores = [], [], [], []
    for row_rt, row_masses, row_intensities, row_qscores in zip(
        _num(surveys, 'rt'), _col(surveys, 'deconv_masses'), _col(surveys, 'deconv_intensities'),
        _col(surveys, 'deconv_qscores'),
    ):
        if not row_masses:
            continue
        groups = row_intensities.split(';') if row_intensities else []
        scores = row_qscores.split(';') if row_qscores else []
        for i, mass in enumerate(row_masses.split(';')):
            try:
                intensity = sum(float(v) for v in groups[i].split(',') if v) if i < len(groups) else 1.0
                mass = float(mass)
            except ValueError:
                continue
            try:
                qscore = float(scores[i]) if i < len(scores) else np.nan
            except ValueError:
                qscore = np.nan
            rt.append(row_rt)
            masses.append(mass)
            intensities.append(intensity)
            qscores.append(qscore)
    return pd.DataFrame({'rt': rt, 'mass': masses, 'intensity': intensities, 'qscore': qscores})


def _count(values: pd.Series) -> pd.Series:
    """Entries of a ';'-separated list per row (empty entries left out); NaN for no list."""
    return values.map(
        lambda v: float(sum(1 for part in v.split(';') if part.strip()))
        if isinstance(v, str) and v.strip() else np.nan
    )


#: Largest mass difference (Da) between a precursor and the survey mass it was
#: chosen from; FLASHIda may pick a neighbouring isotope (±1 Da).
MATCH_TOLERANCE = 1.1


def _matching_intensity(survey: pd.DataFrame, picked: pd.DataFrame) -> np.ndarray:
    """Intensity of the survey mass each precursor was chosen from: the nearest
    mass of the same survey scan (retention time) within MATCH_TOLERANCE; NaN
    where there is none."""
    if not len(survey) or not len(picked):
        return np.full(len(picked), np.nan)
    left = pd.DataFrame({
        'order': np.arange(len(picked)),
        'rt_key': np.round(picked['RT'].to_numpy(float), 6),
        'mass': picked['Mass'].to_numpy(float),
    }).sort_values('mass')
    right = pd.DataFrame({
        'rt_key': np.round(survey['rt'].to_numpy(float), 6),
        'mass': survey['mass'].to_numpy(float),
        'survey_intensity': survey['intensity'].to_numpy(float),
    }).sort_values('mass')
    matched = pd.merge_asof(left, right, on='mass', by='rt_key', direction='nearest',
                            tolerance=MATCH_TOLERANCE)
    return matched.sort_values('order')['survey_intensity'].to_numpy(float)


def decisions(tables: dict) -> dict:
    """What FLASHIda decided in a run.

    - ``scans``: one row per non-survey scan it commanded, with what it was for
      and what came of it;
    - ``precursors``: one row per precursor it chose, with its scans, the
      outcome of explorations and quantification and the proteoform identified;
    - ``points``: the masses of the survey scans and the chosen precursors, for
      the mass/retention time map;
    - ``summary``: counts.
    """
    commands = tables.get('scan_commands', pd.DataFrame())
    results = tables.get('scan_results', pd.DataFrame())
    ident = tables.get('identification', pd.DataFrame())
    pooled = tables.get('pooled_identification', pd.DataFrame())

    rt_of = dict(zip(_col(results, 'tracking_id'), _num(results, 'rt')))
    result_of = results.set_index('tracking_id') if 'tracking_id' in results.columns else None

    surveys = commands[_col(commands, 'scan_type') == 'survey']
    chosen = commands[_col(commands, 'scan_type') != 'survey'].copy()

    # Precursors in the order FLASHIda chose them; `precursor` is the number
    # every view of the viewer selects by
    order = pd.unique(_col(chosen, 'precursor_id'))
    number = {pid: i + 1 for i, pid in enumerate(order)}

    scans = pd.DataFrame({
        'scan': range(len(chosen)),
        'Scan': _col(chosen, 'tracking_id').to_numpy(),
        'precursor': _col(chosen, 'precursor_id').map(number).fillna(0).astype('int64').to_numpy(),
        'Decision': _col(chosen, 'scan_type').str.capitalize().to_numpy(),
        'MS': pd.to_numeric(_col(chosen, 'ms_level'), errors='coerce').to_numpy(),
        'Parent': _col(chosen, 'parent_tracking_id').to_numpy(),
        'Mass': _last(chosen, 'mono_mass').to_numpy(),
        'z': _last(chosen, 'charge').to_numpy(),
        'm/z': _last(chosen, 'precursor_mz').to_numpy(),
        'Width': _last(chosen, 'isolation_width').to_numpy(),
        'QScore': _last(chosen, 'qscore').to_numpy(),
        'Activation': _col(chosen, 'activation').str.split(';').str[-1].to_numpy(),
        'Energy': _last(chosen, 'collision_energy').to_numpy(),
        'Fragment': (_col(chosen, 'ion_type') + _col(chosen, 'ion_index').where(
            _col(chosen, 'ion_type') != '', '')).to_numpy(),
        'Priority': pd.to_numeric(_col(chosen, 'priority'), errors='coerce').to_numpy(),
    })
    acquired = scans['Scan'].isin(rt_of)
    scans['Acquired'] = acquired.to_numpy()
    # The instrument acquires commanded scans in the order they were sent, so
    # a scan sent before the last one that came back will not come back any
    # more (dropped, or the run stopped): it is not waiting
    sent = _col(commands, 'tracking_id')
    came_back = np.flatnonzero(sent.isin(rt_of).to_numpy())
    last_back = came_back[-1] if len(came_back) else -1
    position = np.flatnonzero((_col(commands, 'scan_type') != 'survey').to_numpy())
    scans['Waiting'] = (~acquired.to_numpy()) & (position > last_back)
    scans['RT'] = scans['Scan'].map(rt_of)
    scans['Survey RT'] = scans['Parent'].map(rt_of)

    def result(column, numeric=True):
        if result_of is None or column not in result_of.columns:
            return pd.Series(np.nan if numeric else '', index=scans.index)
        values = scans['Scan'].map(result_of[column]).fillna('')
        return pd.to_numeric(values, errors='coerce') if numeric else values

    scans['Masses'] = result('mass_count')
    scans['Tags'] = result('tag_count').where(lambda s: s >= 0)
    scans['Variant'] = [
        f'{int(v) + 1}/{int(n)}' if v >= 0 and n > 0 else ''
        for v, n in zip(result('variant_index').fillna(-1), result('total_variants').fillna(0))
    ]
    scans['Exploration score'] = result('exploration_score').where(lambda s: s >= 0)
    winners = set(result('winner_tracking_id', numeric=False)) - {''}
    scans['Winner'] = scans['Scan'].isin(winners)
    scans['Quant verdict'] = result('group_verdict', numeric=False)

    if 'tracking_id' in ident.columns and len(ident):
        best = ident.assign(_score=_num(ident, 'flash_extender_score')).sort_values(
            '_score', ascending=False).drop_duplicates('tracking_id').set_index('tracking_id')
        scans['Proteoform'] = scans['Scan'].map(best['proteoform']).fillna('')
        scans['Score'] = scans['Scan'].map(best['_score']).where(lambda s: s >= 0)
        # Fragment ions of the identification, listed ';'-separated by FLASHIda
        scans['Fragment ions'] = scans['Scan'].map(_count(_col(best, 'ms2_fragments')))
    else:
        scans['Proteoform'] = ''
        scans['Score'] = np.nan
        scans['Fragment ions'] = np.nan
    scans['Status'] = np.select(
        [scans['Acquired'], scans['Waiting']], ['acquired', 'waiting'], 'not acquired')

    # One row per precursor
    groups = scans.groupby('precursor', sort=True)
    precursors = pd.DataFrame({
        'precursor': list(groups.groups),
        'RT': groups['Survey RT'].min().to_numpy(),
        'Mass': groups.apply(lambda g: g.loc[g['MS'] == g['MS'].min(), 'Mass'].iloc[0],
                             include_groups=False).to_numpy() if len(scans) else [],
        'z': groups.apply(lambda g: g.loc[g['MS'] == g['MS'].min(), 'z'].iloc[0],
                          include_groups=False).to_numpy() if len(scans) else [],
        'QScore': groups['QScore'].max().to_numpy(),
        'MS2': groups['MS'].apply(lambda s: int((s == 2).sum())).to_numpy(),
        'MS3': groups['MS'].apply(lambda s: int((s == 3).sum())).to_numpy(),
        'Waiting': groups['Waiting'].apply(lambda s: int(s.sum())).to_numpy(),
        'Decisions': groups['Decision'].apply(lambda s: ', '.join(pd.unique(s))).to_numpy(),
        'Best energy': groups.apply(
            lambda g: g.loc[g['Winner'], 'Energy'].iloc[0] if g['Winner'].any() else np.nan,
            include_groups=False).to_numpy() if len(scans) else [],
        'Quant verdict': groups['Quant verdict'].apply(
            lambda s: next((v for v in s if v), '')).to_numpy(),
        # Most masses deconvolved from one of its MS2 scans, and most fragment
        # ions one of its identifications matched (NaN: none acquired/identified)
        'MS2 masses': groups.apply(
            lambda g: g.loc[g['MS'] == 2, 'Masses'].max(), include_groups=False
        ).to_numpy() if len(scans) else [],
        'Fragment ions': groups['Fragment ions'].max().to_numpy(),
    })
    precursors['Proteoform'] = ''
    precursors['Score'] = np.nan
    precursors['Coverage'] = np.nan
    if 'precursor_id' in pooled.columns and len(pooled):
        latest = pooled.assign(_u=_num(pooled, 'update_index')).sort_values('_u').drop_duplicates(
            'precursor_id', keep='last')
        latest = latest.assign(precursor=latest['precursor_id'].map(number)).dropna(
            subset=['precursor']).set_index('precursor')
        precursors['Proteoform'] = precursors['precursor'].map(latest['proteoform']).fillna('')
        precursors['Score'] = precursors['precursor'].map(_num(latest, 'flash_extender_score'))
        precursors['Coverage'] = precursors['precursor'].map(_num(latest, 'coverage_pct') * 100)
    # A precursor identified in a single scan but not pooled
    single = scans[scans['Proteoform'] != ''].sort_values('Score', ascending=False).drop_duplicates(
        'precursor').set_index('precursor')
    missing = precursors['Proteoform'] == ''
    precursors.loc[missing, 'Proteoform'] = precursors.loc[missing, 'precursor'].map(
        single['Proteoform']).fillna('')
    precursors.loc[missing, 'Score'] = precursors.loc[missing, 'precursor'].map(single['Score'])
    precursors['Identified'] = precursors['Proteoform'] != ''

    # Mass/RT map: survey masses, and the precursors on top. A precursor keeps
    # the largest intensity, so downsampling never drops it; `survey_intensity`
    # is its own (that of the survey mass it was chosen from).
    survey = survey_masses(results)
    survey['precursor'] = 0
    survey['Decision'] = 'Deconvolved'
    survey['Layer'] = 'Survey masses'
    survey['survey_intensity'] = survey['intensity']
    survey['ms2_masses'] = np.nan
    survey['fragment_ions'] = np.nan
    picked = precursors.dropna(subset=['RT', 'Mass'])
    top = pd.DataFrame({
        'rt': picked['RT'], 'mass': picked['Mass'],
        'intensity': survey['intensity'].max() if len(survey) else 1.0,
        'qscore': picked['QScore'],
        'precursor': picked['precursor'],
        'Decision': np.where(picked['Identified'], 'Identified', 'Chosen'),
        'Layer': np.where(picked['Identified'], 'Identified precursors', 'Chosen precursors'),
        'survey_intensity': _matching_intensity(survey, picked),
        'ms2_masses': picked['MS2 masses'],
        'fragment_ions': picked['Fragment ions'],
    })
    points = pd.concat([survey, top], ignore_index=True)
    # A precursor's own intensity goes in its label: on the map it has the largest
    points['label'] = [
        f'Precursor {p}: {m:.2f} Da'
        + (f'<br>Survey intensity: {si:.2e}' if np.isfinite(si) else '') if p else f'{m:.2f} Da'
        for p, m, si in zip(points['precursor'], points['mass'], points['survey_intensity'])
    ]

    summary = {
        'surveys': int(len(surveys)),
        'surveys_acquired': int(_col(surveys, 'tracking_id').isin(rt_of).sum()),
        'ms2': int((scans['MS'] == 2).sum()),
        'ms3': int((scans['MS'] == 3).sum()),
        'waiting': int(scans['Waiting'].sum()),
        'not_acquired': int((~scans['Acquired'] & ~scans['Waiting']).sum()),
        'precursors': int(len(precursors)),
        'identified': int(precursors['Identified'].sum()),
        'rt': float(np.nanmax(_num(results, 'rt'))) if len(results) and _num(results, 'rt').notna().any() else None,
    }
    return {'scans': scans, 'precursors': precursors, 'points': points, 'summary': summary}
