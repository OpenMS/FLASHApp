"""The FLASHIda Method editor's layout and its sense checks (Streamlit-free).

A method is edited by tab, one tab per step of an acquisition (`TAB_PARTS`):
Run, MS1, Deconvolution, Precursor selection, MS2, Identification, MS3 and
Quantification. A tab holds one or more parts (`PARTS`): the step itself, the
scan it acquires, and options such as FAIMS or an energy sweep. A part has
groups of key paths; every key of FLASHIda's method schema has one place.

Optional parts are switched on and off in the method (`present`, `switch`),
the way FLASHIda reads it: an omitted `faims.cv_values` means FAIMS on, the
`identification` block's presence means identification, and so on.

`inactive(method, path)` says why FLASHIda would not read a key with the method
as it is (the switch it hangs on is off, another option is chosen), and
`hints(method)` what it would refuse or what does not fit together, per part.
Both follow the descriptions in FLASHIda's method schema. `line(method, part)`
sums a part up in one line, `default_of(method, path)` names what FLASHIda uses
when a key is left out.
"""
import copy
import json
import os.path
from functools import lru_cache

from src import flashida_method as fm

# --- names --------------------------------------------------------------------

LABELS = {
    'global.method_name': 'Method name',
    'global.method_description': 'Description',
    'global.duration': 'Duration',
    'runtime.log_dir': 'Log folder',
    'deconvolution.tqscore_threshold': 'Qscore that closes a species',
    'deconvolution.score_threshold': 'Lowest qscore',
    'deconvolution.min_charge': 'Lowest charge',
    'deconvolution.max_charge': 'Highest charge',
    'deconvolution.min_mass': 'Lowest mass',
    'deconvolution.max_mass': 'Highest mass',
    'deconvolution.tol': 'Tolerance MS1, MS2, MS3',
    'precursor_selection.max_precursors': 'Precursors per survey',
    'precursor_selection.rank_by': 'Rank species by',
    'precursor_selection.consider_all_charges': 'Rank by the best charge',
    'precursor_selection.tie_threshold': 'Tie threshold',
    'precursor_selection.charge_selection': 'Isolate the charge with the best',
    'precursor_selection.precursor_charges': 'Charge states per species',
    'precursor_selection.min_precursor_charge': 'Lowest anchor charge',
    'precursor_selection.snr_threshold': 'Lowest S/N per charge',
    'precursor_selection.min_charge_states': 'Fewest charge states',
    'precursor_selection.max_charge_states': 'Most charge states',
    'precursor_selection.rt_window': 'Exclusion window',
    'precursor_selection.targeting': 'Target list',
    'precursor_selection.strict_inclusion': 'Inclusion list only',
    'precursor_selection.log_target_rt_window': 'Target logs: keep their RT window',
    'precursor_selection.log_target_charges': 'Target logs: charges',
    'precursor_selection.additional_scans': 'Acquire for every precursor',
    'files.inclusion_list': 'Inclusion list',
    'files.target_logs': 'Target logs',
    'identification.reference.fasta': 'Protein database (FASTA)',
    'identification.reference.sequence': 'Protein sequence',
    'identification.method': 'Search',
    'identification.threshold': 'Identified at score',
    'identification.search_threads': 'Search threads',
    'identification.flashtnt.min_length': 'Shortest tag',
    'identification.flashtnt.max_length': 'Longest tag',
    'identification.flashtnt.allow_gap': 'Allow mass gaps in tags',
    'identification.flashtnt.max_aa_in_gap': 'Residues per gap',
    'identification.flashtnt.gapped_tags': 'Gapped tags',
    'identification.flashtnt.max_protein_hits': 'Candidate proteins per MS2',
    'identification.flashtnt.terminal_ladders': 'Terminal ladders',
    'identification.flashtnt.precursor_ladders': 'Precursor ladders',
    'identification.flashtnt.skip_precursor_inference': 'Skip precursor inference',
    'identification.flashtnt.var_mod': 'Variable modifications',
    'identification.flashtnt.max_var_mod_count': 'Variable modifications per proteoform',
    'identification.flashtnt.max_blind_mod_count': 'Unknown modifications',
    'identification.flashtnt.max_mod_mass': 'Largest unknown modification',
    'identification.flashtnt.fixed_mod': 'Fixed modifications',
    'identification.refine.fdr': 'FDR',
    'identification.refine.shuffles': 'Shuffled sequences',
    'identification.refine.grid_da': 'Shift grid step',
    'identification.refine.max_shift_da': 'Largest shift',
    'identification.refine.min_ions': 'Fewest matched ions',
    'identification.refine.max_levels': 'Most shift levels',
    'identification.refine.jump_penalty': 'Jump penalty',
    'identification.characterization.mode': 'MS3 objective',
    'identification.characterization.max_targets': 'Targets per precursor',
    'identification.characterization.per_proteoform': 'One plan per proteoform',
    'identification.characterization.fragment_charges': 'Charge states per fragment',
    'identification.characterization.min_fragment_charge': 'Lowest fragment charge',
    'identification.characterization.min_target_mass': 'Lowest target mass',
    'quantification.enabled': 'Quantify',
    'quantification.labelling': 'Labels',
    'quantification.scan': 'Reporter ions from',
    'quantification.ms3_targets': 'Quant MS3s per identification',
    'quantification.reporter_mz_tol': 'Reporter m/z tolerance',
    'quantification.correction_matrix': 'Isotope correction matrix',
    'quantification.conditions': 'Conditions',
    'quantification.channel_normalization': 'Channel normalization',
    'quantification.test': 'Test',
    'quantification.fold_change_threshold': 'Fold change',
    'quantification.p_value_threshold': 'p-value',
    'quantification.empty_condition': 'An empty condition counts as',
    'quantification.identify': 'Identify the groups that are',
    'quantification.enriched_in': 'Enriched in',
    'quantification.scan_normalize': 'Normalize each scan',
    'quantification.aggregate': 'Combine scans by',
    'quantification.min_usable_scans': 'Fewest usable scans',
    'faims.cv_values': 'FAIMS CVs',
    'faims.max_cv_skip': 'Most cycles a CV is skipped',
    'faims.cv_precursor_threshold': 'Precursors that keep a CV',
    'scheduling.agc_interval_seconds': 'AGC prescan every',
    'scheduling.target_depth': 'Commands in flight',
    'scheduling.cycle_time.enabled': 'Limit the gap between surveys',
    'scheduling.cycle_time.value_ms': 'Longest gap between surveys',
    'scheduling.scan_timeout.enabled': 'Discard commands that wait too long',
    'scheduling.scan_timeout.value_ms': 'Longest wait',
    'ms_settings.additional_ms2': 'Named MS2 scans',
}

# Keys every scan config has, and every sweep level, named by their last part
SCAN_LABELS = {
    'analyzer': 'Analyzer', 'resolution': 'Resolution', 'data_type': 'Data type',
    'scan_rate': 'Ion trap scan rate', 'microscans': 'Microscans',
    'agc_target': 'AGC target', 'max_it': 'Max. injection time',
    'first_mass': 'First m/z', 'last_mass': 'Last m/z',
    'activation': 'Activation', 'collision_energy': 'Collision energy',
    'reaction_time': 'Reaction time', 'reagent_agc_target': 'Reagent AGC target',
    'reagent_max_it': 'Reagent max. injection time',
    'rf_lens': 'RF lens', 'source_cid': 'In-source CID', 'source_cid_scaling': 'Source CID scaling',
}
SWEEP_LABELS = {
    'metric': 'Score pre-scans by', 'activations': 'Activations to sweep',
    'tolerance_ppm': 'Pre-scan tolerance', 'remaining_precursor_target': 'Remaining precursor target',
    'ce_min': 'Lowest energy', 'ce_max': 'Highest energy', 'ce_step': 'Energy step',
    'reaction_time_min': 'Shortest reaction time', 'reaction_time_max': 'Longest reaction time',
    'reaction_time_step': 'Reaction time step',
    'monitor_ms1.enabled': 'MS1 monitor scans while sweeping', 'monitor_ms1.interval_ms': 'Monitor scan every',
    'overrides': 'Pre-scan settings', 'kinetic_model.k2': 'Kinetic model k2',
}

SCANS = {'MS1': 'ms1', 'MS2': 'ms2', 'MS3': 'ms3', 'MS2 quant': 'ms2_quant', 'MS3 quant': 'ms3_quant'}
SWEEPS = {'MS2 sweep': 'precursor_selection.exploration',
          'MS3 sweep': 'identification.characterization.exploration'}


def _context(path: str):
    """(context name, rest of the path) for a scan or sweep key, else None."""
    for name, scan in SCANS.items():
        if path.startswith(f'ms_settings.{scan}.'):
            return name, path[len(f'ms_settings.{scan}.'):], SCAN_LABELS
    for name, prefix in SWEEPS.items():
        if path.startswith(prefix + '.'):
            return name, path[len(prefix) + 1:], SWEEP_LABELS
    return None


def name(path: str, short: bool = False) -> str:
    """The key's name in the dialog. A scan or sweep key names its scan or
    level ('MS2 quant: AGC target') unless `short`, for the views that show
    one scan or level at a time."""
    if path in LABELS:
        return LABELS[path]
    found = _context(path)
    if found:
        context, rest, labels = found
        text = labels.get(rest) or _plain(rest)
        return text if short else f'{context}: {text}'
    return _plain(path.split('.')[-1])


def _plain(key: str) -> str:
    text = key.replace('_', ' ')
    return text[:1].upper() + text[1:]


# --- where each key sits ------------------------------------------------------------

def _scan_view(scan: str) -> dict:
    p = f'ms_settings.{scan}.'
    groups = {
        'Readout': ['analyzer', 'resolution', 'data_type', 'scan_rate', 'microscans'],
        'Injection and range': ['agc_target', 'max_it', 'first_mass', 'last_mass'],
        'Fragmentation': ['activation', 'collision_energy', 'reaction_time',
                          'reagent_agc_target', 'reagent_max_it'],
        'Ion source': ['rf_lens', 'source_cid', 'source_cid_scaling'],
    }
    if scan == 'ms1':
        del groups['Fragmentation']
    return {g: [p + k for k in ks] for g, ks in groups.items()}


def _sweep_view(prefix: str) -> dict:
    p = prefix + '.'
    return {
        'Sweep': [p + 'metric', p + 'activations', p + 'tolerance_ppm', p + 'remaining_precursor_target'],
        'Collision energy grid': [p + 'ce_min', p + 'ce_max', p + 'ce_step'],
        'Reaction time grid': [p + 'reaction_time_min', p + 'reaction_time_max', p + 'reaction_time_step'],
        'MS1 monitor scans': [p + 'monitor_ms1.enabled', p + 'monitor_ms1.interval_ms'],
        'Pre-scans and kinetic model': [p + 'overrides', p + 'kinetic_model.k2'],
    }


TABS = {
    'General': {None: {
        'Method': ['global.method_name', 'global.duration', 'global.method_description', 'runtime.log_dir'],
    }},
    'Deconvolution': {None: {
        'Charges and masses': ['deconvolution.min_charge', 'deconvolution.max_charge',
                               'deconvolution.min_mass', 'deconvolution.max_mass'],
        'Tolerance and qscores': ['deconvolution.tol', 'deconvolution.score_threshold',
                                  'deconvolution.tqscore_threshold'],
    }},
    'Precursor selection': {None: {
        'Per survey': ['precursor_selection.max_precursors', 'precursor_selection.rank_by',
                       'precursor_selection.consider_all_charges', 'precursor_selection.tie_threshold'],
        'Charge states': ['precursor_selection.charge_selection', 'precursor_selection.precursor_charges',
                          'precursor_selection.min_precursor_charge', 'precursor_selection.snr_threshold',
                          'precursor_selection.min_charge_states', 'precursor_selection.max_charge_states'],
        'Exclusion and target lists': ['precursor_selection.rt_window', 'precursor_selection.targeting',
                                       'files.inclusion_list', 'precursor_selection.strict_inclusion',
                                       'files.target_logs', 'precursor_selection.log_target_rt_window',
                                       'precursor_selection.log_target_charges'],
    }},
    'Identification': {None: {
        'Reference and search': ['identification.reference.fasta', 'identification.reference.sequence',
                                 'identification.method', 'identification.threshold',
                                 'identification.search_threads'],
        'FLASHTnT tags and candidates': [
            'identification.flashtnt.min_length', 'identification.flashtnt.max_length',
            'identification.flashtnt.allow_gap', 'identification.flashtnt.max_aa_in_gap',
            'identification.flashtnt.gapped_tags', 'identification.flashtnt.max_protein_hits',
            'identification.flashtnt.terminal_ladders', 'identification.flashtnt.precursor_ladders',
            'identification.flashtnt.skip_precursor_inference'],
        'FLASHTnT modifications': [
            'identification.flashtnt.var_mod', 'identification.flashtnt.max_var_mod_count',
            'identification.flashtnt.max_blind_mod_count', 'identification.flashtnt.max_mod_mass',
            'identification.flashtnt.fixed_mod'],
        'Refine': ['identification.refine.fdr', 'identification.refine.shuffles',
                   'identification.refine.grid_da', 'identification.refine.max_shift_da',
                   'identification.refine.min_ions', 'identification.refine.max_levels',
                   'identification.refine.jump_penalty'],
    }},
    'MS3': {None: {
        'Characterization': ['identification.characterization.mode',
                             'identification.characterization.max_targets',
                             'identification.characterization.per_proteoform',
                             'identification.characterization.fragment_charges',
                             'identification.characterization.min_fragment_charge',
                             'identification.characterization.min_target_mass'],
    }},
    'Quantification': {None: {
        'Labels and reporter ions': ['quantification.enabled', 'quantification.labelling',
                                     'quantification.scan', 'quantification.ms3_targets',
                                     'quantification.reporter_mz_tol', 'quantification.correction_matrix'],
        'Conditions': ['quantification.conditions', 'quantification.channel_normalization'],
        'Differential test': ['quantification.test', 'quantification.fold_change_threshold',
                              'quantification.p_value_threshold', 'quantification.empty_condition',
                              'quantification.identify', 'quantification.enriched_in'],
        'Combining scans': ['quantification.scan_normalize', 'quantification.aggregate',
                            'quantification.min_usable_scans'],
    }},
    'Energy sweep': {level: _sweep_view(prefix) for level, prefix in SWEEPS.items()},
    'Scans': {**{label: _scan_view(scan) for label, scan in SCANS.items()},
              'Additional MS2': {'Named MS2 scans': ['ms_settings.additional_ms2',
                                                     'precursor_selection.additional_scans']}},
    'Instrument': {None: {
        'FAIMS': ['faims.cv_values', 'faims.max_cv_skip', 'faims.cv_precursor_threshold'],
        'Scheduling': ['scheduling.agc_interval_seconds', 'scheduling.target_depth',
                       'scheduling.cycle_time.enabled', 'scheduling.cycle_time.value_ms',
                       'scheduling.scan_timeout.enabled', 'scheduling.scan_timeout.value_ms'],
    }},
}


def placed() -> list:
    """Every key path the tabs show, in tab order."""
    return [p for views in TABS.values() for groups in views.values()
            for paths in groups.values() for p in paths]


# Turning a choice off removes this block, which is how FLASHIda is told it is off
CLEARS = {'identification.characterization.mode': 'identification.characterization'}
# What the empty choice of a key means, where it is not the instrument method's value
EMPTY_CHOICE = {'identification.characterization.mode': 'Off (no MS3)'}
# Choices for keys the schema leaves free text, though FLASHIda knows a few values
ACTIVATIONS = ['HCD', 'CID', 'ETD', 'EThcD', 'ETciD']


DATA_TYPES = ['Centroid', 'Profile']
SCAN_RATES = ['Normal', 'Enhanced', 'Rapid', 'Turbo']  # ion trap scan rates FLASHIda passes on
# Resolutions an Orbitrap Tribrid (Fusion, Lumos, Eclipse: what FLASHIda drives) offers, at m/z 200
ORBITRAP_RESOLUTIONS = [7500, 15000, 30000, 50000, 60000, 120000, 240000, 500000]


def choices(path: str, method=None):
    """The values a key with a fixed set of values (the schema lists them in its
    description only) can take, or None. For a list key, its items."""
    if path.startswith('ms_settings.'):
        leaf = path.rsplit('.', 1)[-1]
        if leaf == 'activation':
            return ACTIVATIONS
        if leaf == 'data_type':
            return DATA_TYPES
        if leaf == 'scan_rate':
            return SCAN_RATES
        if leaf == 'resolution':
            return ORBITRAP_RESOLUTIONS
    if path.endswith('.exploration.activations'):
        return ACTIVATIONS
    if path == 'precursor_selection.additional_scans':
        return sorted((fm.get(method or {}, 'ms_settings.additional_ms2') or {}).keys())
    if path == 'quantification.enriched_in':
        names = [c.get('name') for c in fm.get(method or {}, 'quantification.conditions') or []
                 if isinstance(c, dict) and c.get('name')]
        return ['either'] + names
    return None



# --- when FLASHIda does not read a key ------------------------------------------------

_MISSING = object()


@lru_cache(maxsize=1)
def _keys() -> dict:
    return fm.keys()


def value(method: dict, path: str):
    """The key's value in the method, else its default."""
    found = fm.get(method, path, missing=_MISSING)
    if found is _MISSING:
        key = _keys().get(path)
        return key.default if key is not None else None
    return found


def _sweep_rules(prefix: str, level: str) -> list:
    metric = f'{prefix}.metric'
    return [
        ([prefix + '.'], metric, lambda m: value(m, metric) not in (None, 'none'),
         f'Read only when the {level} sweep scores its pre-scans (Score pre-scans by is not none).'),
        ([f'{prefix}.remaining_precursor_target'], metric,
         lambda m: value(m, metric) == 'remaining_precursor',
         'Read only when pre-scans are scored by remaining_precursor.'),
        ([f'{prefix}.kinetic_model.'], metric, lambda m: value(m, metric) == 'kinetic_model',
         'Read only when pre-scans are scored by kinetic_model.'),
        ([f'{prefix}.monitor_ms1.interval_ms'], f'{prefix}.monitor_ms1.enabled',
         lambda m: bool(value(m, f'{prefix}.monitor_ms1.enabled')),
         'Read only with MS1 monitor scans on.'),
    ]


# (paths or prefixes ending in '.', the key that decides, active?, why not)
RULES = [
    (['precursor_selection.strict_inclusion', 'files.inclusion_list',
      'precursor_selection.log_target_rt_window', 'precursor_selection.log_target_charges'],
     'precursor_selection.targeting',
     lambda m: value(m, 'precursor_selection.targeting') == 'inclusion',
     'Read only when Target list is inclusion.'),
    (['faims.max_cv_skip', 'faims.cv_precursor_threshold'], 'faims.cv_values',
     lambda m: len(value(m, 'faims.cv_values') or []) >= 2,
     'Read only with two or more FAIMS CVs.'),
    (['identification.flashtnt.'], 'identification.method',
     lambda m: value(m, 'identification.method') == 'flashtnt',
     'Read only when Search is flashtnt.'),
    (['identification.refine.'], 'identification.method',
     lambda m: value(m, 'identification.method') == 'refine',
     'Read only when Search is refine.'),
    (['identification.flashtnt.max_aa_in_gap'], 'identification.flashtnt.allow_gap',
     lambda m: bool(value(m, 'identification.flashtnt.allow_gap')),
     'Read only when mass gaps are allowed in tags.'),
    (['identification.characterization.'], 'identification.characterization.mode',
     lambda m: value(m, 'identification.characterization.mode') is not None,
     'Read only with MS3 on: pick an MS3 objective in the MS3 tab.'),
    *_sweep_rules('precursor_selection.exploration', 'MS2'),
    *_sweep_rules('identification.characterization.exploration', 'MS3'),
    (['quantification.'], 'quantification.enabled',
     lambda m: bool(value(m, 'quantification.enabled')),
     'Read only with Quantify on.'),
    (['quantification.ms3_targets'], 'quantification.scan',
     lambda m: value(m, 'quantification.scan') in ('ms3', 'both'),
     'Read only when reporter ions come from ms3 or both.'),
    (['quantification.enriched_in'], 'quantification.identify',
     lambda m: value(m, 'quantification.identify') == 'differential',
     'Read only when differential groups are identified.'),
    (['quantification.p_value_threshold'], 'quantification.test',
     lambda m: 'welch' in str(value(m, 'quantification.test')),
     'Read only by a Welch test.'),
    (['quantification.fold_change_threshold'], 'quantification.test',
     lambda m: 'fold_change' in str(value(m, 'quantification.test')),
     'Read only by a fold change test.'),
    (['scheduling.cycle_time.value_ms'], 'scheduling.cycle_time.enabled',
     lambda m: bool(value(m, 'scheduling.cycle_time.enabled')),
     'Read only when the gap between surveys is limited.'),
    (['scheduling.scan_timeout.value_ms'], 'scheduling.scan_timeout.enabled',
     lambda m: bool(value(m, 'scheduling.scan_timeout.enabled')),
     'Read only when commands that wait too long are discarded.'),
]

# When FLASHIda sends a scan config; shown above that scan's settings
SCAN_USE = {
    'ms1': 'The survey scan.',
    'ms2': 'The MS2 of every selected precursor (with Quantify on, the identification scan '
           'a quantified group buys).',
    'ms3': 'The characterization MS3; sent only with MS3 on.',
    'ms2_quant': 'The reporter ion MS2; sent only with Quantify on.',
    'ms3_quant': 'The reporter ion MS3; sent only with Quantify on and reporter ions from ms3 or both.',
}


def inactive(method: dict, path: str):
    """Why FLASHIda would not read this key with the method as it is, else None."""
    for paths, decider, active, why in RULES:
        if path == decider:
            continue
        if any(path == p or (p.endswith('.') and path.startswith(p)) for p in paths):
            if not active(method):
                return why
    return None


RULES.append((['files.target_logs'], 'precursor_selection.targeting',
              lambda m: value(m, 'precursor_selection.targeting') != 'none',
              'Read only when Target list is not none.'))


# --- keys that depend on another key (from the schema's descriptions) ----------------

CE_ACTIVATIONS = ('HCD', 'CID', 'EThcD', 'ETciD')
ETD_ACTIVATIONS = ('ETD', 'EThcD', 'ETciD')
_SCAN_NAMES = {'ms1': 'MS1', 'ms2': 'MS2', 'ms3': 'MS3', 'ms2_quant': 'Quant MS2', 'ms3_quant': 'Quant MS3'}


def _scan_rules(scan: str) -> list:
    p = f'ms_settings.{scan}.'
    analyzer, activation = p + 'analyzer', p + 'activation'
    found = [
        ([p + 'scan_rate'], analyzer, lambda m: value(m, analyzer) != 'Orbitrap',
         'Ion trap scan rate: not read while Analyzer is Orbitrap.'),
        ([p + 'resolution'], analyzer, lambda m: value(m, analyzer) != 'IonTrap',
         'Orbitrap resolution: not read while Analyzer is IonTrap.'),
        ([p + 'last_mass'], p + 'first_mass', lambda m: (value(m, p + 'first_mass') or 0) > 0,
         'A scan range needs First m/z as well; without it only a first mass is commanded.'),
    ]
    if scan != 'ms1':
        found += [
            ([p + 'collision_energy'], activation,
             lambda m: value(m, activation) in (None, '') or value(m, activation) in CE_ACTIVATIONS,
             'Read only by HCD, CID, EThcD and ETciD.'),
            ([p + 'reaction_time', p + 'reagent_max_it', p + 'reagent_agc_target'], activation,
             lambda m: value(m, activation) in (None, '') or value(m, activation) in ETD_ACTIVATIONS,
             'Read only by ETD, EThcD and ETciD.'),
        ]
    return found


def _swept(m, prefix: str, scan: str) -> list:
    acts = value(m, f'{prefix}.activations')
    return list(acts) if acts else [value(m, f'ms_settings.{scan}.activation')]


def _grid_rules(prefix: str, scan: str) -> list:
    return [
        ([f'{prefix}.ce_min', f'{prefix}.ce_max', f'{prefix}.ce_step'], f'{prefix}.activations',
         lambda m: any(a in CE_ACTIVATIONS or a in (None, '') for a in _swept(m, prefix, scan)),
         'Read only when a swept activation varies collision energy (HCD, CID, EThcD, ETciD).'),
        ([f'{prefix}.reaction_time_min', f'{prefix}.reaction_time_max', f'{prefix}.reaction_time_step'],
         f'{prefix}.activations',
         lambda m: any(a in ETD_ACTIVATIONS or a in (None, '') for a in _swept(m, prefix, scan)),
         'Read only when an ETD-family activation is swept (ETD, EThcD, ETciD).'),
    ]


RULES.extend([
    *[r for scan in _SCAN_NAMES for r in _scan_rules(scan)],
    *_grid_rules('precursor_selection.exploration', 'ms2'),
    *_grid_rules('identification.characterization.exploration', 'ms3'),
    (['precursor_selection.consider_all_charges'], 'precursor_selection.rank_by',
     lambda m: value(m, 'precursor_selection.rank_by') == 'qscore',
     'Applies only when species are ranked by qscore.'),
    (['precursor_selection.tie_threshold'], 'files.inclusion_list',
     lambda m: value(m, 'precursor_selection.targeting') == 'inclusion' and bool(value(m, 'files.inclusion_list')),
     'Read only under Target list inclusion with an inclusion list loaded.'),
    (['precursor_selection.log_target_rt_window', 'precursor_selection.log_target_charges'], 'files.target_logs',
     lambda m: (value(m, 'precursor_selection.targeting') == 'inclusion' and bool(value(m, 'files.target_logs'))
                and not value(m, 'files.inclusion_list')),
     'Read only under Target list inclusion with target logs and no inclusion list.'),
    (['precursor_selection.max_charge_states'], 'precursor_selection.precursor_charges',
     lambda m: value(m, 'precursor_selection.precursor_charges') != 'single',
     'Not read while Charge states per species is single.'),
    (['identification.characterization.min_target_mass'], 'identification.characterization.mode',
     lambda m: value(m, 'identification.characterization.mode') == 'exhaustive',
     'Read only by the exhaustive MS3 objective.'),
    (['quantification.scan_normalize'], 'precursor_selection.precursor_charges',
     lambda m: value(m, 'precursor_selection.precursor_charges') == 'separate',
     'Matters only when charge states are acquired in separate scans (Precursor selection).'),
    (['faims.cv_precursor_threshold'], 'faims.max_cv_skip',
     lambda m: (value(m, 'faims.max_cv_skip') or 0) > 0,
     'Not read while Most cycles a CV is skipped is 0.'),
])


def etd_scans(method) -> list:
    """The scans and sweeps FLASHIda sends with an ETD-family activation (ETD,
    EThcD, ETciD), by name. ETD needs the reagent set up in Tune, which a
    method cannot do, so the Acquisition page warns about these methods."""
    if not isinstance(method, dict):
        return []
    sent = {
        'ms2': True,
        'ms3': value(method, 'identification.characterization.mode') is not None,
        'ms2_quant': bool(value(method, 'quantification.enabled')),
        'ms3_quant': (bool(value(method, 'quantification.enabled'))
                      and value(method, 'quantification.scan') in ('ms3', 'both')),
    }
    found = [_SCAN_NAMES[scan] for scan, on in sent.items()
             if on and value(method, f'ms_settings.{scan}.activation') in ETD_ACTIVATIONS]
    for label, prefix in SWEEPS.items():
        scan = 'ms3' if prefix.startswith('identification.') else 'ms2'
        if ((scan == 'ms2' or sent['ms3'])
                and value(method, f'{prefix}.metric') not in (None, 'none')
                and any(a in ETD_ACTIVATIONS for a in _swept(method, prefix, scan))):
            found.append(label)
    return found


_G = TABS


# --- the parts of an acquisition, one or more per tab ------------------------------------------------------

_Q = _G['Quantification'][None]
# Reporter ions from decides which quant scans are acquired, so it comes first
QUANT_GROUPS = {
    'Labels and reporter ions': ['quantification.scan', 'quantification.labelling',
                                 'quantification.reporter_mz_tol', 'quantification.ms3_targets',
                                 'quantification.correction_matrix'],
    **{g: ps for g, ps in _Q.items() if g != 'Labels and reporter ions'},
}
GLOBAL_GROUPS = {**_G['General'][None], 'Scheduling': _G['Instrument'][None]['Scheduling']}

# part: (title, icon, kind, settings groups, short names, what it does)
PARTS = {
    'faims': ('FAIMS', '🌀', 'option', {'FAIMS': _G['Instrument'][None]['FAIMS']}, False,
              'Compensation voltages. Each survey uses the next CV; its MS2 scans keep that CV.'),
    'ms1': ('MS1 survey scan', '🔬', 'scan', _scan_view('ms1'), True, SCAN_USE['ms1']),
    'deconv': ('Deconvolution', '🧮', 'step', _G['Deconvolution'][None], False,
               'Turns each survey into neutral masses with charges and qscores.'),
    'select': ('Precursor selection', '🎯', 'step', _G['Precursor selection'][None], False,
               'Picks the species of each survey that are acquired, and bars them afterwards.'),
    'ms2q': ('Quant MS2 scan', '🏷️', 'scan', _scan_view('ms2_quant'), True,
             'Settings `ms_settings.ms2_quant`. With quantification on, fired once per selected precursor; '
             'with Reporter ions from ms2 or both, its reporter ions decide whether the identification MS2 is bought.'),
    'quant': ('Quantification', '⚖️', 'step', QUANT_GROUPS, False,
              'Tests the reporter ions of each group and decides which groups buy the identification MS2.'),
    'sweep2': ('Energy sweep MS2', '📈', 'option', _sweep_view('precursor_selection.exploration'), True,
               'Replaces the MS2 of every selected precursor with pre-scans over an energy grid, '
               'then acquires at the best energy.'),
    'ms2': ('MS2 scan', '🔬', 'scan', _scan_view('ms2'), True, SCAN_USE['ms2']),
    'ms2x': ('Additional MS2 scans', '➕', 'option', _G['Scans']['Additional MS2'], False,
             'Named extra MS2 configs, acquired for every selected precursor after the first MS2.'),
    'ident': ('Identification', '🧬', 'step', _G['Identification'][None], False,
              'Searches every MS2 in real time and closes out identified species.'),
    'ms3q': ('Quant MS3 scan', '🏷️', 'scan', _scan_view('ms3_quant'), True,
             'Settings `ms_settings.ms3_quant`. Shown with Reporter ions from ms3 or both: reads reporter ions '
             'from labelled fragments of an identified species (needs identification).'),
    'ms3plan': ('MS3 characterization', '🧩', 'step', _G['MS3'][None], False,
                'Plans MS3 scans on fragments of identified proteoforms.'),
    'sweep3': ('Energy sweep MS3', '📈', 'option', _sweep_view('identification.characterization.exploration'),
               True, 'Sweeps every MS3 target over an energy grid before it is acquired.'),
    'ms3': ('MS3 scan', '🔬', 'scan', _scan_view('ms3'), True, SCAN_USE['ms3']),
}
REQUIRED = ['ms1', 'deconv', 'select', 'ms2']
# what a part needs, and what it cannot be combined with, before it can be added
NEEDS = {'ms3plan': ['ident'], 'sweep3': ['ms3plan'], 'ms3q': ['ident', 'quant']}
CONFLICTS = {'sweep2': ['ms2x'], 'ms2x': ['sweep2']}
# parts that come with another one and are not added on their own
COMES_WITH = {'ms2q': 'quant', 'ms3': 'ms3plan'}



def present(method) -> set:
    """The parts the method as it is has."""
    on = set(REQUIRED)
    if value(method, 'faims.cv_values'):
        on.add('faims')
    if value(method, 'precursor_selection.exploration.metric') not in (None, 'none'):
        on.add('sweep2')
    if fm.get(method, 'ms_settings.additional_ms2') or fm.get(method, 'precursor_selection.additional_scans'):
        on.add('ms2x')
    if 'identification' in method:
        on.add('ident')
        if value(method, 'identification.characterization.mode'):
            on |= {'ms3plan', 'ms3'}
            if value(method, 'identification.characterization.exploration.metric') not in (None, 'none'):
                on.add('sweep3')
    if value(method, 'quantification.enabled'):
        on |= {'quant', 'ms2q'}
        if value(method, 'quantification.scan') in ('ms3', 'both'):
            on.add('ms3q')
    return on


def switch(method, node: str, on: bool) -> None:
    """Switch a part of the acquisition on or off in the method."""
    if node == 'faims':  # an omitted key means FAIMS on at -50, so off is an empty list
        fm.put(method, 'faims.cv_values', [-50] if on else [])
    elif node == 'sweep2':
        if on:
            fm.put(method, 'precursor_selection.exploration.metric', 'mass_count')
        else:
            fm.drop(method, 'precursor_selection.exploration')
    elif node == 'ms2x':
        if on:
            fm.put(method, 'ms_settings.additional_ms2', {'second_hcd': {'activation': 'HCD', 'collision_energy': 35}})
            fm.put(method, 'precursor_selection.additional_scans', ['second_hcd'])
        else:
            fm.drop(method, 'ms_settings.additional_ms2')
            fm.drop(method, 'precursor_selection.additional_scans')
    elif node == 'ident':
        if on:
            fm.put(method, 'identification.method', 'flashtnt')
        else:
            fm.drop(method, 'identification')
            if value(method, 'quantification.scan') != 'ms2':
                fm.put(method, 'quantification.scan', 'ms2')
    elif node in ('ms3plan', 'ms3'):
        if on:
            fm.put(method, 'identification.characterization.mode', 'coverage')
        else:
            fm.drop(method, 'identification.characterization')
    elif node == 'sweep3':
        if on:
            fm.put(method, 'identification.characterization.exploration.metric', 'mass_count')
        else:
            fm.drop(method, 'identification.characterization.exploration')
    elif node in ('quant', 'ms2q'):
        fm.put(method, 'quantification.enabled', on)
    elif node == 'ms3q':
        fm.put(method, 'quantification.scan', 'both' if on else 'ms2')


def why_not_add(node: str, have: set):
    missing = [PARTS[m][0] for m in NEEDS.get(node, []) if m not in have]
    if missing:
        return 'Needs ' + ' and '.join(missing) + '.'
    clash = [PARTS[m][0] for m in CONFLICTS.get(node, []) if m in have]
    if clash:
        return 'Cannot be combined with ' + ', '.join(clash) + '.'
    return None


# --- does the method make sense: hints per part ---------------------------------------

OWNER = {p: n for n, spec in PARTS.items() for ps in spec[3].values() for p in ps}
OWNER.update({p: 'global' for ps in GLOBAL_GROUPS.values() for p in ps})


def owner(path: str) -> str:
    if path in OWNER:
        return OWNER[path]
    for p, n in OWNER.items():
        if p.startswith(path + '.') or path.startswith(p + '.'):
            return n
    return 'global'


def hints(method) -> dict:
    """{part: [problem, ...]}: what FLASHIda would refuse (schema) and what does not fit together."""
    have = present(method)
    found = {}

    def add(node, text):
        found.setdefault(node, []).append(text)

    for problem in fm.check(copy.deepcopy(method)):
        path = problem.split(':', 1)[0]
        add(owner(path), problem.split(':', 1)[1].strip() + f' (`{path}`)')
    for scan, node in (('ms2', 'ms2'), ('ms3', 'ms3'), ('ms2_quant', 'ms2q'), ('ms3_quant', 'ms3q')):
        if node in have and not value(method, f'ms_settings.{scan}.activation'):
            add(node, 'Pick an activation; FLASHIda refuses an MSn scan without one.')
    for scan, node in (('ms2', 'ms2'), ('ms3', 'ms3'), ('ms2_quant', 'ms2q'), ('ms3_quant', 'ms3q')):
        if node in have and value(method, f'ms_settings.{scan}.activation') in ETD_ACTIVATIONS \
                and not (value(method, f'ms_settings.{scan}.reaction_time') or 0) > 0:
            add(node, 'An ETD-family scan needs a Reaction time above 0; the instrument rejects 0.')
    if 'ident' in have:
        fasta, seq = value(method, 'identification.reference.fasta'), value(method, 'identification.reference.sequence')
        if not fasta and not seq:
            add('ident', 'Needs a protein database (FASTA) or a protein sequence.')
        if fasta and seq:
            add('ident', 'Set either a protein database or a sequence, not both.')
    if 'ms3plan' in have and value(method, 'identification.characterization.mode') == 'discriminating' \
            and value(method, 'identification.method') != 'refine':
        add('ms3plan', 'The discriminating objective needs Search: refine.')
    if 'quant' in have and not value(method, 'quantification.conditions'):
        add('quant', 'Needs the two conditions to compare.')
    if 'ms3q' in have:
        for ok, text in ((value(method, 'identification.threshold') is not None, 'an Identified at score'),
                         (value(method, 'quantification.labelling') == 'tmt6plex', 'labels tmt6plex'),
                         (value(method, 'identification.method') == 'flashtnt', 'Search: flashtnt')):
            if not ok:
                add('ms3q', f'Reporter ions from the MS3 need {text}.')
        if value(method, 'ms_settings.ms2.activation') == 'HCD':
            add('ms2', 'HCD on the identification MS2 releases the reporters before the reporter MS3; '
                       'use CID or ETD.')
    if 'sweep2' in have and 'ms2x' in have:
        add('sweep2', 'An MS2 sweep cannot be combined with additional MS2 scans.')
    for node, prefix in (('sweep2', 'precursor_selection.exploration'),
                         ('sweep3', 'identification.characterization.exploration')):
        if node not in have:
            continue
        metric = value(method, f'{prefix}.metric')
        if metric == 'fragment_count' and 'ident' not in have:
            add(node, 'Scoring by fragment_count needs Identification with a reference.')
        if metric == 'kinetic_model' and len(value(method, f'{prefix}.activations') or []) != 1:
            add(node, 'The kinetic model sweeps exactly one activation, HCD or CID.')
    targeting = value(method, 'precursor_selection.targeting')
    if targeting == 'inclusion' and not value(method, 'files.inclusion_list') \
            and not value(method, 'files.target_logs'):
        add('select', 'Inclusion targeting needs an inclusion list or target logs.')
    if targeting in ('in_depth', 'exclusion_masses') and not value(method, 'files.target_logs'):
        add('select', f'Targeting {targeting} reads target logs; none are set.')
    return found


# --- the line that sums up a part ----------------------------------------------

def _num(v, unit=''):
    if v in (None, '', 0):
        return '–'
    if isinstance(v, float) and v.is_integer():
        v = int(v)
    return f'{v:,}'.replace(',', ' ') + unit if isinstance(v, (int, float)) else f'{v}{unit}'


def _scan_line(method, scan):
    g = lambda k: value(method, f'ms_settings.{scan}.{k}')
    parts = []
    if scan != 'ms1':
        act = g('activation') or 'no activation'
        ce = g('collision_energy')
        parts.append(f'{act} {ce} NCE' if ce and act in ('HCD', 'CID', 'EThcD', 'ETciD') else act)
    parts.append(f"{g('analyzer') or 'instrument'} {('R ' + _num(g('resolution'))) if g('resolution') else ''}".strip())
    if g('max_it'):
        parts.append(f"{_num(g('max_it'))} ms")
    if g('first_mass') and g('last_mass'):
        parts.append(f"{_num(g('first_mass'))}–{_num(g('last_mass'))} m/z")
    return ' · '.join(parts)


def _sweep_line(method, prefix):
    g = lambda k: value(method, f'{prefix}.{k}')
    acts = ', '.join(g('activations') or []) or 'scan activation'
    return f"{g('metric')} · {acts} · CE {_num(g('ce_min'))}–{_num(g('ce_max'))} / {_num(g('ce_step'))}"


def line(method, node: str) -> str:
    v = lambda p: value(method, p)
    if node in ('ms1', 'ms2', 'ms3'):
        return _scan_line(method, node)
    if node == 'ms2q':
        return _scan_line(method, 'ms2_quant')
    if node == 'ms3q':
        return _scan_line(method, 'ms3_quant')
    if node == 'faims':
        return ', '.join(str(c) for c in v('faims.cv_values')) + ' V'
    if node == 'deconv':
        return (f"z {v('deconvolution.min_charge')}–{v('deconvolution.max_charge')} · "
                f"{_num(v('deconvolution.min_mass'))}–{_num(v('deconvolution.max_mass'))} Da · "
                f"{'/'.join(str(t) for t in v('deconvolution.tol') or [])} ppm")
    if node == 'select':
        text = (f"{v('precursor_selection.max_precursors')} per survey · by {v('precursor_selection.rank_by')} · "
                f"{_num(v('precursor_selection.rt_window'))} s exclusion")
        t = v('precursor_selection.targeting')
        return text + (f' · {t}' if t != 'none' else '')
    if node == 'quant':
        return (f"{v('quantification.labelling')} · reporters from {v('quantification.scan')} · "
                f"{v('quantification.test')}")
    if node == 'ident':
        ref = v('identification.reference.fasta')
        ref = os.path.basename(ref) if ref else ('sequence' if v('identification.reference.sequence') else 'no reference')
        thr = v('identification.threshold')
        return f"{v('identification.method')} · {ref}" + (f' · closes ≥ {_num(thr)}' if thr is not None else '')
    if node == 'ms3plan':
        return f"{v('identification.characterization.mode')} · {v('identification.characterization.max_targets')} targets"
    if node == 'sweep2':
        return _sweep_line(method, 'precursor_selection.exploration')
    if node == 'sweep3':
        return _sweep_line(method, 'identification.characterization.exploration')
    if node == 'ms2x':
        return ', '.join(v('precursor_selection.additional_scans') or []) or 'none named'
    return ''



def _fmt(v) -> str:
    if isinstance(v, bool):
        return 'on' if v else 'off'
    if isinstance(v, (list, dict)):
        return json.dumps(v) if v else 'none'
    if isinstance(v, float) and v.is_integer():
        v = int(v)
    if isinstance(v, int):
        return f'{v:,}'.replace(',', ' ')
    return '(empty)' if v == '' else str(v)


def default_value(path):
    """The value of the key in FLASHIda's defaults (the schema's default, else
    the scan settings of FLASHIda defaults.json), or None."""
    key = _keys()[path]
    if key.default is not None:
        return key.default
    v = fm.get(fm.default_method(), path)
    return None if v in (None, '') else v


def default_of(method, path) -> str:
    """The default of a key, as the editor shows it: FLASHIda's default value,
    for a scan key MS1's value where the key takes MS1's, else what leaving it
    out means."""
    key = _keys()[path]
    unit = f' {key.unit}' if key.unit else ''
    text = key.default_text or ''
    if text.startswith('inherits ms1'):
        ms1 = 'ms_settings.ms1.' + path.split('.')[-1]
        v = fm.get(method, ms1)
        v = v if v not in (None, '') else default_value(ms1)
        if v not in (None, ''):
            return f'{_fmt(v)}{unit} (from MS1)'
    v = default_value(path)
    if v is not None:
        return _fmt(v) + (unit if not isinstance(v, (bool, list, dict, str)) else '')
    if text.startswith('absent'):
        rest = text.split(':', 1)[1].strip().replace('`', '') if ':' in text else ''
        return rest or 'not set'
    if text.startswith('required'):
        return 'none, ' + text
    if text == 'instrument method default':
        return 'not set'
    return text or 'not set'


TAB_PARTS = {
    'Global': ('📋', [('global', 'main')]),
    'MS1': ('🔬', [('ms1', 'main'), ('faims', 'option')]),
    'Deconvolution': ('🧮', [('deconv', 'main')]),
    'Precursor selection': ('🎯', [('select', 'main')]),
    'MS2': ('🔬', [('ms2', 'main'), ('sweep2', 'option'), ('ms2x', 'option')]),
    'Identification': ('🧬', [('ident', 'switch')]),
    'MS3': ('🧩', [('ms3plan', 'switch'), ('ms3', 'with'), ('sweep3', 'option')]),
    'Quantification': ('⚖️', [('quant', 'switch'), ('ms2q', 'with'), ('ms3q', 'with')]),
}
PARTS['global'] = ('Global', '📋', 'global', GLOBAL_GROUPS, False,
                   'Name, duration and log folder of the run, and how the host schedules commands.')
TITLES = {'ms1': 'MS1 survey scan', 'ms2': 'MS2 scan', 'ms3': 'MS3 scan', 'faims': 'FAIMS',
          'sweep2': 'Energy sweep', 'sweep3': 'Energy sweep', 'ms2x': 'Additional MS2 scans',
          'ms2q': 'Quant MS2 scan', 'ms3q': 'Quant MS3 scan'}
SWITCH_TEXT = {'ident': 'Identify in real time', 'ms3plan': 'Acquire MS3 on identified proteoforms',
               'quant': 'Quantify with isobaric labels'}


def run_line(method):
    g = lambda p: value(method, p)
    return f"{g('global.method_name') or 'Untitled'} · {_num(g('global.duration'))} min"
