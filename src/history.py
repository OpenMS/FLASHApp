"""The history of a result: how it came about, its settings and its processing report.

Every start stores, with each run it makes and with its comparison, a record
of the start (``HISTORY_KEY``: when, which files, the analysis, the run whose
settings or results it followed) and keeps the start's processing report in
the result's folder (``keep``, the template's ``processing_report.keep_report``),
so a later start no longer writes over it.

``timeline`` follows those records back: a run with identification added to
it starts with the run whose deconvolution it reused, a comparison with the
starts of all the runs it compares. ``changed_settings`` lists the tool
parameters a result was made with that differ from the tools' defaults.
Runs made before the history was recorded show what their stored settings
tell (``analysis_settings``, ``FD_parameters``, ``FTnT_parameters``).

Streamlit-free: ``keep`` and ``record`` run in the workflow process.
"""

import re
import time
import xml.etree.ElementTree as ET
from datetime import datetime
from os.path import basename
from pathlib import Path

from src import analysis
from src.quant import isobaric
from src.workflow.processing_report import keep_report, report_for, to_markdown
from src.workflow.parameter_names import display_name

HISTORY_KEY = 'history'
# The kept processing report, in the result's folder of the cache
REPORT_DIR = 'processing-report'
REPORT_TAG = 'processing_report_md'
REPORT_NAME = 'processing-report.md'
SETTINGS_KEY = 'analysis_settings'

# The tool parameters a stored run keeps (``TopDownWorkflow._store_parameters``)
STORED_PARAMETERS = (('FLASHDeconv', 'FD_parameters'), ('FLASHTnT', 'FTnT_parameters'))
# What the history leaves out of the settings: set by the workflow, not by the user
NOT_SETTINGS = {'threads', 'in', 'out', 'in_id', 'out_id', 'in_feature', 'database', 'fasta', 'ida_log'}

STARTED, PROCESSED, IDENTIFICATION, FILES_ADDED, COMPARED = (
    'Started', 'Processed', 'Identification added', 'Files added', 'Compared',
)

_STAMP = re.compile(r'_(\d{8}-\d{6})$')


# -- recorded by the workflow ----------------------------------------------

def start_record(start, files, identify, quant, plex, follows=None, source=None):
    """What ``HISTORY_KEY`` stores for a result of one start.

    Args:
        start: When the start began (seconds), the same for all its results.
        files: The input files of this result (names).
        identify, quant, plex: The analysis.
        follows: The start's follow-up (``follow_up_info()``): the result whose
            settings it started from, the runs it added files to.
        source: The run whose deconvolution this run reused (identification added).
    """
    follows = follows or {}
    return {
        'start': float(start),
        'files': [basename(str(f).rstrip('/\\')) for f in files],
        'analysis': analysis.name(identify, quant, plex if quant == analysis.ISOBARIC else None),
        'follows': {
            'name': follows.get('name') or '',
            'note': follows.get('note') or '',
            'adds_to': list(follows.get('adds_to') or []),
        },
        'source': source,
    }


def record(file_manager, dataset_id, entry):
    file_manager.store_data(dataset_id, HISTORY_KEY, entry)


def report_dir(file_manager, dataset_id):
    """Where the processing report of the start that made ``dataset_id`` is kept."""
    return Path(file_manager.cache_path, 'files', dataset_id, REPORT_DIR)


def keep(workflow_dir, file_manager, dataset_ids, settings=None):
    """Keep the processing report of the start that just ended with each of
    its results, and its Markdown among their files (exported with them)."""
    for dataset_id in dataset_ids:
        folder = report_dir(file_manager, dataset_id)
        if not keep_report(workflow_dir, folder):
            continue
        try:
            markdown = to_markdown(report_for(folder, settings=settings))
        except Exception:  # a report that cannot be written leaves the result as it is
            continue
        path = Path(folder, REPORT_NAME)
        path.write_text(markdown, encoding='utf-8')
        file_manager.store_file(dataset_id, REPORT_TAG, path, file_name=REPORT_NAME)


# -- read on the pages --------------------------------------------------------

def _data(file_manager, dataset_id, name):
    if not file_manager.result_exists(dataset_id, name):
        return None
    try:
        return file_manager.get_results(dataset_id, [name])[name]
    except Exception:
        return None


def history_of(file_manager, dataset_id):
    entry = _data(file_manager, dataset_id, HISTORY_KEY)
    return entry if isinstance(entry, dict) else None


def has_report(file_manager, dataset_id):
    return Path(report_dir(file_manager, dataset_id), 'logs', 'progress.json').exists()


def stamp_time(dataset_id):
    """The time in a dataset id (``<name>_YYYYmmdd-HHMMSS``), in seconds, else None."""
    found = _STAMP.search(str(dataset_id))
    if not found:
        return None
    try:
        return time.mktime(datetime.strptime(found.group(1), '%Y%m%d-%H%M%S').timetuple())
    except ValueError:
        return None


def format_time(seconds):
    if not seconds:
        return ''
    when = datetime.fromtimestamp(seconds)
    return f"{when.day} {when.strftime('%b %Y, %H:%M')}"


def is_comparison(file_manager, dataset_id):
    return file_manager.result_exists(dataset_id, 'quant_manifest')


def compared_runs(file_manager, dataset_id):
    """The runs a comparison was made from (dataset ids)."""
    manifest = _data(file_manager, dataset_id, 'quant_manifest') or {}
    runs = (manifest.get('sources') or {}).get('runs') or manifest.get('fractions') or []
    return [r for r in runs if isinstance(r, str)]


def comparisons_of(file_manager, dataset_id):
    """The comparisons that compare run ``dataset_id``, newest first."""
    found = sorted(file_manager.get_results_list(['quant_manifest']), reverse=True)
    return [c for c in found if dataset_id in compared_runs(file_manager, c)]


def _event(what, when, files, detail='', group=None, dataset=None):
    return {'what': what, 'when': when, 'files': list(files), 'detail': detail,
            'group': group, 'datasets': [dataset] if dataset else []}


def _run_files(file_manager, dataset_id, entry):
    if entry and entry.get('files'):
        return entry['files']
    record_ = _data(file_manager, dataset_id, SETTINGS_KEY)
    if isinstance(record_, dict) and record_.get('mzml'):
        return [basename(str(record_['mzml']))]
    return []


def _run_analysis(file_manager, dataset_id):
    record_ = _data(file_manager, dataset_id, SETTINGS_KEY)
    if not isinstance(record_, dict):
        return ''
    return analysis.name(record_.get('identify'), record_.get('quantification'), record_.get('plex'))


def _run_events(file_manager, dataset_id, seen):
    if dataset_id in seen:
        return []
    seen.add(dataset_id)
    entry = history_of(file_manager, dataset_id)
    files = _run_files(file_manager, dataset_id, entry)
    if entry is None:
        # Made before the history was recorded: what its stored settings tell
        return [_event(PROCESSED, stamp_time(dataset_id), files, _run_analysis(file_manager, dataset_id),
                       dataset=dataset_id)]
    group = entry.get('start')
    follows = entry.get('follows') or {}
    source = entry.get('source')
    if source:
        events = _run_events(file_manager, source, seen) if file_manager.result_exists(source, SETTINGS_KEY) else []
        detail = f"deconvolution of {file_manager.get_display_name(source)} reused; {entry.get('analysis', '')}"
        return events + [_event(IDENTIFICATION, entry.get('start'), files, detail, group, dataset_id)]
    if follows.get('adds_to'):
        earlier = ', '.join(file_manager.get_display_name(d) for d in follows['adds_to'])
        detail = f"to {earlier}, with its settings; {entry.get('analysis', '')}"
        return [_event(FILES_ADDED, entry.get('start'), files, detail, group, dataset_id)]
    detail = entry.get('analysis', '')
    if follows.get('name'):
        detail += f"; settings of {follows['name']}"
    return [_event(STARTED, entry.get('start'), files, detail, group, dataset_id)]


def _comparison_detail(file_manager, dataset_id):
    manifest = _data(file_manager, dataset_id, 'quant_manifest') or {}
    if manifest.get('kind') == 'isobaric':
        what = f"{isobaric.label(manifest.get('plex', ''))}, {len(manifest.get('runs') or [])} samples"
    else:
        what = f"label-free, {len(compared_runs(file_manager, dataset_id))} runs"
    design = _data(file_manager, dataset_id, 'quant_design')
    try:
        conditions = [c for c in dict.fromkeys(design['condition'].to_list()) if c]
    except Exception:
        conditions = []
    if conditions:
        what += ', conditions ' + ' / '.join(f'*{c}*' for c in conditions)
    return what


def timeline(file_manager, dataset_id):
    """How ``dataset_id`` came about, oldest first: one event per start
    (``what``, ``when`` in seconds, ``files``, ``detail``, ``datasets``)."""
    seen = set()
    if is_comparison(file_manager, dataset_id):
        events = []
        for run in compared_runs(file_manager, dataset_id):
            if file_manager.result_exists(run, SETTINGS_KEY) or history_of(file_manager, run):
                events += _run_events(file_manager, run, seen)
        entry = history_of(file_manager, dataset_id) or {}
        events.append(_event(COMPARED, entry.get('start') or stamp_time(dataset_id), [],
                             _comparison_detail(file_manager, dataset_id), dataset=dataset_id))
    else:
        events = _run_events(file_manager, dataset_id, seen)
    # The runs of one start are one event
    merged, by_group = [], {}
    for event in events:
        key = (event['what'], event['group']) if event['group'] is not None else None
        if key is not None and key in by_group:
            kept = by_group[key]
            kept['files'] += [f for f in event['files'] if f not in kept['files']]
            kept['datasets'] += event['datasets']
            continue
        if key is not None:
            by_group[key] = event
        merged.append(event)
    return sorted(merged, key=lambda e: (e['what'] == COMPARED, e['when'] or 0))


def event_markdown(event):
    """One event as a line of Markdown."""
    text = f"**{event['what']}**"
    if event['when']:
        text += f" {format_time(event['when'])}"
    if event['files']:
        names = event['files']
        shown = ', '.join(f'`{n}`' for n in names[:4])
        if len(names) > 4:
            shown += f' and {len(names) - 4} more'
        text += f", on {shown}"
    if event['detail']:
        text += f" ({event['detail']})"
    return text


# -- settings ---------------------------------------------------------------------

def ini_defaults(ini_path):
    """parameter path -> default value (text) of a TOPP tool's .ini, ``{}`` without one."""
    try:
        root = ET.parse(ini_path).getroot()
    except (ET.ParseError, OSError):
        return {}
    found = {}

    def tag(el):
        return el.tag.rsplit('}', 1)[-1] if isinstance(el.tag, str) else ''

    def walk(el, parts):
        for child in el:
            kind = tag(child)
            if kind == 'NODE':
                walk(child, parts + [child.get('name') or ''])
            elif kind == 'ITEM':
                found[':'.join(parts[2:] + [child.get('name') or ''])] = child.get('value') or ''
            elif kind == 'ITEMLIST':
                values = [i.get('value') or '' for i in child if tag(i) == 'LISTITEM']
                found[':'.join(parts[2:] + [child.get('name') or ''])] = ', '.join(values)

    # <PARAMETERS><NODE name="Tool"><NODE name="1"> ...: the path starts below those two
    walk(root, [])
    return found


def _text(value):
    if isinstance(value, bool):
        return 'true' if value else 'false'
    if isinstance(value, (list, tuple)):
        return ', '.join(str(v) for v in value)
    return str(value).replace('\n', ', ')


def _same(value, default):
    if value == default:
        return True
    try:
        return float(value) == float(default)
    except (TypeError, ValueError):
        return False


def _label(tool, parameter):
    entry = display_name(tool, parameter)
    return entry['name'] if entry else parameter


def changed_settings(file_manager, dataset_id, workflow_dir=None):
    """[(tool, parameter, default, value)] the result was made with that differ
    from the tools' defaults: from its kept processing report, else from the
    tool parameters stored with the run (``workflow_dir``'s .ini for defaults)."""
    folder = report_dir(file_manager, dataset_id)
    per_tool = {}
    if has_report(file_manager, dataset_id):
        try:
            steps = report_for(folder, settings={})['steps']
        except Exception:
            steps = []
        for step in steps:
            if step.get('kind') != 'topp':
                continue
            tool = step.get('program') or step.get('tool')
            values = per_tool.setdefault(tool, {})
            for key, value in (step.get('params') or {}).items():
                values.setdefault(key, [])
                if value not in values[key]:
                    values[key].append(value)
    for tool, name in STORED_PARAMETERS:
        # A run whose deconvolution was reused has no FLASHDeconv call in its report
        if tool in per_tool:
            continue
        stored = _data(file_manager, dataset_id, name)
        if isinstance(stored, dict) and stored:
            per_tool[tool] = {k: [_text(v)] for k, v in stored.items()}
    rows = []
    for tool, values in per_tool.items():
        defaults = ini_defaults(Path(folder, 'ini', f'{tool}.ini'))
        if not defaults and workflow_dir is not None:
            defaults = ini_defaults(Path(workflow_dir, 'ini', f'{tool}.ini'))
        for key, found in values.items():
            if key in NOT_SETTINGS or isinstance(found, dict):
                continue
            value = ' | '.join(found)
            default = defaults.get(key)
            if default is not None and len(found) == 1 and _same(found[0], default):
                continue
            rows.append((tool, _label(tool, key), '' if default is None else default, value))
    return rows
