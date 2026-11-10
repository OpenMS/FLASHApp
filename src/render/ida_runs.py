"""The runs table of the FLASHIda acquisition page and its Method dialog.

Each run is a raw file name and the method FLASHIda gets for it. The method
starts from a preset (FLASHIda's defaults, the presets shipped in
assets/flashida/presets.json, or the user's own, kept once per instrument PC)
and is edited by section, key by key, or as JSON. While the acquisition runs,
the table shows which runs the newest sequence (.sld) schedules and what
FLASHIda is doing for each.
"""

import copy
import json
import time
from pathlib import Path

import pandas as pd
import streamlit as st

from src import flashida
from src import flashida_method as fm
from src import flashida_sequence as fs

# Human names for the keys the sections show (the rest keep FLASHIda's name).
LABELS = {
    'global.method_name': 'Method name',
    'global.method_description': 'Description',
    'global.duration': 'Duration',
    'deconvolution.tqscore_threshold': 'Qscore that closes a species',
    'deconvolution.score_threshold': 'Lowest qscore',
    'deconvolution.min_charge': 'Lowest charge',
    'deconvolution.max_charge': 'Highest charge',
    'deconvolution.min_mass': 'Lowest mass',
    'deconvolution.max_mass': 'Highest mass',
    'deconvolution.tol': 'Tolerance MS1, MS2, MS3',
    'precursor_selection.max_precursors': 'Precursors per survey',
    'precursor_selection.rank_by': 'Rank species by',
    'precursor_selection.charge_selection': 'Isolate the charge with the best',
    'precursor_selection.precursor_charges': 'Charge states per species',
    'precursor_selection.consider_all_charges': 'Rank by the best charge',
    'precursor_selection.rt_window': 'Exclusion window',
    'precursor_selection.tie_threshold': 'Tie threshold',
    'precursor_selection.targeting': 'Target list',
    'precursor_selection.strict_inclusion': 'Inclusion list only',
    'identification.reference.fasta': 'Protein database (FASTA)',
    'identification.threshold': 'Identified at score',
    'identification.method': 'Search',
    'identification.flashtnt.min_length': 'Shortest tag',
    'identification.flashtnt.max_length': 'Longest tag',
    'identification.flashtnt.max_blind_mod_count': 'Unknown modifications',
    'identification.flashtnt.max_mod_mass': 'Largest unknown modification',
    'identification.flashtnt.var_mod': 'Variable modifications',
    'quantification.enabled': 'Quantify',
    'quantification.labelling': 'Labels',
    'quantification.scan': 'Reporter ions from',
    'quantification.reporter_mz_tol': 'Reporter m/z tolerance',
    'quantification.fold_change_threshold': 'Fold change',
    'faims.cv_values': 'FAIMS CVs',
    'scheduling.agc_interval_seconds': 'AGC prescan every',
    'files.inclusion_list': 'Inclusion list',
    'files.target_logs': 'Target logs',
    'runtime.log_dir': 'Log folder',
}

_SCAN_KEYS = ('analyzer', 'activation', 'collision_energy', 'resolution', 'agc_target',
              'max_it', 'first_mass', 'last_mass', 'microscans', 'data_type')

SECTIONS = {
    'Run': ['global.method_name', 'global.duration', 'global.method_description',
            'runtime.log_dir', 'files.inclusion_list', 'faims.cv_values'],
    'Deconvolution': ['deconvolution.min_charge', 'deconvolution.max_charge',
                      'deconvolution.min_mass', 'deconvolution.max_mass',
                      'deconvolution.tol', 'deconvolution.tqscore_threshold'],
    'Precursor selection': ['precursor_selection.max_precursors', 'precursor_selection.rank_by',
                            'precursor_selection.charge_selection',
                            'precursor_selection.precursor_charges',
                            'precursor_selection.rt_window', 'precursor_selection.targeting'],
    'Identification': ['identification.reference.fasta', 'identification.threshold',
                       'identification.flashtnt.min_length', 'identification.flashtnt.max_length',
                       'identification.flashtnt.max_blind_mod_count',
                       'identification.flashtnt.max_mod_mass'],
    'Quantification': ['quantification.enabled', 'quantification.labelling',
                       'quantification.scan', 'quantification.reporter_mz_tol',
                       'quantification.fold_change_threshold'],
    'MS1 scan': [f'ms_settings.ms1.{k}' for k in _SCAN_KEYS if k not in ('activation', 'collision_energy')],
    'MS2 scan': [f'ms_settings.ms2.{k}' for k in _SCAN_KEYS],
}

_VERSION = 'ida-method-version'
OPEN_KEY = 'ida-open-method'
DIALOG_KEY = 'ida-open-dialog'  # 'import' or 'pattern'
CLIPBOARD_KEY = 'ida-clipboard'  # {name, preset, method} of the copied run


def _keys() -> dict:
    if '_ida_keys' not in st.session_state:
        st.session_state['_ida_keys'] = fm.keys()
    return st.session_state['_ida_keys']


def label(path: str) -> str:
    key = _keys().get(path)
    name = LABELS.get(path)
    if name is None:
        name = path.split('.')[-1].replace('_', ' ')
        name = name[:1].upper() + name[1:]
        if path.startswith('ms_settings.'):
            scan = path.split('.')[1].upper().replace('_QUANT', ' quant')
            name = f'{name}'
            if not any(path in v for v in SECTIONS.values()):
                name = f'{scan} {name.lower()}'
    if key is not None and key.unit:
        name += f' ({key.unit})'
    return name


def _help(key: fm.Key) -> str:
    text = f'`{key.path}`\n\n{key.description}'
    if key.enum_help:
        text += '\n\n' + '\n\n'.join(f'**{v or "(none)"}**: {h}' for v, h in zip(key.enum, key.enum_help))
    if key.default is None and key.default_text:
        text += f'\n\nLeft empty: {key.default_text}.'
    return text


def _version() -> int:
    return st.session_state.get(_VERSION, 0)


def _bump() -> None:
    st.session_state[_VERSION] = _version() + 1


def key_widget(method: dict, path: str, wkey: str, disabled: bool = False) -> bool:
    """One widget for one key of a method; True when it changed the method.
    A cleared field removes the key, so FLASHIda uses its own default."""
    key = _keys().get(path)
    if key is None:
        return False
    present = fm.get(method, path, missing=KeyError) is not KeyError
    value = fm.get(method, path) if present else key.default
    help_ = _help(key)
    name = label(path)

    if key.enum:
        options = list(key.enum)
        if not present and key.default is None:
            options = [None] + options
        new = st.selectbox(
            name, options, index=options.index(value) if value in options else 0, key=wkey,
            help=help_, disabled=disabled,
            format_func=lambda v: ('(instrument method)' if v in (None, '') else str(v)),
        )
    elif key.type == 'boolean':
        new = st.toggle(name, value=bool(value), key=wkey, help=help_, disabled=disabled)
    elif key.type in ('integer', 'number'):
        is_int = key.type == 'integer'
        shown = value if isinstance(value, (int, float)) and not isinstance(value, bool) else None
        if shown is not None:
            shown = int(shown) if is_int else float(shown)
        new = st.number_input(
            name, value=shown, step=1 if is_int else None, key=wkey, help=help_,
            disabled=disabled, placeholder=key.default_text or 'FLASHIda default',
            format=None if is_int else '%g',
        )
        if new is not None and not is_int and float(new).is_integer() and isinstance(value, int):
            new = int(new)
    elif key.type == 'array':
        if key.items_type == 'object':
            st.text_input(name, value=json.dumps(value) if value is not None else '',
                          key=wkey, help=help_, disabled=True)
            st.caption('Edit this one in the JSON tab.')
            return False
        text = ', '.join(str(v) for v in (value or []))
        raw = st.text_input(name, value=text, key=wkey, help=help_ + '\n\nComma-separated.',
                            disabled=disabled)
        if raw == text:
            return False
        parts = [p.strip() for p in raw.split(',') if p.strip()]
        try:
            if key.items_type == 'integer':
                new = [int(p) for p in parts]
            elif key.items_type == 'number':
                new = [int(p) if p.lstrip('-').isdigit() else float(p) for p in parts]
            else:
                new = parts
        except ValueError:
            st.error(f'{name}: numbers only, separated by commas.')
            return False
    else:
        text = '' if value is None else str(value)
        new = st.text_input(name, value=text, key=wkey, help=help_, disabled=disabled,
                            placeholder=key.default_text or '')
        if new == text:
            return False
        if new == '' and key.default in (None, ''):
            new = None

    if new == value and present:
        return False
    if new is None:
        if not present:
            return False
        fm.drop(method, path)
    else:
        if new == value and not present:
            return False
        fm.put(method, path, new)
    return True


# --- the dialog -------------------------------------------------------------------

def _store(workflow_dir, runs) -> None:
    fm.save_runs(workflow_dir, runs)


def _apply_preset(workflow_dir, presets_dir, index, select_key) -> None:
    runs = fm.load_runs(workflow_dir)
    name = st.session_state[select_key]
    presets = fm.all_presets(presets_dir)
    if name in presets:
        runs[index].update(preset=name, method=copy.deepcopy(presets[name]))
        _store(workflow_dir, runs)
        _bump()


# Closing it with ✕ reruns the page too, so Run analysis sees the method as it is now
@st.dialog('Method', width='large', on_dismiss='rerun')
def method_dialog(workflow_dir, presets_dir, index: int, read_only: bool = False) -> None:
    runs = fm.load_runs(workflow_dir)
    if index >= len(runs):
        st.rerun()
    run = runs[index]
    method = run['method']
    v = _version()
    st.markdown(f"**{run['name'] or 'Unnamed run'}**")

    presets = fm.all_presets(presets_dir)
    users = fm.user_presets(presets_dir)
    names = list(presets)
    c1, c2, c3 = st.columns([3, 1.3, 1.3], vertical_alignment='bottom')
    select_key = f'ida-preset-{index}-{v}'
    current = run.get('preset') if run.get('preset') in presets else None
    c1.selectbox(
        'Start from preset', names, index=names.index(current) if current else None,
        key=select_key, placeholder='Pick a preset', disabled=read_only,
        on_change=_apply_preset, args=(workflow_dir, presets_dir, index, select_key),
        help='Picking a preset replaces the settings below. Your own presets are '
             'kept on this computer, for every workspace.',
    )
    with c2.popover('💾 Save as preset', use_container_width=True, disabled=read_only):
        new_name = st.text_input('Preset name', value=run.get('preset') if run.get('preset') in users else '',
                                 key=f'ida-preset-name-{index}-{v}')
        builtin = new_name.strip() in fm.builtin_presets()
        if builtin:
            st.caption('Built-in presets keep their name; pick another one.')
        if st.button('Save', type='primary', disabled=not new_name.strip() or builtin,
                     key=f'ida-preset-save-{index}-{v}'):
            fm.save_user_preset(presets_dir, new_name.strip(), method)
            run['preset'] = new_name.strip()
            _store(workflow_dir, runs)
            _bump()
            st.rerun(scope='fragment')
    deletable = run.get('preset') in users and run.get('preset') not in fm.builtin_presets()
    if c3.button('🗑️ Delete preset', use_container_width=True, disabled=read_only or not deletable,
                 key=f'ida-preset-del-{index}-{v}',
                 help=f"Delete your preset '{run.get('preset')}' from this computer."
                 if deletable else 'Only your own presets can be deleted.'):
        fm.delete_user_preset(presets_dir, run['preset'])
        _bump()
        st.rerun(scope='fragment')

    if current:
        diff = fm.changes(presets[current], method)
        if diff:
            shown = ', '.join(f'{label(p).split(" (")[0].lower()} {_short(a)} → {_short(b)}'
                              for p, a, b in diff[:4])
            more = f' and {len(diff) - 4} more' if len(diff) > 4 else ''
            st.caption(f'✏️ Changed from **{current}**: {shown}{more}')

    keys = _keys()
    tabs = st.tabs(list(SECTIONS) + [f'All {len(keys)} keys', 'JSON'])
    changed = False
    for tab, (section, paths) in zip(tabs, SECTIONS.items()):
        with tab:
            cols = st.columns(2)
            for n, path in enumerate(paths):
                with cols[n % 2]:
                    changed |= key_widget(method, path, f'ida-w-{index}-{v}-{path}', read_only)

    with tabs[-2]:
        query = st.text_input('Find a key', key=f'ida-q-{index}',
                              placeholder='Part of a key or its description, e.g. faims, agc, exploration')
        q = query.strip().lower()
        matching = [p for p, k in keys.items()
                    if not q or q in p.lower() or q in k.description.lower() or q in label(p).lower()]
        if matching:
            pick = st.selectbox('Key', matching, key=f'ida-pick-{index}-{q}',
                                format_func=lambda p: f'{p}  ·  {label(p)}')
            changed |= key_widget(method, pick, f'ida-w-{index}-{v}-all-{pick}', read_only)
            rows = [{'Key': p, 'Value': _short(fm.get(method, p), 60) if fm.get(method, p, KeyError) is not KeyError else '',
                     'Default': _short(keys[p].default, 40) if keys[p].default is not None else keys[p].default_text,
                     'Unit': keys[p].unit} for p in matching]
            st.dataframe(pd.DataFrame(rows), hide_index=True, height=260, use_container_width=True)
        else:
            st.info('No key matches.')

    with tabs[-1]:
        text_key = f'ida-json-{index}-{v}'
        text = st.text_area('Method as JSON', value=json.dumps(method, indent=2), height=380,
                            key=text_key, disabled=read_only,
                            help='The file FLASHIda gets for this run. Apply checks it against '
                                 "FLASHIda's method schema.")
        a, b, c = st.columns([1, 1, 1.2], vertical_alignment='bottom')
        if a.button('Apply JSON', type='primary', use_container_width=True, disabled=read_only,
                    key=f'ida-json-apply-{index}-{v}'):
            try:
                parsed = json.loads(text)
            except ValueError as e:
                st.error(f'Not valid JSON: {e}')
            else:
                if not isinstance(parsed, dict):
                    st.error('A method is a JSON object ({...}).')
                else:
                    run['method'] = parsed
                    _store(workflow_dir, runs)
                    _bump()
                    st.rerun(scope='fragment')
        b.download_button('⬇️ Download', json.dumps(method, indent=2), use_container_width=True,
                          file_name=fm.method_file_name(run['name'] or 'method'),
                          mime='application/json', key=f'ida-json-dl-{index}-{v}')
        with c.popover('📥 Import .json file', use_container_width=True, disabled=read_only):
            upload = st.file_uploader('Method file', type=['json'], key=f'ida-json-up-{index}-{v}')
            if upload is not None:
                try:
                    parsed = json.loads(upload.getvalue().decode('utf-8-sig'))
                    assert isinstance(parsed, dict)
                except (ValueError, AssertionError):
                    st.error('This file is not a JSON method.')
                else:
                    run['method'] = parsed
                    _store(workflow_dir, runs)
                    _bump()
                    st.rerun(scope='fragment')

    if changed:
        _store(workflow_dir, runs)

    problems = fm.check(method)
    if problems:
        with st.container(border=True):
            st.markdown(f'⚠️ **This FLASHApp expects FLASHIda to refuse this method** ({len(problems)} '
                        f'problem{"s" if len(problems) != 1 else ""}):')
            st.markdown('\n'.join(f'- `{p}`' for p in problems[:8]))
            st.caption('The check uses the method schema bundled with FLASHApp, which may be older than '
                       'your FLASHIda. **Run anyway** on the page starts the run with this method as it is.')
    unknown = fm.unknown_keys(method)
    if unknown:
        with st.container(border=True):
            st.markdown(f'ℹ️ **{len(unknown)} key{"s" if len(unknown) != 1 else ""} this FLASHApp does not know** '
                        '(a newer FLASHIda may). They are kept and passed to FLASHIda as they are, which checks '
                        'them when the run starts; edit or remove them in the JSON tab.')
            st.markdown(', '.join(f'`{k}`' for k in unknown[:20]) + (' …' if len(unknown) > 20 else ''))

    st.divider()
    a, b = st.columns([1.5, 1])
    if a.button('Use this method for all runs', use_container_width=True, disabled=read_only,
                key=f'ida-all-{index}-{v}'):
        for other in runs:
            other.update(preset=run.get('preset'), method=copy.deepcopy(method))
        _store(workflow_dir, runs)
        st.toast(f'All {len(runs)} runs use this method.')
    if b.button('Done', type='primary', use_container_width=True, key=f'ida-done-{index}-{v}'):
        st.rerun()


def _short(value, limit: int = 24) -> str:
    text = json.dumps(value) if not isinstance(value, str) else (value or '""')
    return text if len(text) <= limit else text[:limit - 1] + '…'


# --- import, export and runs by pattern ---------------------------------------------

_MODES = {'append': 'Append to the table', 'replace': 'Replace the table'}
_COPIED = '📋 Copied method'


def _method_choices(presets: dict) -> dict:
    """{label: (preset, method)} a dialog offers for new runs: the copied
    method first (when there is one), then every preset."""
    choices = {}
    clip = st.session_state.get(CLIPBOARD_KEY)
    if clip:
        choices[f"{_COPIED} ({clip['name'] or 'unnamed run'})"] = (clip.get('preset'), clip['method'])
    for name, method in presets.items():
        choices[name] = (name, method)
    return choices


def _mode_choice(runs: list, key: str) -> str:
    return st.radio('Into the table', list(_MODES), key=key, horizontal=True,
                    format_func=lambda m: _MODES[m] if m == 'append' or not runs
                    else f'{_MODES[m]} ({len(runs)} run{"s" if len(runs) != 1 else ""} now)')


def _duplicates_note(duplicates: list) -> None:
    if duplicates:
        st.warning(f'{len(duplicates)} name{"s are" if len(duplicates) != 1 else " is"} in the table '
                   'more than once: ' + ', '.join(f'`{d}`' for d in duplicates[:10])
                   + (' …' if len(duplicates) > 10 else '')
                   + '. The runs are added and marked; rename or remove them before the run can start.')


def _commit(workflow_dir, runs: list, added: int, duplicates: list) -> None:
    fm.save_runs(workflow_dir, runs)
    _bump()
    st.session_state['ida-added'] = (added, duplicates)
    st.rerun()


@st.dialog('Import runs', width='large', on_dismiss='rerun')
def import_dialog(workflow_dir, presets_dir) -> None:
    runs = fm.load_runs(workflow_dir)
    source = st.radio('Source', ['json', 'sld'], horizontal=True, key='ida-imp-source',
                      format_func=lambda s: {'json': 'Sequence file (JSON)',
                                             'sld': 'Xcalibur sequence (.sld)'}[s])
    upload = st.file_uploader('File', type=['json'] if source == 'json' else ['sld'],
                              key=f'ida-imp-file-{source}',
                              help='A file exported here with **Export sequence**.' if source == 'json' else
                              'Any Xcalibur sequence. Only the raw file names are taken.')
    mode = _mode_choice(runs, 'ida-imp-mode')

    new = []
    choices = _method_choices(fm.all_presets(presets_dir))
    pick = st.selectbox('Method for the new runs', list(choices), key='ida-imp-method') \
        if source == 'sld' else None
    if upload is not None:
        try:
            if source == 'json':
                new = fs.parse_runs(upload.getvalue())
            else:
                preset, method = choices[pick]
                new = fs.runs_from_names(fs.sld_names(upload.getvalue()), preset, method)
        except ValueError as e:
            st.error(str(e))
            new = []
        else:
            if not new:
                st.info('The file holds no runs.')
            else:
                st.caption(f'{len(new)} run{"s" if len(new) != 1 else ""} in the file'
                           + (', each with its own method.' if source == 'json' else
                              '. Only the names are taken; every run gets the method picked above.'))
                st.dataframe(pd.DataFrame({'Raw file name': [r['name'] for r in new],
                                           'Preset': [r['preset'] or 'Own method' for r in new]}),
                             hide_index=True, use_container_width=True, height=min(36 + 35 * len(new), 250))
    combined, duplicates = fs.combine(runs, new, mode)
    if new:
        _duplicates_note(duplicates)
    if st.button(f'Import {len(new)} run{"s" if len(new) != 1 else ""}' if new else 'Import',
                 type='primary', disabled=not new, key='ida-imp-go'):
        _commit(workflow_dir, combined, len(new), duplicates)


@st.dialog('Add runs by pattern', width='large', on_dismiss='rerun')
def pattern_dialog(workflow_dir, presets_dir) -> None:
    runs = fm.load_runs(workflow_dir)
    st.caption('Name = prefix_date_sample_comment_replicate_suffix; empty parts are left out.')
    a, b, c = st.columns(3)
    prefix = a.text_input('Prefix', key='ida-pat-prefix')
    date = b.text_input('Date', value=fs.today(), key='ida-pat-date')
    suffix = c.text_input('Suffix', key='ida-pat-suffix')
    comment = a.text_input('Comment', key='ida-pat-comment', help='Used for every sample without its own.')
    replicates = b.number_input('Replicates', min_value=1, max_value=99, value=1, step=1,
                                key='ida-pat-reps', help='Used for every sample without its own.')
    rfmt = c.selectbox('Replicate as', list(fs.REPLICATE_FORMATS), key='ida-pat-rfmt',
                       format_func=fs.REPLICATE_FORMATS.get)
    left, right = st.columns([1.3, 1])
    samples_text = left.text_area(
        'Samples, one per line', key='ida-pat-samples', height=130,
        placeholder='BSA\nEcoli; comment=top3; replicates=3\nHeLa',
        help='A line may set its own values after the name: `; comment=…` and `; replicates=…`.')
    order = right.radio('Order', list(fs.ORDERS), key='ida-pat-order', format_func=fs.ORDERS.get)
    choices = _method_choices(fm.all_presets(presets_dir))
    labels = list(choices)
    pick = right.selectbox('Method for the samples', labels, key='ida-pat-method')

    with st.container(border=True):
        w1, w2, w3, w4 = st.columns([0.9, 1, 2, 1.8], vertical_alignment='bottom')
        wash = w1.checkbox('Wash', key='ida-pat-wash')
        every = w2.number_input('Every … runs', min_value=1, max_value=99, value=3, step=1,
                                key='ida-pat-every', disabled=not wash)
        wash_name = w3.text_input('Wash name', value=fs.DEFAULT_WASH_NAME, key='ida-pat-wname',
                                  disabled=not wash, help='Takes {prefix}, {date}, {suffix} and {n}, the wash number.')
        wash_pick = w4.selectbox('Wash method', labels, disabled=not wash, key='ida-pat-wmethod',
                                 index=labels.index(fs.WASH_PRESET) if fs.WASH_PRESET in labels else 0)

    mode = _mode_choice(runs, 'ida-pat-mode')

    names, error = [], ''
    try:
        samples = fs.parse_samples(samples_text, comment.strip(), int(replicates))
    except ValueError as e:
        error = str(e)
    else:
        if wash and fs.unknown_placeholders(wash_name):
            error = 'Unknown placeholder in the wash name: ' + ', '.join(
                '{' + p + '}' for p in fs.unknown_placeholders(wash_name))
        else:
            names = fs.pattern_names(samples, prefix.strip(), date.strip(), suffix.strip(), rfmt, order,
                                     int(every) if wash else 0, wash_name)
    if error:
        st.error(error)
    new = []
    for n in names:
        preset, method = choices[wash_pick if n['wash'] else pick]
        new.append({'name': n['name'], 'preset': preset, 'method': copy.deepcopy(method)})
    combined, duplicates = fs.combine(runs, new, mode)
    dup = {d.casefold() for d in duplicates}
    if names:
        washes = sum(n['wash'] for n in names)
        start = len(runs) if mode == 'append' else 0
        st.markdown(f'**Preview:** {len(names) - washes} sample run{"s" if len(names) - washes != 1 else ""}'
                    + (f' and {washes} wash{"es" if washes != 1 else ""}' if washes else ''))
        st.dataframe(pd.DataFrame({
            '#': range(start + 1, start + len(names) + 1),
            'Raw file name': [n['name'] for n in names],
            'Method': ['Wash: ' + wash_pick if n['wash'] else pick for n in names],
            '': ['⚠️ duplicate' if n['name'].casefold() in dup else '' for n in names],
        }), hide_index=True, use_container_width=True, height=min(36 + 35 * len(names), 300))
        _duplicates_note(duplicates)
    elif not error:
        st.caption('Enter at least one sample to see the runs.')
    if st.button(f'Add {len(new)} run{"s" if len(new) != 1 else ""}' if new else 'Add runs',
                 type='primary', disabled=not new, key='ida-pat-go'):
        _commit(workflow_dir, combined, len(new), duplicates)


def _copy_method(run: dict) -> None:
    st.session_state[CLIPBOARD_KEY] = {'name': run['name'], 'preset': run.get('preset'),
                                       'method': copy.deepcopy(run['method'])}


def _clipboard_bar(workflow_dir, runs: list, selected: list) -> None:
    clip = st.session_state.get(CLIPBOARD_KEY)
    with st.container(border=True):
        a, b = st.columns([5, 1.9], vertical_alignment='center')
        if clip:
            a.markdown(f"📋 Copied: the method of **{clip['name'] or 'an unnamed run'}** "
                       f"({clip.get('preset') or 'own method'}). Tick runs, then paste.")
        else:
            a.caption('📋 No method copied. Copy one with the 📋 button of a run, tick the runs that '
                      'should get it, then paste.')
        if b.button(f'Paste into selection ({len(selected)})', use_container_width=True,
                    disabled=not clip or not selected, key='ida-paste'):
            for i in selected:
                runs[i].update(preset=clip.get('preset'), method=copy.deepcopy(clip['method']))
            fm.save_runs(workflow_dir, runs)
            _bump()
            st.toast(f'Pasted the method into {len(selected)} run{"s" if len(selected) != 1 else ""}.')
            st.rerun(scope='fragment')


# --- the table --------------------------------------------------------------------

def _new_run(presets_dir, runs, name='') -> dict:
    """A new run takes the preset and method of the last run, else the first preset."""
    if runs:
        last = runs[-1]
        return {'name': name, 'preset': last.get('preset'), 'method': copy.deepcopy(last['method'])}
    presets = fm.all_presets(presets_dir)
    first = next(iter(presets))
    return {'name': name, 'preset': first, 'method': copy.deepcopy(presets[first])}


def _rename(workflow_dir, index, wkey) -> None:
    runs = fm.load_runs(workflow_dir)
    if index < len(runs):
        runs[index]['name'] = st.session_state[wkey].strip()
        fm.save_runs(workflow_dir, runs)


def _sequence_now(raw_dir) -> dict:
    """The newest sequence in the raw file folder, read for the table while
    nothing runs (cached by its file and time)."""
    sld = flashida.newest_sequence(raw_dir) if raw_dir else None
    if sld is None:
        return {}
    try:
        stamp = (str(sld), sld.stat().st_mtime)
        cache = st.session_state.get('_ida_sld')
        if not cache or cache[0] != stamp:
            cache = (stamp, flashida.sequence_rows(sld))
            st.session_state['_ida_sld'] = cache
    except OSError:
        return {}
    return {'name': sld.name, 'modified': stamp[1], 'rows': cache[1]}


def _when(t) -> str:
    return time.strftime('%H:%M:%S', time.localtime(t)) if t else ''


_STATE = {
    'waiting': '⏳ Waiting for its raw file',
    'skipped': '⏭️ Raw file was there before the start',
    'done': '✅ Done',
    'failed': '❌ FLASHIda failed (see the log)',
}


def runs_table(workflow_dir, presets_dir, raw_dir: str, running: bool, footer=None) -> None:
    """The runs table. While running it refreshes itself every two seconds.

    ``footer`` (a function) is drawn below the table in the same fragment: the
    Run analysis bar, whose checks must see every change made in the table."""
    if OPEN_KEY in st.session_state:
        method_dialog(workflow_dir, presets_dir, st.session_state.pop(OPEN_KEY), read_only=running)
    dialog = st.session_state.pop(DIALOG_KEY, None)
    if dialog == 'import' and not running:
        import_dialog(workflow_dir, presets_dir)
    elif dialog == 'pattern' and not running:
        pattern_dialog(workflow_dir, presets_dir)
    if 'ida-added' in st.session_state:
        added, duplicates = st.session_state.pop('ida-added')
        st.toast(f'Added {added} run{"s" if added != 1 else ""}.'
                 + (f' {len(duplicates)} duplicate name{"s" if len(duplicates) != 1 else ""}, marked in the table.'
                    if duplicates else ''), icon='⚠️' if duplicates else '✅')

    @st.fragment(run_every=2 if running else None)
    def table():
        runs = fm.load_runs(workflow_dir)
        status = flashida.read_status(workflow_dir) if running else {}
        if running and status:
            sequence = status.get('sequence')
            planned = status.get('scheduled') or {}
            extra = [r['name'] + f" (row {r['row']})" for r in status.get('extra') or []]
        else:
            seq = _sequence_now(raw_dir)
            sequence = seq or None
            planned = flashida.scheduled(seq.get('rows', []), [r['name'] for r in runs]) if seq else {}
            names = {r['name'].casefold() for r in runs}
            extra = [r['name'] + f" (row {r['row']})" for r in seq.get('rows', [])
                     if r['name'].casefold() not in names] if seq else []

        if sequence:
            line = (f"Newest sequence **{sequence['name']}** (changed {_when(sequence.get('modified'))}): "
                    f"**{len(planned)} of {len(runs)} runs scheduled**.")
            if running:
                st.info('👀 Watching the raw file folder. ' + line + ' FLASHIda starts for a run as '
                        'soon as its raw file appears, also for a run the sequence does not list.')
            else:
                st.caption(line)
        elif raw_dir:
            (st.info if running else st.caption)(
                f'No sequence (.sld) in the raw file folder yet. ' +
                ('FLASHIda still starts for a run whose raw file appears.' if running else ''))

        presets = fm.all_presets(presets_dir)
        selected = [i for i in range(len(runs)) if st.session_state.get(f'ida-sel-{i}-{_version()}')]
        if not running and runs:
            _clipboard_bar(workflow_dir, runs, selected)
        widths = [0.35, 0.35, 4.4, 2.3, 1.7, 0.5, 1.8, 2.4, 0.5]
        head = st.columns(widths, vertical_alignment='center')
        for col, title in zip(head, ['', '', 'Raw file name', 'Preset', '', '', 'In sequence',
                                     'Status' if running else '', '']):
            col.markdown(f'<span style="font-size:.8rem;opacity:.65;font-weight:600">{title}</span>',
                         unsafe_allow_html=True)
        for i, run in enumerate(runs):
            row = st.columns(widths, vertical_alignment='center')
            row[0].checkbox('Select', key=f'ida-sel-{i}-{_version()}', label_visibility='collapsed',
                            disabled=running, help='Select for pasting a copied method')
            if row[5].button('📋', key=f'ida-copy-{i}', help='Copy the method of this run'):
                _copy_method(run)
                st.rerun(scope='fragment')
            # number, name, preset, method, in sequence, status, remove
            cols = [row[1], row[2], row[3], row[4], row[6], row[7], row[8]]
            cols[0].markdown(f'**{i + 1}**')
            wkey = f'ida-name-{i}-{_version()}'
            cols[1].text_input('Raw file name', value=run['name'], key=wkey, label_visibility='collapsed',
                               disabled=running, placeholder='Raw file name, without .raw',
                               on_change=_rename, args=(workflow_dir, i, wkey))
            problem = fm.name_problem(run['name'], [r['name'] for r in runs[:i]])
            preset = run.get('preset')
            edited = preset in presets and fm.changes(presets[preset], run['method'])
            text = preset if preset in presets else (f'{preset} (deleted)' if preset else 'Own method')
            badges = ' :orange-badge[✏️ changed]' if edited else ''
            if fm.check(run['method']):
                badges += ' :red-badge[⚠️ fails the method check]'
            elif fm.unknown_keys(run['method']):
                badges += ' :blue-badge[ℹ️ new keys]'
            cols[2].markdown(text + badges)
            if cols[3].button('⚙️ Method', key=f'ida-method-{i}', use_container_width=True):
                # Opened by the page, outside this fragment: the table refreshes
                # itself while running and would close a dialog of its own.
                st.session_state[OPEN_KEY] = i
                st.rerun()
            if problem and run['name']:
                cols[4].markdown(f':red-badge[{problem}]')
            elif not sequence:
                cols[4].markdown(':gray-badge[no sequence]')
            elif run['name'] in planned:
                cols[4].markdown(f':green-badge[✓ row {planned[run["name"]]}]')
            else:
                cols[4].markdown(':orange-badge[⚠️ not in sequence]',
                                 help='Not in the newest sequence. FLASHIda still starts if this raw file appears.')
            if running:
                state = (status.get('runs') or {}).get(run['name'], {})
                if state.get('state') == 'running':
                    minutes = int((time.time() - state.get('since', time.time())) // 60)
                    duration = fm.get(run['method'], 'global.duration')
                    cols[5].markdown(f'🟢 FLASHIda running · {minutes}'
                                     + (f' of {duration:g} min' if isinstance(duration, (int, float)) else ' min'))
                else:
                    cols[5].markdown(_STATE.get(state.get('state'), '⏳ Waiting for its raw file'))
            if cols[6].button('🗑️', key=f'ida-del-{i}', disabled=running, help='Remove this run'):
                del runs[i]
                fm.save_runs(workflow_dir, runs)
                _bump()
                st.rerun(scope='fragment')

        if not runs:
            st.caption('No runs yet. Add one per acquisition, named as the sequence names its raw file.')
        if extra:
            st.caption('Also in the sequence, not in the table (FLASHIda is not started for these): '
                       + ', '.join(f'`{e}`' for e in extra))

        if not running:
            a, b, c, d, e = st.columns([1.2, 2.4, 1.9, 1.2, 1.9])
            if c.button('🔢 Add runs by pattern', use_container_width=True,
                        help='Names made from prefix, date, sample, comment, replicate and suffix, '
                             'with a wash every few runs.'):
                st.session_state[DIALOG_KEY] = 'pattern'
                st.rerun()
            if d.button('📥 Import', use_container_width=True,
                        help='Runs from a sequence file exported here, or the raw file names of any '
                             'Xcalibur sequence (.sld).'):
                st.session_state[DIALOG_KEY] = 'import'
                st.rerun()
            e.download_button('⬇️ Export sequence', fs.export_runs(runs), use_container_width=True,
                              file_name=f'flashida-runs-{fs.today()}.json', mime='application/json',
                              disabled=not runs, key=f'ida-export-{_version()}',
                              help='Every run with its preset and method, as JSON. Import it again here.')
            if a.button('➕ Add run', use_container_width=True):
                runs.append(_new_run(presets_dir, runs))
                fm.save_runs(workflow_dir, runs)
                _bump()
                st.rerun(scope='fragment')
            seq_rows = _sequence_now(raw_dir).get('rows', []) if raw_dir else []
            have = {r['name'].casefold() for r in runs}
            new = [r['name'] for r in seq_rows if r['name'].casefold() not in have]
            if b.button('📋 Add names from the newest sequence', use_container_width=True,
                        disabled=not new,
                        help='Adds a run for each raw file the newest .sld in the raw file folder '
                             'writes that the table does not have yet.' if new else
                             'No sequence in the raw file folder, or the table has all its names.'):
                for name in new:
                    runs.append(_new_run(presets_dir, runs, name))
                fm.save_runs(workflow_dir, runs)
                _bump()
                st.rerun(scope='fragment')

        if footer is not None:
            footer()

    table()


def start_problems(workflow_dir, params: dict) -> list:
    """What keeps Run analysis from starting."""
    problems = []
    exe = str(params.get('executable') or '').strip()
    if not exe:
        problems.append('Enter the FLASHIda executable.')
    elif not (exe.lower().endswith('.exe') and Path(exe).is_file()):
        problems.append('The FLASHIda executable is not an existing .exe file.')
    raw = str(params.get('raw-files') or '').strip()
    if not raw or not Path(raw).is_dir():
        problems.append('Enter an existing raw file folder.')
    runs = fm.load_runs(workflow_dir)
    if not runs:
        problems.append('Add the runs to acquire.')
    names = []
    for i, run in enumerate(runs):
        problem = fm.name_problem(run['name'], names)
        if problem:
            problems.append(f'Run {i + 1}: {problem}')
        names.append(run['name'])
    return problems


def method_warnings(workflow_dir) -> list:
    """Runs whose method fails FLASHApp's check. They do not keep a run from
    starting: the user may override the check (Run anyway), and FLASHIda
    checks the method itself."""
    return [f"Run {i + 1}: its method fails FLASHApp's check, FLASHIda may refuse it (open ⚙️ Method)."
            for i, run in enumerate(fm.load_runs(workflow_dir)) if fm.check(run['method'])]
