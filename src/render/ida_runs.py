"""The runs table of the FLASHIda acquisition page and its Method dialog.

Each run is a raw file name and the method FLASHIda gets for it, a file in the
method folder (src/flashida_folder.py) picked in the run's dropdown. The
Method dialog makes a new method from any method in the folder or changes one,
by section, key by key, or as JSON. The built-in presets (FLASHIda's defaults,
those shipped in assets/flashida/presets.json and the user's own) are saved
into the folder once. While the acquisition runs, the table shows which runs
the newest sequence (.sld) schedules and what FLASHIda is doing for each.
"""

import json
import time
from pathlib import Path

import pandas as pd
import streamlit as st

from src import flashida
from src import flashida_folder as ff
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
OPEN_KEY = 'ida-open-method'  # the method the Method dialog opens with
DIALOG_KEY = 'ida-open-dialog'  # 'import', 'pattern' or 'method'
CLIPBOARD_KEY = 'ida-clipboard'  # {name, method_file} of the copied run


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


# --- method choices -----------------------------------------------------------------

def _short(value, limit: int = 24) -> str:
    text = json.dumps(value) if not isinstance(value, str) else (value or '""')
    return text if len(text) <= limit else text[:limit - 1] + '…'


def _folder_methods(folder) -> tuple:
    """(methods in the folder with the favorites first, favorites)."""
    methods = ff.list_methods(folder)
    favs = ff.favorites(folder)
    return ff.ordered(methods, favs), favs


def _method_label(favs, methods=None):
    def fmt(m):
        if m is None:
            return 'Pick a method'
        text = ('★ ' if m in favs else '') + m
        return text + ' (missing)' if methods is not None and m not in methods else text
    return fmt


def _default_method(methods: list, favs: list, runs: list):
    """The method a new run starts with: the last run's, else the first
    favorite, else the first built-in preset in the folder, else the first."""
    if runs and runs[-1].get('method_file') in methods:
        return runs[-1]['method_file']
    for m in list(favs) + ['LC-MS real-time ID, MS2 only', fm.DEFAULTS_PRESET]:
        if m in methods:
            return m
    return methods[0] if methods else None


# --- the method dialog -----------------------------------------------------------------

EDIT_KEY = 'ida-edit'  # the method being edited: {key, method, v}


@st.dialog('Method', width='large', on_dismiss='rerun')
def method_dialog(workflow_dir, folder) -> None:
    """New method (starting from any method in the folder) or edit one; saved
    as a file in the method folder."""
    methods, favs = _folder_methods(folder)
    if not methods:
        st.info('The method folder has no methods.')
        return
    mode = st.radio('Method', ['new', 'edit'], horizontal=True, key='ida-ed-mode',
                    label_visibility='collapsed',
                    format_func={'new': '➕ New method', 'edit': '✏️ Edit a method'}.get)
    a, b, c = st.columns([3, 3, 0.7], vertical_alignment='bottom')
    # Set once by the page (the picked run's method, or the one just saved)
    start = st.session_state.pop('ida-ed-start', None)
    if start in methods:
        st.session_state['ida-ed-src'] = start  # opened for a picked run, or after a save
    elif st.session_state.get('ida-ed-src') not in methods:
        st.session_state['ida-ed-src'] = methods[0]
    src = a.selectbox('Start from' if mode == 'new' else 'Method to edit', methods, key='ida-ed-src',
                      format_func=_method_label(favs),
                      help='Type to search. ★ marks favorites, listed first.')
    work = st.session_state.get(EDIT_KEY)
    if not work or work['key'] != (mode, src):
        try:
            method = ff.read(folder, src)
        except ValueError as e:
            st.error(str(e))
            return
        work = {'key': (mode, src), 'method': method, 'v': (work or {}).get('v', 0) + 1}
        st.session_state[EDIT_KEY] = work
    method, v = work['method'], work['v']
    if mode == 'new':
        new_name = b.text_input('Name of the new method', key='ida-ed-name', placeholder='Histone_MS2',
                                help='Saved as <name>.json in the method folder.').strip()
    else:
        new_name = src
        b.text_input('File', value=f'{src}.json', disabled=True, key=f'ida-ed-file-{v}')
    if c.button('★' if src in favs else '☆', key='ida-ed-fav', use_container_width=True,
                help=f"{'Remove' if src in favs else 'Add'} '{src}' {'from' if src in favs else 'to'} the favorites"):
        ff.toggle_favorite(folder, src)
        st.rerun(scope='fragment')

    keys = _keys()
    tabs = st.tabs(list(SECTIONS) + [f'All {len(keys)} keys', 'JSON'])
    for tab, (section, paths) in zip(tabs, SECTIONS.items()):
        with tab:
            cols = st.columns(2)
            for n, path in enumerate(paths):
                with cols[n % 2]:
                    key_widget(method, path, f'ida-w-{v}-{path}')
    with tabs[-2]:
        query = st.text_input('Find a key', key='ida-q',
                              placeholder='Part of a key or its description, e.g. faims, agc, exploration')
        q = query.strip().lower()
        matching = [p for p, k in keys.items()
                    if not q or q in p.lower() or q in k.description.lower() or q in label(p).lower()]
        if matching:
            pick = st.selectbox('Key', matching, key=f'ida-pick-{q}',
                                format_func=lambda p: f'{p}  ·  {label(p)}')
            key_widget(method, pick, f'ida-w-{v}-all-{pick}')
            rows = [{'Key': p, 'Value': _short(fm.get(method, p), 60) if fm.get(method, p, KeyError) is not KeyError else '',
                     'Default': _short(keys[p].default, 40) if keys[p].default is not None else keys[p].default_text,
                     'Unit': keys[p].unit} for p in matching]
            st.dataframe(pd.DataFrame(rows), hide_index=True, height=260, use_container_width=True)
        else:
            st.info('No key matches.')
    with tabs[-1]:
        text = st.text_area('Method as JSON', value=json.dumps(method, indent=2), height=380,
                            key=f'ida-json-{v}', help="Apply checks it against FLASHIda's method schema.")
        j1, j2, j3 = st.columns([1, 1, 1.2], vertical_alignment='bottom')
        if j1.button('Apply JSON', use_container_width=True, key=f'ida-json-apply-{v}'):
            try:
                parsed = json.loads(text)
            except ValueError as e:
                st.error(f'Not valid JSON: {e}')
            else:
                if not isinstance(parsed, dict):
                    st.error('A method is a JSON object ({...}).')
                else:
                    work.update(method=parsed, v=v + 1)
                    st.rerun(scope='fragment')
        j2.download_button('⬇️ Download', json.dumps(method, indent=2), use_container_width=True,
                           file_name=f'{new_name or src}.json', mime='application/json', key=f'ida-json-dl-{v}')
        with j3.popover('📥 Import .json file', use_container_width=True):
            upload = st.file_uploader('Method file', type=['json'], key=f'ida-json-up-{v}')
            if upload is not None:
                try:
                    parsed = json.loads(upload.getvalue().decode('utf-8-sig'))
                    assert isinstance(parsed, dict)
                except (ValueError, AssertionError):
                    st.error('This file is not a JSON method.')
                else:
                    work.update(method=parsed, v=v + 1)
                    st.rerun(scope='fragment')

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
    if mode == 'new':
        problem = ff.name_problem(new_name) if new_name else ''
        if not problem and new_name and ff.stem(new_name).casefold() in {m.casefold() for m in methods}:
            problem = f"'{ff.stem(new_name)}.json' is in the method folder already. Pick another name, or edit it."
        if problem:
            st.error(problem)
        save = st.button('💾 Save to method folder', type='primary', disabled=not new_name or bool(problem),
                         key='ida-ed-save')
    else:
        users = [r['name'] or f'run {i + 1}' for i, r in enumerate(fm.load_runs(workflow_dir))
                 if r.get('method_file') == src]
        st.caption(f"Saving overwrites `{src}.json`. "
                   + (f"{len(users)} run{'s use' if len(users) != 1 else ' uses'} it and get"
                      f"{'' if len(users) != 1 else 's'} the change when FLASHIda starts for "
                      f"{'them' if len(users) != 1 else 'it'}." if users else 'No run of the table uses it.'))
        save = st.button('💾 Save changes', type='primary', key='ida-ed-save')
    if save:
        name = ff.stem(new_name)
        try:
            ff.write(folder, name, method, overwrite=mode == 'edit')
        except (OSError, FileExistsError) as e:
            st.error(f'Could not save the method: {e}')
        else:
            # Widget values cannot be set once drawn; dropping them lets the
            # dialog start on the saved method the next time it opens.
            for key in (EDIT_KEY, 'ida-ed-name', 'ida-ed-src'):
                st.session_state.pop(key, None)
            st.session_state['ida-ed-start'] = name
            st.session_state['ida-saved'] = name
            _bump()
            st.rerun()


# --- import, export and runs by pattern ---------------------------------------------

_MODES = {'append': 'Append to the table', 'replace': 'Replace the table'}


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


def _method_pick(label_: str, methods: list, favs: list, key: str, default=None, **kw):
    index = methods.index(default) if default in methods else 0
    return st.selectbox(label_, methods, index=index, key=key, format_func=_method_label(favs),
                        help='Type to search. ★ marks favorites, listed first.', **kw)


@st.dialog('Import runs', width='large', on_dismiss='rerun')
def import_dialog(workflow_dir, folder) -> None:
    runs = fm.load_runs(workflow_dir)
    methods, favs = _folder_methods(folder)
    source = st.radio('Source', ['json', 'sld'], horizontal=True, key='ida-imp-source',
                      format_func=lambda s: {'json': 'Sequence file (JSON)',
                                             'sld': 'Xcalibur sequence (.sld)'}[s])
    upload = st.file_uploader('File', type=['json'] if source == 'json' else ['sld'],
                              key=f'ida-imp-file-{source}',
                              help='A file exported here with **Export sequence**.' if source == 'json' else
                              'Any Xcalibur sequence. Only the raw file names are taken.')
    mode = _mode_choice(runs, 'ida-imp-mode')
    pick = _method_pick('Method for the new runs', methods, favs, 'ida-imp-method',
                        _default_method(methods, favs, runs)) if source == 'sld' else None

    parsed = []
    if upload is not None:
        try:
            parsed = fs.parse_runs(upload.getvalue()) if source == 'json' else \
                fs.runs_from_names(fs.sld_names(upload.getvalue()), pick)
        except ValueError as e:
            st.error(str(e))
            parsed = []
        else:
            if not parsed:
                st.info('The file holds no runs.')
            else:
                def shown(r):
                    m = r.get('method_file')
                    if m and m in methods:
                        return m
                    if m:
                        return f'{m} (saved to the method folder)' if isinstance(r.get('method'), dict) \
                            else f'{m} (not in the method folder)'
                    return f"{r.get('preset') or r['name']} (saved to the method folder)"
                st.caption(f'{len(parsed)} run{"s" if len(parsed) != 1 else ""} in the file'
                           + ('.' if source == 'json' else
                              '. Only the names are taken; every run gets the method picked above.'))
                st.dataframe(pd.DataFrame({'Raw file name': [r['name'] for r in parsed],
                                           'Method': [shown(r) for r in parsed]}),
                             hide_index=True, use_container_width=True, height=min(36 + 35 * len(parsed), 250))
    combined, duplicates = fs.combine(runs, parsed, mode)
    if parsed:
        _duplicates_note(duplicates)
    if st.button(f'Import {len(parsed)} run{"s" if len(parsed) != 1 else ""}' if parsed else 'Import',
                 type='primary', disabled=not parsed, key='ida-imp-go'):
        new, written, missing = fs.into_folder(parsed, folder) if source == 'json' else (parsed, [], [])
        combined, duplicates = fs.combine(runs, new, mode)
        if written:
            st.session_state['ida-written'] = written
        if missing:
            st.session_state['ida-missing'] = missing
        _commit(workflow_dir, combined, len(new), duplicates)


@st.dialog('Add runs by pattern', width='large', on_dismiss='rerun')
def pattern_dialog(workflow_dir, folder) -> None:
    runs = fm.load_runs(workflow_dir)
    methods, favs = _folder_methods(folder)
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
    with right:
        pick = _method_pick('Method for the samples', methods, favs, 'ida-pat-method',
                            _default_method(methods, favs, runs))

    with st.container(border=True):
        w1, w2, w3, w4 = st.columns([0.9, 1, 2, 1.8], vertical_alignment='bottom')
        wash = w1.checkbox('Wash', key='ida-pat-wash')
        every = w2.number_input('Every … runs', min_value=1, max_value=99, value=3, step=1,
                                key='ida-pat-every', disabled=not wash)
        wash_name = w3.text_input('Wash name', value=fs.DEFAULT_WASH_NAME, key='ida-pat-wname',
                                  disabled=not wash, help='Takes {prefix}, {date}, {suffix} and {n}, the wash number.')
        with w4:
            wash_pick = _method_pick('Wash method', methods, favs, 'ida-pat-wmethod', fs.WASH_PRESET,
                                     disabled=not wash)

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
    if not methods:
        st.warning('Set the method folder under Instrument PC first.')
    new = [{'name': n['name'], 'method_file': wash_pick if n['wash'] else pick} for n in names]
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
            'Method': [f'{wash_pick} (wash run)' if n['wash'] else pick for n in names],
            '': ['⚠️ duplicate' if n['name'].casefold() in dup else '' for n in names],
        }), hide_index=True, use_container_width=True, height=min(36 + 35 * len(names), 300))
        _duplicates_note(duplicates)
    elif not error:
        st.caption('Enter at least one sample to see the runs.')
    if st.button(f'Add {len(new)} run{"s" if len(new) != 1 else ""}' if new else 'Add runs',
                 type='primary', disabled=not new or not methods, key='ida-pat-go'):
        _commit(workflow_dir, combined, len(new), duplicates)


def _clipboard_bar(workflow_dir, runs: list, selected: list) -> None:
    clip = st.session_state.get(CLIPBOARD_KEY)
    with st.container(border=True):
        a, b = st.columns([5, 1.9], vertical_alignment='center')
        if clip:
            a.markdown(f"📋 Copied: the method **{clip['method_file']}** of {clip['name'] or 'an unnamed run'}. "
                       'Tick runs, then paste.')
        else:
            a.caption('📋 No method copied. Copy one with the 📋 button of a run, tick the runs that '
                      'should get it, then paste.')
        if b.button(f'Paste into selection ({len(selected)})', use_container_width=True,
                    disabled=not clip or not selected, key='ida-paste'):
            for i in selected:
                runs[i] = {'name': runs[i].get('name', ''), 'method_file': clip['method_file']}
            fm.save_runs(workflow_dir, runs)
            _bump()
            st.toast(f'Pasted the method into {len(selected)} run{"s" if len(selected) != 1 else ""}.')
            st.rerun(scope='fragment')


# --- the table --------------------------------------------------------------------

def _rename(workflow_dir, index, wkey) -> None:
    runs = fm.load_runs(workflow_dir)
    if index < len(runs):
        runs[index]['name'] = st.session_state[wkey].strip()
        fm.save_runs(workflow_dir, runs)


def _set_method(workflow_dir, index, wkey) -> None:
    runs = fm.load_runs(workflow_dir)
    if index < len(runs) and st.session_state.get(wkey):
        runs[index] = {'name': runs[index].get('name', ''), 'method_file': st.session_state[wkey]}
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

# Raw file names are long: a smaller monospace font lets them fit their column.
_NAME_CSS = ("<style>[class*='st-key-ida-name-'] input{font-size:.7rem;"
             "font-family:'JetBrains Mono',ui-monospace,monospace;letter-spacing:-.01em;"
             "padding-left:.45rem;padding-right:.45rem}</style>")


def _prepare_folder(workflow_dir, presets_dir, folder) -> None:
    """Once per folder and session: the presets saved into the folder, and runs
    that still keep their method in the table moved to method files."""
    if st.session_state.get('_ida_prepared') == folder:
        return
    try:
        seeded = ff.seed_presets(folder, fm.all_presets(presets_dir))
        runs = fm.load_runs(workflow_dir)
        adopted, written = ff.adopt_embedded(runs, folder)
        if adopted != runs:
            fm.save_runs(workflow_dir, adopted)
            _bump()
    except OSError as e:
        st.error(f'Could not write to the method folder: {e}')
        return
    st.session_state['_ida_prepared'] = folder
    if seeded:
        st.toast(f'Saved {len(seeded)} built-in preset{"s" if len(seeded) != 1 else ""} to the method folder: '
                 + ', '.join(seeded))
    if written:
        st.toast(f'Moved the methods of the table to the method folder: {", ".join(written)}')


def _duration_line(runs, folder, running, status) -> str:
    if not runs:
        return ''
    if running and status:
        left, unknown = 0.0, 0
        for run in runs:
            state = (status.get('runs') or {}).get(run['name'], {})
            if state.get('state') in ('done', 'failed', 'skipped'):
                continue
            d = ff.duration(ff.run_method(run, folder))
            if d is None:
                unknown += 1
                continue
            if state.get('state') == 'running':
                d = max(0.0, d - (time.time() - state.get('since', time.time())) / 60)
            left += d
        text = f'⏱️ About **{ff.format_minutes(left)}** left · ends about **{ff.format_end(left)}**'
    else:
        total, unknown = ff.total_minutes(runs, folder)
        text = (f'⏱️ {len(runs)} run{"s" if len(runs) != 1 else ""} · **{ff.format_minutes(total)}** · '
                f'ends about **{ff.format_end(total)}** if started now')
    if unknown:
        text += f' ({unknown} run{"s" if unknown != 1 else ""} without a duration not counted)'
    return text


def runs_table(workflow_dir, presets_dir, raw_dir: str, running: bool, footer=None, folder: str = '') -> None:
    """The runs table. While running it refreshes itself every two seconds.

    ``footer`` (a function) is drawn below the table in the same fragment: the
    Run analysis bar, whose checks must see every change made in the table.
    ``folder`` is the method folder."""
    folder = (folder or '').strip()
    usable = ff.usable(folder)
    if usable and not running:
        _prepare_folder(workflow_dir, presets_dir, folder)
    dialog = st.session_state.pop(DIALOG_KEY, None)
    if not running and usable:
        if dialog == 'import':
            import_dialog(workflow_dir, folder)
        elif dialog == 'pattern':
            pattern_dialog(workflow_dir, folder)
        elif dialog == 'method':
            if OPEN_KEY in st.session_state:
                st.session_state['ida-ed-start'] = st.session_state.pop(OPEN_KEY)
            method_dialog(workflow_dir, folder)
    if 'ida-added' in st.session_state:
        added, duplicates = st.session_state.pop('ida-added')
        st.toast(f'Added {added} run{"s" if added != 1 else ""}.'
                 + (f' {len(duplicates)} duplicate name{"s" if len(duplicates) != 1 else ""}, marked in the table.'
                    if duplicates else ''), icon='⚠️' if duplicates else '✅')
    if 'ida-written' in st.session_state:
        st.toast('Saved to the method folder: ' + ', '.join(st.session_state.pop('ida-written')))
    if 'ida-missing' in st.session_state:
        st.toast('Not in the method folder (pick another method for these runs): '
                 + ', '.join(st.session_state.pop('ida-missing')), icon='⚠️')
    if 'ida-saved' in st.session_state:
        st.toast(f"Saved '{st.session_state.pop('ida-saved')}' to the method folder.", icon='💾')
    st.html(_NAME_CSS)

    @st.fragment(run_every=2 if running else None)
    def table():
        runs = fm.load_runs(workflow_dir)
        methods, favs = _folder_methods(folder)
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
        if not usable:
            st.info('📁 Set the **method folder** under Instrument PC: each run picks its method from there.')

        v = _version()
        selected = [i for i in range(len(runs)) if st.session_state.get(f'ida-sel-{i}-{v}')]
        if not running and runs:
            _clipboard_bar(workflow_dir, runs, selected)
        # select, #, name, method, favorite, copy, up, down, duplicate, sequence, [status], remove
        widths = [0.3, 0.35, 3.0, 3.6, 0.5, 0.5, 0.45, 0.45, 0.5, 0.9] + ([2.4] if running else []) + [0.5]
        titles = ['', '', 'Raw file name', 'Method', '', '', '', '', '', 'Seq.'] + (['Status'] if running else []) + ['']
        head = st.columns(widths, vertical_alignment='center')
        for col, title in zip(head, titles):
            col.markdown(f'<span style="font-size:.8rem;opacity:.65;font-weight:600">{title}</span>',
                         unsafe_allow_html=True)
        for i, run in enumerate(runs):
            c = st.columns(widths, vertical_alignment='center')
            c[0].checkbox('Select', key=f'ida-sel-{i}-{v}', label_visibility='collapsed',
                          disabled=running, help='Select for pasting a copied method')
            c[1].markdown(f'**{i + 1}**')
            wkey = f'ida-name-{i}-{v}'
            c[2].text_input('Raw file name', value=run['name'], key=wkey, label_visibility='collapsed',
                            disabled=running, placeholder='Raw file name, without .raw',
                            on_change=_rename, args=(workflow_dir, i, wkey), help=run['name'] or None)
            current = run.get('method_file')
            options = list(methods)
            if current and current not in options:
                options = [current] + options
            mkey = f'ida-m-{i}-{v}'
            c[3].selectbox('Method', options, index=options.index(current) if current in options else None,
                           key=mkey, label_visibility='collapsed', placeholder='Pick a method',
                           format_func=_method_label(favs, methods), disabled=running or not methods,
                           on_change=_set_method, args=(workflow_dir, i, mkey),
                           help=f'{current}.json in the method folder' if current else None)
            fav = current in favs
            if c[4].button('★' if fav else '☆', key=f'ida-fav-{i}', disabled=not current or current not in methods,
                           help=f"{'Remove' if fav else 'Add'} '{current}' {'from' if fav else 'to'} the "
                                'favorites, listed first in every method choice' if current else None):
                ff.toggle_favorite(folder, current)
                st.rerun(scope='fragment')
            if c[5].button('📋', key=f'ida-copy-{i}', disabled=not current,
                           help='Copy the method of this run'):
                st.session_state[CLIPBOARD_KEY] = {'name': run['name'], 'method_file': current}
                st.rerun(scope='fragment')
            for col, step, icon in ((c[6], -1, '↑'), (c[7], 1, '↓')):
                if col.button(icon, key=f'ida-mv{step}-{i}', disabled=running or not 0 <= i + step < len(runs),
                              help='Move up' if step < 0 else 'Move down'):
                    fm.save_runs(workflow_dir, ff.move(runs, i, step))
                    _bump()
                    st.rerun(scope='fragment')
            if c[8].button('⧉', key=f'ida-dup-{i}', disabled=running,
                           help='Duplicate this run below, with the next free name (R1 → R2)'):
                fm.save_runs(workflow_dir, ff.duplicate(runs, i))
                _bump()
                st.rerun(scope='fragment')

            problem = fm.name_problem(run['name'], [r['name'] for r in runs[:i]])
            method = ff.run_method(run, folder)
            if problem and run['name']:
                c[9].markdown(':red-badge[⚠️]', help=problem)
            elif current and usable and current not in methods:
                c[9].markdown(':red-badge[⚠️]', help=f"'{current}.json' is not in the method folder.")
            elif method is not None and fm.check(method):
                c[9].markdown(':orange-badge[⚠️]', help="The method fails FLASHApp's check; FLASHIda may refuse it.")
            elif not sequence:
                c[9].markdown(':gray-badge[–]', help='No sequence in the raw file folder.')
            elif run['name'] in planned:
                c[9].markdown(f':green-badge[✓ {planned[run["name"]]}]',
                              help=f'Row {planned[run["name"]]} of the newest sequence.')
            else:
                c[9].markdown(':orange-badge[–]',
                              help='Not in the newest sequence. FLASHIda still starts if this raw file appears.')
            if running:
                state = (status.get('runs') or {}).get(run['name'], {})
                if state.get('state') == 'running':
                    minutes = int((time.time() - state.get('since', time.time())) // 60)
                    d = ff.duration(method)
                    c[10].markdown(f'🟢 FLASHIda running · {minutes}' + (f' of {d:g} min' if d is not None else ' min'))
                else:
                    c[10].markdown(_STATE.get(state.get('state'), '⏳ Waiting for its raw file'))
            if c[-1].button('🗑️', key=f'ida-del-{i}', disabled=running, help='Remove this run'):
                del runs[i]
                fm.save_runs(workflow_dir, runs)
                _bump()
                st.rerun(scope='fragment')

        if not runs:
            st.caption('No runs yet. Add one per acquisition, named as the sequence names its raw file.')
        else:
            st.markdown(_duration_line(runs, folder, running, status))
        if extra:
            st.caption('Also in the sequence, not in the table (FLASHIda is not started for these): '
                       + ', '.join(f'`{e}`' for e in extra))

        if not running:
            a, b, c3, d, e, f, g = st.columns([1.1, 2.3, 1.8, 1.0, 1.25, 1.1, 1.9])
            if a.button('➕ Add run', use_container_width=True):
                runs.append({'name': '', 'method_file': _default_method(methods, favs, runs)})
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
                method_file = _default_method(methods, favs, runs)
                runs.extend({'name': n, 'method_file': method_file} for n in new)
                fm.save_runs(workflow_dir, runs)
                _bump()
                st.rerun(scope='fragment')
            if c3.button('🔢 Add runs by pattern', use_container_width=True, disabled=not usable,
                         help='Names made from prefix, date, sample, comment, replicate and suffix, '
                              'with a wash every few runs.'):
                st.session_state[DIALOG_KEY] = 'pattern'
                st.rerun()
            if d.button('📥 Import', use_container_width=True, disabled=not usable,
                        help='Runs from a sequence file exported here, or the raw file names of any '
                             'Xcalibur sequence (.sld).'):
                st.session_state[DIALOG_KEY] = 'import'
                st.rerun()
            e.download_button('⬇️ Export JSON', fs.export_runs(runs, folder), use_container_width=True,
                              file_name=f'flashida-runs-{fs.today()}.json', mime='application/json',
                              disabled=not runs, key=f'ida-export-{v}',
                              help='Every run with its method, as JSON. Import it again here; methods the '
                                   'method folder lacks are saved there.')
            f.download_button('⬇️ CSV', fs.export_csv(runs), use_container_width=True,
                              file_name=f'flashida-runs-{fs.today()}.csv', mime='text/csv',
                              disabled=not runs, key=f'ida-export-csv-{v}',
                              help='Row, raw file name and method of every run, as CSV.')
            if g.button('📝 New or edit method', use_container_width=True, disabled=not usable,
                        help='Make a new method from any method in the folder, or change one. '
                             'Saved as a file in the method folder.'):
                st.session_state[DIALOG_KEY] = 'method'
                if selected and runs[selected[0]].get('method_file'):
                    st.session_state[OPEN_KEY] = runs[selected[0]]['method_file']
                st.rerun()

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
    folder = str(params.get('method-folder') or '').strip()
    if not ff.usable(folder):
        problems.append('Enter an existing method folder.')
    runs = fm.load_runs(workflow_dir)
    if not runs:
        problems.append('Add the runs to acquire.')
    names = []
    for i, run in enumerate(runs):
        problem = fm.name_problem(run['name'], names)
        if problem:
            problems.append(f'Run {i + 1}: {problem}')
        names.append(run['name'])
    if ff.usable(folder):
        problems += [f'Run {i + 1}: {text}' for i, text in ff.method_problems(runs, folder)]
    return problems


def method_warnings(workflow_dir, folder: str = '') -> list:
    """Runs whose method fails FLASHApp's check. They do not keep a run from
    starting: the user may override the check (Run anyway), and FLASHIda
    checks the method itself."""
    out = []
    for i, run in enumerate(fm.load_runs(workflow_dir)):
        method = ff.run_method(run, folder)
        if method is not None and fm.check(method):
            out.append(f"Run {i + 1}: its method fails FLASHApp's check, FLASHIda may refuse it "
                       '(📝 New or edit method).')
    return out
