"""The runs table of the FLASHIda acquisition page, and the Method editor page.

Each run is a raw file name and the method FLASHIda gets for it, a file in the
method folder (src/flashida_folder.py) picked in the run's dropdown. The
Method editor page makes a new method from any method in the folder or changes
one, by step of the acquisition (the tabs and parts of src/render/ida_sections.py),
key by key, or in a JSON tree (streamlit-extras' json_editor). The
built-in presets (FLASHIda's defaults, the methods shipped as one file each in
assets/flashida/methods, and the user's own) are saved into the folder once. While the acquisition runs, the table shows which runs the newest
sequence (.sld) schedules and what FLASHIda is doing for each.
"""

import copy
import hashlib
import html
import random
import json
import time
from pathlib import Path

import pandas as pd
import streamlit as st
from streamlit_extras.json_editor import json_editor

from src import flashida
from src import flashida_compare as fc
from src import flashida_folder as ff
from src import flashida_method as fm
from src import flashida_sequence as fs
from src.render import ida_sections as sections
from src.render.run_order import run_order

_VERSION = 'ida-method-version'
DIALOG_KEY = 'ida-open-dialog'  # 'import' or 'pattern'
CLIPBOARD_KEY = 'ida-clipboard'  # {name, method_file} of the copied run


def _keys() -> dict:
    if '_ida_keys' not in st.session_state:
        st.session_state['_ida_keys'] = fm.keys()
    return st.session_state['_ida_keys']


def label(path: str, short: bool = False) -> str:
    """The key's name with its unit; `short` leaves out the scan or sweep level."""
    key = _keys().get(path)
    text = sections.name(path, short)
    if key is not None and key.unit:
        text += f' ({key.unit})'
    return text


def _help(key: fm.Key, default: str) -> str:
    text = f'**Default:** {default}\n\n`{key.path}`\n\n{key.description}'
    if key.enum_help:
        text += '\n\n' + '\n\n'.join(f'**{v or "(none)"}**: {h}' for v, h in zip(key.enum, key.enum_help))
    return text


def _version() -> int:
    return st.session_state.get(_VERSION, 0)


def _bump() -> None:
    st.session_state[_VERSION] = _version() + 1


def key_widget(method: dict, path: str, wkey: str, disabled: bool = False,
               short: bool = False, why_not: str | None = None) -> bool:
    """One widget for one key of a method; True when it changed the method.
    A cleared field removes the key, so FLASHIda uses its own default.
    `why_not` says why FLASHIda does not read the key now; it is then drawn
    greyed out."""
    key = _keys().get(path)
    if key is None:
        return False
    present = fm.get(method, path, missing=KeyError) is not KeyError
    value = fm.get(method, path) if present else key.default
    default = sections.default_of(method, path)
    help_ = _help(key, default)
    if why_not:
        help_ = f'**Not read now.** {why_not}\n\n' + help_
        disabled = True
    name = label(path, short)
    choices = list(key.enum) if key.enum else sections.choices(path, method)
    marked = sections.default_value(path)  # marked '(default)' in a list of values

    if choices and key.type == 'array':
        # A list of fixed values (activations to sweep, additional scans): pick any of them
        current = [v for v in (value or [])]
        options = list(choices) + [v for v in current if v not in choices]
        new = st.multiselect(name, options, default=current, key=wkey, help=help_, disabled=disabled,
                             placeholder=f'default: {default}')
        if new == current:
            return False
        if not new and key.default is None:
            new = None
    elif choices:
        options = list(choices)
        if value is not None and value not in options:
            options.append(value)
        # A part switched on and off elsewhere (MS3) has no empty choice here
        if None not in options and '' not in options and (key.default is None or not key.enum) \
                and path not in sections.CLEARS:
            options = [None] + options
        empty = f'default: {default}'
        # Not set shows as the empty choice ('' where the schema lists one), never as a value
        shown = '' if value is None and '' in options else value
        new = st.selectbox(
            name, options, index=options.index(shown) if shown in options else 0, key=wkey,
            help=help_, disabled=disabled,
            format_func=lambda v: (empty if v in (None, '') else
                                   f'{v}  (default)' if v == marked else str(v)),
        )
        if new == '' and not present:
            return False
    elif key.type == 'object' or (key.type == 'array' and key.items_type == 'object'):
        text = json.dumps(value, indent=1) if value not in (None, {}, []) else ''
        raw = st.text_area(name, value=text, key=wkey, disabled=disabled, height=110,
                           help=help_ + '\n\nJSON; leave empty to remove the key.',
                           placeholder=f'default: {default}')
        if raw == text:
            return False
        if not raw.strip():
            new = None
        else:
            try:
                new = json.loads(raw)
            except ValueError as e:
                st.error(f'{name}: not valid JSON ({e}).')
                return False
    elif key.type == 'boolean':
        new = st.toggle(name, value=bool(value), key=wkey, help=help_, disabled=disabled)
    elif key.type in ('integer', 'number'):
        is_int = key.type == 'integer'
        shown = value if isinstance(value, (int, float)) and not isinstance(value, bool) else None
        if shown is not None:
            shown = int(shown) if is_int else float(shown)
        new = st.number_input(
            name, value=shown, step=1 if is_int else None, key=wkey, help=help_,
            disabled=disabled, placeholder=f'default: {default}',
            format=None if is_int else '%g',
        )
        if new is not None and not is_int and float(new).is_integer() and isinstance(value, int):
            new = int(new)
    elif key.type == 'array':
        text = ', '.join(str(v) for v in (value or []))
        raw = st.text_input(name, value=text, key=wkey, help=help_ + '\n\nComma-separated.',
                            disabled=disabled, placeholder=f'default: {default}')
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
                            placeholder=f'default: {default}')
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


# --- the Method editor page --------------------------------------------------------------

EDIT_KEY = 'ida-edit'  # the method being edited: {key, method, v}
EDITOR_PAGE = 'content/FLASHIda/FLASHIdaMethodEditor.py'


def open_in_editor(method_file=None) -> None:
    """Switch to the Method editor page, on this method of the folder."""
    if method_file:
        st.session_state['ida-ed-start'] = method_file
        st.session_state['ida-ed-mode'] = 'edit'
    st.switch_page(EDITOR_PAGE)


def _toggle(work: dict, part: str, wkey: str) -> None:
    sections.switch(work['method'], part, st.session_state[wkey])
    work['v'] += 1  # fresh widgets for the part's keys


def _part(work: dict, part: str, role: str, have: set, found: dict) -> None:
    """One part of a tab: its name, an On switch when it is optional, what it
    does and sums up to now, its problems, then its keys by group."""
    method, v = work['method'], work['v']
    title = sections.TITLES.get(part, sections.PARTS[part][0])
    icon, what, groups, short = (sections.PARTS[part][i] for i in (1, 5, 3, 4))
    on = part in have
    if role == 'with' and not on:
        return
    with st.container(border=role != 'switch'):
        head, right = st.columns([3, 1.4], vertical_alignment='center')
        head.markdown(f'#### {icon} {title}' if role == 'switch' else f'**{title}**' + ('' if on else ' · off'))
        if role in ('option', 'switch'):
            reason = None if on else sections.why_not_add(part, have)
            wkey = f'ida-on-{v}-{part}'
            right.toggle('On' if role == 'option' else sections.SWITCH_TEXT[part], value=on, key=wkey,
                         disabled=bool(reason), help=reason or what,
                         on_change=_toggle, args=(work, part, wkey))
        st.caption(f'{what}  \n**Now:** {sections.line(method, part)}' if on and part != 'global' else what)
        for text in found.get(part, []):
            st.warning(text, icon='⚠️')
        if not on:
            if role == 'switch':
                st.info('Off: FLASHIda skips this step. Switch it on to set it up.', icon='⏸')
            return
        changed = False
        for gtitle, paths in groups.items():
            st.markdown(f'**{gtitle}**')
            reasons = {sections.inactive(method, p) for p in paths}
            if len(reasons) == 1 and None not in reasons:
                st.caption(f'⏸ {reasons.pop()}')
            cols = st.columns(2)
            for n, path in enumerate(paths):
                with cols[n % 2]:
                    changed |= key_widget(method, path, f'ida-w-{v}-{path}', short=short,
                                          why_not=sections.inactive(method, path))
        if changed:
            st.rerun()  # the summaries and the line above the tabs follow the change


def _status(have: set, found: dict) -> str:
    """One line above the tabs: the tabs with something to look at, and the
    steps that are off. (In the tab names it would change them, and st.tabs
    then starts again on its first tab.)"""
    warn, off = [], []
    for name, (icon, parts) in sections.TAB_PARTS.items():
        main, role = parts[0]
        if role == 'switch' and main not in have:
            off.append(name)
        n = sum(len(found.get(p, [])) for p, _ in parts if p in have)
        if n:
            warn.append(f'{icon} **{name}** {n}')
    text = ('⚠️ ' + ' · '.join(warn)) if warn else '✅ Nothing FLASHIda would refuse'
    return text + (f"  —  off: {', '.join(off)}" if off else '')


SHIPPED = 'shipped:'  # a Start from choice that is a method FLASHApp ships, not a folder file
EMPTY = 'empty:'  # a Start from choice: a method with no keys at all


def _start_label(favs):
    folder_label = _method_label(favs)

    def fmt(m):
        if m == EMPTY:
            return '∅ Empty'
        if isinstance(m, str) and m.startswith(SHIPPED):
            return f'📦 {m[len(SHIPPED):]}'
        return folder_label(m)
    return fmt


# --- the Compare mode of the Method editor -------------------------------------------------

def open_compare(method_a: str, method_b: str) -> None:
    """Switch to the Method editor comparing two methods of the folder."""
    st.session_state.update({'ida-ed-mode': 'compare', 'ida-cmp-a-start': method_a,
                             'ida-cmp-b-start': method_b})
    st.switch_page(EDITOR_PAGE)


def _to_compare(src: str) -> None:
    """From Edit: compare the working copy (B) with the saved file (A)."""
    work = st.session_state.get(EDIT_KEY)
    if work:
        work['key'] = ('compare', src)
    st.session_state.update({'ida-ed-mode': 'compare', 'ida-cmp-a-start': src, 'ida-cmp-b-start': src})


def _to_edit(src: str) -> None:
    """From Compare: edit B further, the changes taken over kept."""
    work = st.session_state.get(EDIT_KEY)
    if work:
        work['key'] = ('edit', src)
    st.session_state.update({'ida-ed-mode': 'edit', 'ida-ed-start': src})


def _shown(value) -> str:
    return '–' if value is None else _short(value, 40)


def _compare(workflow_dir, folder, methods: list, favs: list) -> None:
    """A (any method) next to B (a method of the folder): the settings in which
    B acquires differently, taken over from A one by one, by step or all at
    once, into B's working copy, which is then saved like an edited method."""
    shipped = fm.builtin_presets()
    a_choices = ([m for m in methods if m in favs] + [SHIPPED + n for n in shipped]
                 + [m for m in methods if m not in favs])
    for key, choices, first in (('ida-cmp-a', a_choices, 0), ('ida-cmp-b', methods, 1)):
        start = st.session_state.pop(f'{key}-start', None)
        if start in choices:
            st.session_state[key] = start
        elif st.session_state.get(key) not in choices:
            st.session_state[key] = choices[min(first, len(choices) - 1)]
    c1, c2 = st.columns(2)
    a_src = c1.selectbox('A: take settings from', a_choices, key='ida-cmp-a', format_func=_start_label(favs),
                         help='Any method of the folder, or one FLASHApp ships (📦). It is only read.')
    b_src = c2.selectbox('B: method to change', methods, key='ida-cmp-b', format_func=_method_label(favs),
                         help='A method of the folder. Settings taken over from A change its working copy; '
                              'nothing is written before you save.')
    try:
        a = copy.deepcopy(shipped[a_src[len(SHIPPED):]]) if a_src.startswith(SHIPPED) else ff.read(folder, a_src)
    except (ValueError, KeyError) as e:
        st.error(str(e))
        return
    a, notes = fc.migrate(a)
    work = st.session_state.get(EDIT_KEY)
    if not work or work['key'] != ('compare', b_src):
        try:
            b = ff.read(folder, b_src)
        except ValueError as e:
            st.error(str(e))
            return
        work = {'key': ('compare', b_src), 'method': b, 'v': (work or {}).get('v', 0) + 1,
                'saved': copy.deepcopy(b)}
        st.session_state[EDIT_KEY] = work
    b, v = work['method'], work['v']
    for note in notes:
        st.info(note, icon='🔁')

    def take(new: dict) -> None:
        work.update(method=new, v=v + 1)
        st.rerun()

    diff = fm.changes(work['saved'], b)
    rows = fc.compare(a, b)
    with st.container():  # always there, so what follows keeps its place
        if diff:
            d1, d2 = st.columns([5, 1.2], vertical_alignment='center')
            d1.caption(f'✏️ B not saved yet: {len(diff)} key{"s" if len(diff) != 1 else ""} taken over or changed.')
            if d2.button('↺ Reset B', use_container_width=True, key=f'ida-cmp-reset-{v}',
                         help=f'Back to {b_src}.json as it is saved'):
                take(copy.deepcopy(work['saved']))
        if not rows:
            st.success('B acquires exactly as A.' if a_src != b_src or diff else
                       'Pick two methods; B acquires exactly as A.', icon='✅')

    if rows:
        st.markdown('**Settings that differ**')
        table = pd.DataFrame({
            'A → B': [False] * len(rows),
            'Step': [r.tab for r in rows],
            'Setting': [fc.switch_label(r.path) or (r.path if r.tab == fc.UNKNOWN else label(r.path))
                        for r in rows],
            'A': [_shown(r.a) for r in rows],
            'B': [_shown(r.b) for r in rows],
        })
        picked = st.data_editor(table, hide_index=True, use_container_width=True, key=f'ida-cmp-table-{v}',
                                disabled=['Step', 'Setting', 'A', 'B'],
                                height=min(38 + 35 * len(rows), 460),
                                column_config={'A → B': st.column_config.CheckboxColumn(width='small')})
        chosen = [rows[i].path for i, on in enumerate(picked['A → B']) if on]
        t1, t2, _ = st.columns([1.6, 1.6, 3])
        if t1.button(f'Ticked A → B ({len(chosen)})', type='primary', disabled=not chosen,
                     use_container_width=True, key=f'ida-cmp-take-{v}'):
            take(fc.adopt(b, a, chosen))
        if t2.button(f'All A → B ({len(rows)})', use_container_width=True, key=f'ida-cmp-all-{v}',
                     help='B then acquires exactly as A'):
            take(fc.adopt(b, a, [r.path for r in rows]))

    with st.expander('All settings'):
        every = fc.all_settings(a, b)
        st.dataframe(pd.DataFrame({
            '': ['≠' if d else '' for *_, d in every],
            'Step': [tab for tab, *_ in every],
            'Setting': [label(p) for _, p, *_ in every],
            'A': [_shown(x) for _, _, x, _, _ in every],
            'B': [_shown(y) for _, _, _, y, _ in every],
        }), hide_index=True, use_container_width=True, height=420,
            column_config={'': st.column_config.TextColumn(width='small')})
        st.caption('The values FLASHIda works with: a setting left out shows its default. '
                   '≠ marks the settings that differ above.')

    found = sections.hints(b)
    if found:
        st.caption('⚠️ B: ' + ' · '.join(f'{sections.PARTS[p][0] if p in sections.PARTS else p}: {t}'
                                          for p, ts in found.items() for t in ts[:1]))

    with st.container(border=True):
        users = [r['name'] or f'run {i + 1}' for i, r in enumerate(fm.load_runs(workflow_dir))
                 if r.get('method_file') == b_src]
        s1, s2, s3 = st.columns([3.4, 1.5, 1.5], vertical_alignment='center')
        s1.caption(f"**Save changes** overwrites `{b_src}.json`; "
                   + (f"{len(users)} run{'s use' if len(users) != 1 else ' uses'} it. "
                      if users else 'no run of the table uses it. ')
                   + '**Save as new method** leaves it as it is.')
        s2.button('✏️ Edit B', use_container_width=True, key=f'ida-cmp-edit-{v}', on_click=_to_edit,
                  args=(b_src,), help='Edit B further in the tabs; what was taken over stays')
        if s3.button('💾 Save changes', type='primary', use_container_width=True, disabled=not diff,
                     key=f'ida-cmp-save-{v}'):
            try:
                ff.write(folder, ff.stem(b_src), b, overwrite=True)
            except OSError as e:
                st.error(f'Could not save the method: {e}')
            else:
                st.session_state.pop(EDIT_KEY, None)
                st.session_state['ida-saved'] = b_src
                st.rerun()
        n1, n2 = st.columns([3.4, 3.0], vertical_alignment='bottom')
        new_name = n1.text_input('Name of the new method', key=f'ida-cmp-new-{v}', placeholder=f'{b_src}_2',
                                 help='Saved as <name>.json in the method folder; B stays as it is.').strip()
        problem = ff.name_problem(new_name) if new_name else ''
        if not problem and new_name and ff.stem(new_name).casefold() in {m.casefold() for m in methods}:
            problem = f"'{ff.stem(new_name)}.json' is in the method folder already. Pick another name."
        if problem:
            st.error(problem)
        if n2.button('💾 Save as new method', use_container_width=True, disabled=not new_name or bool(problem),
                     key=f'ida-cmp-save-new-{v}'):
            name = ff.stem(new_name)
            try:
                ff.write(folder, name, b, overwrite=False)
            except (OSError, FileExistsError) as e:
                st.error(f'Could not save the method: {e}')
            else:
                # The new method is B from now on; the old B stays as it was saved
                st.session_state.pop(EDIT_KEY, None)
                st.session_state['ida-cmp-b-start'] = name
                st.session_state['ida-saved'] = name
                st.rerun()


def method_editor(workflow_dir, folder) -> None:
    """The Method editor page: a new method (starting from any method in the
    folder) or a change to one, by step of the acquisition; saved as a file in
    the method folder."""
    methods, favs = _folder_methods(folder)
    if not methods:
        st.info('The method folder has no methods.')
        return
    if 'ida-saved' in st.session_state:
        st.toast(f"Saved '{st.session_state.pop('ida-saved')}' to the method folder.", icon='💾')
    mode = st.radio('Method', ['new', 'edit', 'compare'], horizontal=True, key='ida-ed-mode',
                    label_visibility='collapsed',
                    format_func={'new': '➕ New method', 'edit': '✏️ Edit a method',
                                 'compare': '⇄ Compare'}.get)
    if mode == 'compare':
        _compare(workflow_dir, folder, methods, favs)
        return
    a, b, c = st.columns([3, 3, 0.7], vertical_alignment='bottom')
    # A new method starts from a favorite, a method FLASHApp ships (as it is in
    # assets/flashida/methods, so a git pull updates it) or any method in the folder
    shipped = fm.builtin_presets()
    choices = methods
    if mode == 'new':
        choices = ([EMPTY] + [m for m in methods if m in favs] + [SHIPPED + n for n in shipped]
                   + [m for m in methods if m not in favs])
    # Set once: the picked run's method (from the Acquisition page), or the one just saved
    start = st.session_state.pop('ida-ed-start', None)
    if start in choices:
        st.session_state['ida-ed-src'] = start
    elif st.session_state.get('ida-ed-src') not in choices:
        st.session_state['ida-ed-src'] = next((c for c in choices if c != EMPTY), choices[0])
    src = a.selectbox('Start from' if mode == 'new' else 'Method to edit', choices, key='ida-ed-src',
                      format_func=_start_label(favs),
                      help='Type to search. ∅ Empty starts with no settings at all; ★ marks favorites; '
                           '📦 the methods FLASHApp ships (as they are in the app, also when the folder keeps '
                           'an older copy).')
    work = st.session_state.get(EDIT_KEY)
    if not work or work['key'] != (mode, src):
        try:
            if src == EMPTY:
                method = {}
            elif src.startswith(SHIPPED):
                method = copy.deepcopy(shipped[src[len(SHIPPED):]])
            else:
                method = ff.read(folder, src)
        except (ValueError, KeyError) as e:
            st.error(str(e))
            return
        work = {'key': (mode, src), 'method': method, 'v': (work or {}).get('v', 0) + 1,
                'saved': copy.deepcopy(method)}
        st.session_state[EDIT_KEY] = work
    method, v = work['method'], work['v']
    if mode == 'new':
        new_name = b.text_input('Name of the new method', key='ida-ed-name', placeholder='Histone_MS2',
                                help='Saved as <name>.json in the method folder.').strip()
    else:
        new_name = src
        b.text_input('File', value=f'{src}.json', disabled=True, key=f'ida-ed-file-{v}')
    in_folder = src in methods
    if c.button('★' if src in favs else '☆', key='ida-ed-fav', use_container_width=True, disabled=not in_folder,
                help=(f"{'Remove' if src in favs else 'Add'} '{src}' {'from' if src in favs else 'to'} the favorites"
                      if in_folder else 'Favorites are methods of the method folder.')):
        ff.toggle_favorite(folder, src)
        st.rerun()

    found = sections.hints(method)
    have = sections.present(method) | {'global'}
    keys = _keys()
    diff = fm.changes(work['saved'], method)
    # st.tabs starts again on its first tab when its names change, or when an
    # element above it comes or goes: so the names are fixed, and what changes
    # with the method sits in one container that is always there.
    with st.container():
        if diff:
            shown = ', '.join(f'{label(p).split(" (")[0].lower()} {_short(x)} → {_short(y)}'
                              for p, x, y in diff[:4])
            more = f' and {len(diff) - 4} more' if len(diff) > 4 else ''
            n1, n2 = st.columns([5, 1.6], vertical_alignment='center')
            n1.caption(f'✏️ Not saved yet, changed from **{_start_label(favs)(src)}**: {shown}{more}')
            if mode == 'edit':
                n2.button('⇄ Compare with saved', use_container_width=True, key=f'ida-ed-cmp-{v}',
                          on_click=_to_compare, args=(src,),
                          help='Every setting that now acquires differently from the saved file')
        st.markdown(_status(have, found))
    tabs = st.tabs([f'{icon} {n}' for n, (icon, _) in sections.TAB_PARTS.items()]
                   + [f'All {len(keys)} keys', 'JSON'])
    for tab, (name_, (icon, parts)) in zip(tabs, sections.TAB_PARTS.items()):
        with tab:
            tab_on = parts[0][0] in have
            if name_ == 'Global':
                st.caption(f'**Now:** {sections.run_line(method)}')
            for part, role in parts:
                if role == 'option' and not tab_on:
                    continue
                _part(work, part, role, have, found)
    with tabs[-2]:
        query = st.text_input('Find a key', key='ida-q',
                              placeholder='Part of a key or its description, e.g. faims, agc, exploration')
        q = query.strip().lower()
        matching = [p for p, k in keys.items()
                    if not q or q in p.lower() or q in k.description.lower() or q in label(p).lower()]
        if matching:
            pick = st.selectbox('Key', matching, key=f'ida-pick-{q}',
                                format_func=lambda p: f'{p}  ·  {label(p)}')
            if key_widget(method, pick, f'ida-w-{v}-all-{pick}', why_not=sections.inactive(method, pick)):
                st.rerun()
            rows = [{'Key': p, 'Value': _short(fm.get(method, p), 60) if fm.get(method, p, KeyError) is not KeyError else '',
                     'Default': sections.default_of(method, p), 'Unit': keys[p].unit} for p in matching]
            st.dataframe(pd.DataFrame(rows), hide_index=True, height=260, use_container_width=True)
        else:
            st.info('No key matches.')
    with tabs[-1]:
        st.caption('The method as JSON. Click a value to change it, ➕ on an object to add a '
                   'key, 🗑 to remove one; a change is kept at once and checked against '
                   "FLASHIda's method schema.")
        tree_key = f"ida-json-{v}-{hashlib.sha1(json.dumps(method, sort_keys=True).encode()).hexdigest()[:12]}"
        edited = json_editor(method, name='method', collapsed=1, editable=True,
                             display_data_types=False, sort_keys=False, key=tree_key)['data']
        # The tree returns the method as it was when it was drawn; the other tabs
        # may have added keys since. So a change is what differs from the tree's
        # own previous value, and only that is applied to the method.
        before = st.session_state.get(f'{tree_key}-before')
        st.session_state[f'{tree_key}-before'] = copy.deepcopy(edited)
        if before is not None and edited != before:
            try:
                new = apply_tree_edit(method, before, edited)
            except ValueError as e:
                st.error(str(e))
            else:
                work.update(method=new, v=v + 1)
                st.rerun()
        j2, j3 = st.columns([1, 1.2], vertical_alignment='bottom')
        j2.download_button('⬇️ Download', json.dumps(method, indent=2), use_container_width=True,
                           file_name=f'{new_name or src.removeprefix(SHIPPED)}.json', mime='application/json', key=f'ida-json-dl-{v}')
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
                    st.rerun()

    if fm.check(method):
        st.caption('⚠️ in a tab: this FLASHApp expects FLASHIda to refuse the method there. The check uses the '
                   'method schema bundled with FLASHApp, which may be older than your FLASHIda; **Run anyway** '
                   'on the Acquisition page starts the run with the method as it is.')
    unknown = fm.unknown_keys(method)
    if unknown:
        with st.container(border=True):
            st.markdown(f'ℹ️ **{len(unknown)} key{"s" if len(unknown) != 1 else ""} this FLASHApp does not know** '
                        '(a newer FLASHIda may). They are kept and passed to FLASHIda as they are, which checks '
                        'them when the run starts; edit or remove them in the JSON tab.')
            st.markdown(', '.join(f'`{k}`' for k in unknown[:20]) + (' …' if len(unknown) > 20 else ''))

    with st.container(border=True):
        if mode == 'new':
            problem = ff.name_problem(new_name) if new_name else ''
            if not problem and new_name and ff.stem(new_name).casefold() in {m.casefold() for m in methods}:
                problem = f"'{ff.stem(new_name)}.json' is in the method folder already. Pick another name, or edit it."
            if problem:
                st.error(problem)
            s1, s2 = st.columns([4, 1.6], vertical_alignment='center')
            s1.caption('Saves the method as a new file in the method folder.' if new_name
                       else 'Name the new method above to save it.')
            save = s2.button('💾 Save to method folder', type='primary', disabled=not new_name or bool(problem),
                             key='ida-ed-save', use_container_width=True)
        else:
            users = [r['name'] or f'run {i + 1}' for i, r in enumerate(fm.load_runs(workflow_dir))
                     if r.get('method_file') == src]
            s1, s2 = st.columns([4, 1.6], vertical_alignment='center')
            s1.caption(f"Saving overwrites `{src}.json`. "
                       + (f"{len(users)} run{'s use' if len(users) != 1 else ' uses'} it and get"
                          f"{'' if len(users) != 1 else 's'} the change when FLASHIda starts for "
                          f"{'them' if len(users) != 1 else 'it'}." if users else 'No run of the table uses it.'))
            save = s2.button('💾 Save changes', type='primary', key='ida-ed-save', disabled=not diff,
                             use_container_width=True)
    if save:
        name = ff.stem(new_name)
        try:
            ff.write(folder, name, method, overwrite=mode == 'edit')
        except (OSError, FileExistsError) as e:
            st.error(f'Could not save the method: {e}')
        else:
            # Widget values cannot be set once drawn; dropping them lets the
            # page start on the saved method, to be edited further.
            for key in (EDIT_KEY, 'ida-ed-name', 'ida-ed-src', 'ida-ed-mode'):
                st.session_state.pop(key, None)
            st.session_state['ida-ed-start'] = name
            st.session_state['ida-ed-mode'] = 'edit'
            st.session_state['ida-saved'] = name
            _bump()
            st.rerun()


def apply_tree_edit(method: dict, before, after) -> dict:
    """The method with what changed in the JSON tree from `before` to `after`:
    values set or added, keys removed. Keys the tree did not show (added to the
    method after it was drawn) stay. A method is a JSON object, so a root that
    is not one raises ValueError."""
    if not isinstance(after, dict):
        raise ValueError('A method is a JSON object ({...}).')

    def merge(target, old, new):
        out = copy.deepcopy(target) if isinstance(target, dict) else {}
        old = old if isinstance(old, dict) else {}
        for key in old.keys() - new.keys():
            out.pop(key, None)
        for key, value in new.items():
            if key in old and old[key] == value:
                continue
            if isinstance(value, dict) and isinstance(old.get(key), dict):
                out[key] = merge(out.get(key), old[key], value)
            else:
                out[key] = copy.deepcopy(value)
        return out

    return merge(method, before, after)


# --- import, export and runs by pattern ---------------------------------------------

_MODES = {'append': 'Append to the table', 'replace': 'Replace the table'}


def _mode_choice(runs: list, key: str) -> str:
    return st.radio('Into the table', list(_MODES), key=key, horizontal=True,
                    format_func=lambda m: _MODES[m] if m == 'append' or not runs
                    else f'{_MODES[m]} ({len(runs)} run{"s" if len(runs) != 1 else ""} now)')


def _duplicates_note(duplicates: list, before: list = ()) -> None:
    """Names the table would have more than once: those these runs bring in,
    and apart from them those the table has twice already (``before``)."""
    def listed(names):
        return ', '.join(f'`{d}`' for d in names[:10]) + (' …' if len(names) > 10 else '')
    old = {d.casefold() for d in before}
    new = [d for d in duplicates if d.casefold() not in old]
    kept = [d for d in duplicates if d.casefold() in old]
    if new:
        st.warning(f'{len(new)} of these name{"s are" if len(new) != 1 else " is"} in the table already: '
                   + listed(new) + '. The runs are added and marked; rename or remove them before the '
                   'run can start.')
    if kept:
        st.info(f'The table already has {len(kept)} name{"s" if len(kept) != 1 else ""} more than once, '
                'not caused by these runs: ' + listed(kept) + '. Rename or remove one of each before the '
                'run can start.')


def _commit(workflow_dir, runs: list, added: int, duplicates: list) -> None:
    fm.save_runs(workflow_dir, runs)
    _bump()
    st.session_state['ida-added'] = (added, duplicates)
    st.rerun()


def _method_pick(label_: str, methods: list, favs: list, key: str, default=None, **kw):
    index = methods.index(default) if default in methods else 0
    return st.selectbox(label_, methods, index=index, key=key, format_func=_method_label(favs),
                        help='Type to search. ★ marks favorites, listed first.', **kw)


@st.dialog('Import sequence', width='large', on_dismiss='rerun')
def import_dialog(workflow_dir, folder) -> None:
    runs = fm.load_runs(workflow_dir)
    methods, favs = _folder_methods(folder)
    upload = st.file_uploader('Sequence', type=['json', 'sld', 'csv'], key='ida-imp-file',
                              help='A sequence exported here with **Export sequence** (.json: names and '
                                   'methods; .csv: the table), or an Xcalibur sequence (.sld, or .csv exported '
                                   'from Xcalibur): only the raw file names are taken from those.')
    # The file says what it is: FLASHApp's sequence (.json), Xcalibur's (.sld), or a CSV of either
    suffix = Path(upload.name).suffix.lower() if upload is not None else '.json'
    source = {'.sld': 'sld', '.csv': 'csv'}.get(suffix, 'json')
    csv_parsed, csv_methods = [], False
    if source == 'csv':
        try:
            csv_parsed, csv_methods = fs.csv_runs(upload.getvalue())
        except ValueError as e:
            st.error(str(e))
    names_only = source == 'sld' or (source == 'csv' and not csv_methods)
    mode = _mode_choice(runs, 'ida-imp-mode')
    pick = _method_pick('Method for the new runs', methods, favs, 'ida-imp-method',
                        _default_method(methods, favs, runs)) if names_only else None

    parsed = []
    if upload is not None:
        try:
            if source == 'json':
                parsed = fs.parse_runs(upload.getvalue())
            elif source == 'sld':
                parsed = fs.runs_from_names(fs.sld_names(upload.getvalue()), pick)
            else:
                parsed = (csv_parsed if csv_methods
                          else fs.runs_from_names([r['name'] for r in csv_parsed], pick))
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
                           + ('. Only the names are taken; every run gets the method picked above.'
                              if names_only else '.'))
                st.dataframe(pd.DataFrame({'Raw file name': [r['name'] for r in parsed],
                                           'Method': [shown(r) for r in parsed]}),
                             hide_index=True, use_container_width=True, height=min(36 + 35 * len(parsed), 250))
    combined, duplicates = fs.combine(runs, parsed, mode)
    if parsed:
        _duplicates_note(duplicates, fs.duplicate_names(runs) if mode == 'append' else [])
    if st.button(f'Import {len(parsed)} run{"s" if len(parsed) != 1 else ""}' if parsed else 'Import',
                 type='primary', disabled=not parsed, key='ida-imp-go'):
        if source == 'json':
            new, written, missing = fs.into_folder(parsed, folder)
        else:  # a table exported here names its methods; those the folder lacks are reported
            new, written = parsed, []
            missing = sorted({r['method_file'] for r in parsed if r.get('method_file') not in methods} - {None})
        combined, duplicates = fs.combine(runs, new, mode)
        if written:
            st.session_state['ida-written'] = written
        if missing:
            st.session_state['ida-missing'] = missing
        _commit(workflow_dir, combined, len(new), duplicates)


@st.dialog('Reorder runs', width='large', on_dismiss='rerun')
def order_dialog(workflow_dir, folder) -> None:
    """The runs as a list to drag into a new order, shuffle, or give washes;
    Save writes the table."""
    runs = fm.load_runs(workflow_dir)
    methods, favs = _folder_methods(folder)
    inserts = []
    if methods:
        st.caption('Runs with these methods are put between the samples: 🎲 Randomize leaves them in '
                   'place, and Insert adds or removes them. — leaves a kind out.')
        cols = st.columns(len(fs.INSERT_KINDS))
        for col, (kind, (label_, _)) in zip(cols, fs.INSERT_KINDS.items()):
            options = [None] + methods
            default = fs.WASH_PRESET if kind == 'wash' and fs.WASH_PRESET in methods else None
            method = col.selectbox(f'{label_} method', options, index=options.index(default),
                                   key=f'ida-order-{kind}', format_func=lambda m: '—' if m is None
                                   else _method_label(favs)(m),
                                   help=f'The method of the {label_.lower()} runs. Type to search.')
            if method:
                inserts.append({'kind': kind, 'label': label_, 'method': method,
                                'base': fs.join([fs.today(), {'qc': 'QC'}.get(kind, kind)]),
                                'every': 3 if kind == 'wash' else 0})
    used = [i['method'] for i in inserts]
    if len(set(used)) != len(used):
        st.error('Two kinds have the same method: give each its own, or leave one out (—).')
        return
    # By the kind a run was inserted as, not by its method alone: a sample run
    # that uses the wash method is still a sample run
    in_use = {i['kind']: i['method'] for i in inserts}
    kinds = [fs.insert_kind(r, in_use) for r in runs]
    order = run_order([{'name': r.get('name'), 'method': r.get('method_file') or 'no method', 'kind': k}
                       for r, k in zip(runs, kinds)],
                      key='ida-run-order', inserts=inserts)
    if order is not None:
        # Every inserted run keeps its kind from now on, also one found by its name
        new = [dict(runs[e], **({'kind': kinds[e]} if kinds[e] else {})) if isinstance(e, int)
               else {'name': e['name'], 'method_file': e['method'], 'kind': e['kind']} for e in order]
        fm.save_runs(workflow_dir, new)
        _bump()
        added = sum(not isinstance(e, int) for e in order)
        kept = {e for e in order if isinstance(e, int)}
        removed = sum(i not in kept for i in range(len(runs)))
        st.toast('Saved the new order of the runs'
                 + (f', {added} run{"s" if added != 1 else ""} added' if added else '')
                 + (f', {removed} run{"s" if removed != 1 else ""} removed' if removed else '') + '.',
                 icon='↕️')
        st.rerun()


def _kind_word(kind: str) -> str:
    """'wash', 'QC standard', 'blank': a kind of inserted run inside a sentence."""
    label_ = fs.INSERT_KINDS[kind][0]
    return label_ if label_.startswith('QC') else label_.lower()


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
    # A random order stays the same while the dialog is edited; 🎲 draws a new one
    seed = st.session_state.setdefault('ida-pat-seed', random.randrange(1 << 30))
    if order == 'random' and right.button('🎲 Shuffle again', key='ida-pat-shuffle',
                                          help='A new random order. Washes, QC standards and blanks '
                                               'keep their places between the sample runs.'):
        st.session_state['ida-pat-seed'] = seed = random.randrange(1 << 30)
    with right:
        pick = _method_pick('Method for the samples', methods, favs, 'ida-pat-method',
                            _default_method(methods, favs, runs))

    with st.container(border=True):
        st.caption('Runs between the samples: each kind after every few sample runs, none after the '
                   'last. Names take {prefix}, {date}, {suffix} and {n}, the number of that kind.')
        inserts, picks = [], {}
        for kind, (label_, name_default) in fs.INSERT_KINDS.items():
            # The wash keeps its former keys, so a dialog opened before keeps its values
            k = {'wash': ('ida-pat-wash', 'ida-pat-every', 'ida-pat-wname', 'ida-pat-wmethod')}.get(
                kind, tuple(f'ida-pat-{kind}-{x}' for x in ('on', 'every', 'name', 'method')))
            w1, w2, w3, w4 = st.columns([1.1, 1, 2, 1.8], vertical_alignment='bottom')
            on = w1.checkbox(label_, key=k[0])
            every = w2.number_input('Every … runs', min_value=1, max_value=99,
                                    value=3 if kind == 'wash' else 10, step=1, key=k[1], disabled=not on)
            name = w3.text_input(f'{label_} name', value=name_default, key=k[2], disabled=not on)
            with w4:
                picks[kind] = _method_pick(f'{label_} method', methods, favs, k[3],
                                           fs.WASH_PRESET if kind == 'wash' else None, disabled=not on)
            if on:
                inserts.append({'kind': kind, 'every': int(every), 'name': name})

    mode = _mode_choice(runs, 'ida-pat-mode')

    names, error = [], ''
    try:
        samples = fs.parse_samples(samples_text, comment.strip(), int(replicates))
    except ValueError as e:
        error = str(e)
    else:
        unknown = [(fs.INSERT_KINDS[i['kind']][0], p) for i in inserts for p in fs.unknown_placeholders(i['name'])]
        if unknown:
            error = 'Unknown placeholder: ' + ', '.join(f'{{{p}}} in the {lab.lower()} name' for lab, p in unknown)
        else:
            names = fs.pattern_names(samples, prefix.strip(), date.strip(), suffix.strip(), rfmt, order,
                                     inserts=inserts, seed=seed,
                                     taken=[r.get('name') or '' for r in runs] if mode == 'append' else ())
    if error:
        st.error(error)
    if not methods:
        st.warning('Set the method folder under Instrument PC first.')
    new = [{'name': n['name'], 'method_file': pick} if n['kind'] == 'sample'
           else {'name': n['name'], 'method_file': picks[n['kind']], 'kind': n['kind']} for n in names]
    combined, duplicates = fs.combine(runs, new, mode)
    dup = {d.casefold() for d in duplicates}
    if names:
        start = len(runs) if mode == 'append' else 0
        counts = {kind: sum(n['kind'] == kind for n in names) for kind in ('sample', *fs.INSERT_KINDS)}
        parts = [f"{counts[kind]} {_kind_word(kind)}{('es' if kind == 'wash' else 's') if counts[kind] != 1 else ''}"
                 for kind in fs.INSERT_KINDS if counts[kind]]
        st.markdown(f"**Preview:** {counts['sample']} sample run{'s' if counts['sample'] != 1 else ''}"
                    + (', ' + ', '.join(parts) if parts else ''))
        st.dataframe(pd.DataFrame({
            '#': range(start + 1, start + len(names) + 1),
            'Raw file name': [n['name'] for n in names],
            'Method': [pick if n['kind'] == 'sample'
                       else f"{picks[n['kind']]} ({_kind_word(n['kind'])})" for n in names],
            '': ['⚠️ duplicate' if n['name'].casefold() in dup else '' for n in names],
        }), hide_index=True, use_container_width=True, height=min(36 + 35 * len(names), 300))
        _duplicates_note(duplicates, fs.duplicate_names(runs) if mode == 'append' else [])
    elif not error:
        st.caption('Enter at least one sample to see the runs.')
    if st.button(f'Add {len(new)} run{"s" if len(new) != 1 else ""}' if new else 'Add runs',
                 type='primary', disabled=not new or not methods, key='ida-pat-go'):
        _commit(workflow_dir, combined, len(new), duplicates)


def _clipboard_bar(workflow_dir, runs: list, selected: list) -> None:
    clip = st.session_state.get(CLIPBOARD_KEY)
    with st.container(border=True):
        a, c, d, b = st.columns([3.2, 1.3, 1.3, 1.9], vertical_alignment='center')
        if clip:
            a.markdown(f"📋 Copied: the method **{clip['method_file']}** of {clip['name'] or 'an unnamed run'}. "
                       'Tick runs, then paste.')
        else:
            a.caption('📋 No method copied. Tick one run and copy its method, tick the runs that '
                      'should get it, then paste.')
        one = runs[selected[0]] if len(selected) == 1 else None
        if c.button('📋 Copy method', use_container_width=True, key='ida-copy',
                    disabled=not (one and one.get('method_file')),
                    help='Copies the method of the ticked run' if one else 'Tick one run to copy its method'):
            st.session_state[CLIPBOARD_KEY] = {'name': one['name'], 'method_file': one['method_file']}
            _bump()  # untick it, ready to tick the runs that get the method
            st.rerun(scope='fragment')
        if d.button(f'⧉ Duplicate ({len(selected)})', use_container_width=True, key='ida-dup',
                    disabled=not selected,
                    help='Adds a copy of each ticked run right below it, with the next free name (R1 → R2)'):
            for i in sorted(selected, reverse=True):
                runs = ff.duplicate(runs, i)
            fm.save_runs(workflow_dir, runs)
            _bump()
            st.rerun(scope='fragment')
        if len(selected) == 2 and all(runs[i].get('method_file') for i in selected):
            ma, mb = (runs[i]['method_file'] for i in selected)
            if b.button('⇄ Compare their methods', use_container_width=True, key='ida-compare-runs',
                        disabled=ma == mb, help=f'{ma} (A) next to {mb} (B) in the Method editor'
                        if ma != mb else 'Both runs use the same method.'):
                open_compare(ma, mb)
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

def _failed_notice(status: dict) -> None:
    """Runs FLASHIda failed for, at the top of the table, and a toast the
    first time this page sees each one: a long acquisition hides a ❌ row."""
    failed = [name for name, s in (status.get('runs') or {}).items() if s.get('state') == 'failed']
    if not failed:
        return
    st.error(f"**FLASHIda failed for {len(failed)} run{'s' if len(failed) != 1 else ''}:** "
             + ', '.join(f'`{n}`' for n in failed) + '. The log below says why; the other runs go on.',
             icon='❌')
    seen = st.session_state.setdefault('ida-failed-seen', set())
    started = status.get('started')
    for name in failed:
        if (started, name) not in seen:
            seen.add((started, name))
            st.toast(f"FLASHIda failed for '{name}'.", icon='❌')


# Raw file names are long: a smaller monospace font lets them fit their column.
_NAME_CSS = ("<style>[class*='st-key-ida-name-'] input{font-size:.8rem;"
             "font-family:'JetBrains Mono',ui-monospace,monospace;letter-spacing:-.01em;"
             "padding-left:.45rem;padding-right:.45rem}</style>")

# Badges of the Seq. column, each with its own tooltip shown on hover. The
# colors are mid tones, readable on the light and the dark theme.
_BADGE_CSS = """<style>
.ida-badges{display:flex;gap:.3rem;align-items:center;flex-wrap:nowrap}
.ida-badge{position:relative;display:inline-block;padding:.05rem .4rem;border-radius:.4rem;
  font-size:.95rem;font-weight:600;line-height:1.4;cursor:help;white-space:nowrap}
.ida-badge.green{background:rgba(33,195,84,.16);color:#1f9d4c}
.ida-badge.orange{background:rgba(255,164,33,.18);color:#d97706}
.ida-badge.red{background:rgba(255,43,43,.14);color:#e53935}
.ida-badge.gray{background:rgba(128,132,149,.18);color:#808495}
.ida-badge .tip{visibility:hidden;opacity:0;transition:opacity .12s;position:absolute;z-index:1000;
  right:0;top:calc(100% + 6px);width:22rem;padding:.55rem .7rem;border-radius:.5rem;
  background:#262730;color:#fafafa;font-size:.8rem;font-weight:400;line-height:1.45;
  white-space:pre-line;text-align:left;box-shadow:0 4px 14px rgba(0,0,0,.25);pointer-events:none}
.ida-badge:hover .tip{visibility:visible;opacity:1}
.ida-badge .tip b{color:#fff}
</style>"""


def _badge(color: str, text: str, title: str, body: str) -> str:
    esc = html.escape
    return (f'<span class="ida-badge {color}">{esc(text)}'
            f'<span class="tip"><b>{esc(title)}</b>\n{esc(body)}</span></span>')


def _seq_badges(run, method_file, method, problem, sequence, planned, present,
                missing_file, folder) -> str:
    """The Seq. column of a run: its place in the newest sequence, always,
    and beside it what keeps the run from being acquired as planned."""
    name = run.get('name') or ''
    raw = f'{name}.raw' if name else 'its raw file'
    with_method = f" with the method '{method_file}'" if method_file else ''
    if not sequence:
        out = [_badge('gray', '–', 'No sequence yet',
                      'There is no Xcalibur sequence (.sld) in the raw file folder, so FLASHApp cannot '
                      'tell which runs the instrument will acquire.\n'
                      f'FLASHIda still starts{with_method} as soon as {raw} appears in the folder.')]
    elif name in planned:
        out = [_badge('green', f'✓ {planned[name]}', f"Row {planned[name]} of the newest sequence",
                      f"The newest sequence '{sequence['name']}' writes {raw} in row {planned[name]}.\n"
                      f'When the file appears in the raw file folder, FLASHIda starts{with_method}.')]
    else:
        out = [_badge('orange', '–', 'Not in the newest sequence',
                      f"The newest sequence '{sequence['name']}' has no row that writes {raw}, so the "
                      'instrument will probably not acquire this run.\n'
                      'Check that the name matches the sequence (case does not matter). If the file '
                      f'appears anyway, FLASHIda still starts{with_method}.')]
    if name and name.casefold() in present:
        out.append(_badge('gray', '⏭️', 'Raw file already there: skipped',
                          f'{raw} is already in the raw file folder. FLASHIda starts only for raw files '
                          'that appear after the acquisition starts, so it is not started for this run.\n'
                          'To acquire it, give the run (and its sequence row) a new name, or move the '
                          'old file out of the folder.'))
    if problem and name:
        out.append(_badge('red', '⚠️', 'Raw file name not usable',
                          f'{problem}\nFLASHIda is not started for this run until the name is fixed; '
                          'Run analysis stays off.'))
    elif missing_file:
        out.append(_badge('red', '⚠️', 'Method file missing',
                          f"'{method_file}.json' is not in the method folder '{folder}'. It was renamed, "
                          'moved or deleted.\nPick another method for this run, or put the file back.'))
    elif method is not None:
        issues = fm.check(method)
        if issues:
            listed = '\n'.join(f'• {p}' for p in issues[:6])
            more = f'\n… and {len(issues) - 6} more' if len(issues) > 6 else ''
            out.append(_badge('orange', '⚠️', "Method fails FLASHApp's check",
                              f"'{method_file or 'The method'}' does not fit FLASHIda's method schema:\n"
                              f'{listed}{more}\n'
                              'FLASHIda checks the method itself when the run starts and may refuse it. '
                              'Fix it with 📝 New or edit method, or start with ⚠️ Run anyway.'))
    return _BADGE_CSS + '<div class="ida-badges">' + ''.join(out) + '</div>'


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
    if not running and dialog == 'order':
        order_dialog(workflow_dir, folder)
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
        if running:
            _failed_notice(status)
        if extra:
            shown = extra[:8] + ([f'and {len(extra) - 8} more'] if len(extra) > 8 else [])
            many = len(extra) != 1
            them = 'them' if many else 'it'
            st.warning(f"**{len(extra)} raw file{'s' if many else ''} of the sequence "
                       f"{'are' if many else 'is'} not in the table**, so FLASHIda is not started "
                       f"for {them}: " + ', '.join(f'`{e}`' for e in shown)
                       + (f'. Stop the acquisition to add {them}.' if running else
                          f'. Add {them} with **📋 Add from sequence** below the table.'),
                       icon='⚠️')

        # Raw files already in the folder: the start leaves their runs alone
        present = set() if running or not raw_dir else {n.casefold() for n in flashida.raw_names(raw_dir)}

        v = _version()
        fav_stamp = hashlib.md5('\n'.join(sorted(favs)).encode()).hexdigest()[:8]
        selected = [i for i in range(len(runs)) if st.session_state.get(f'ida-sel-{i}-{v}')]
        if not running and runs:
            _clipboard_bar(workflow_dir, runs, selected)
        # select, #, name, method, favorite, sequence, [status], remove
        widths = [0.3, 0.35, 3.6, 4.1, 0.5, 1.2] + ([2.4] if running else []) + [0.5]
        titles = ['', '', 'Raw file name', 'Method', '', 'Sequence'] + (['Status'] if running else []) + ['']
        tips = {'Raw file name': 'The name of the raw file the instrument writes for this run, without .raw. '
                                 'FLASHIda starts for the run when this file appears in the raw file folder.',
                'Method': 'The FLASHIda method (a .json file in the method folder) the run is acquired with.',
                'Sequence': "The run's row in the newest Xcalibur sequence (.sld) of the raw file folder, "
                            'and beside it anything that keeps the run from being acquired as planned. '
                            'Hover a badge for details.'}
        head = st.columns(widths, vertical_alignment='center')
        for col, title in zip(head, titles):
            if not title:
                continue
            col.markdown(f'<span title="{html.escape(tips.get(title, ""))}" style="font-size:.95rem;'
                         f'opacity:.8;font-weight:700;cursor:{"help" if title in tips else "default"}">'
                         f'{title}</span>', unsafe_allow_html=True)
        if len(runs) > 1 and not running:
            if head[1].button('↕', key='ida-reorder', type='tertiary',
                              help='**Reorder runs**: drag the runs into a new order in a list, shuffle '
                                   'them (🎲 Randomize) or put a wash after every few runs. The order is '
                                   'the order of the exported sequence and of the duration estimate.'):
                st.session_state[DIALOG_KEY] = 'order'
                st.rerun()
        if runs and not running:
            def tick_all(keys=[f'ida-sel-{i}-{v}' for i in range(len(runs))], all_key=f'ida-sel-all-{v}'):
                for k in keys:
                    st.session_state[k] = st.session_state[all_key]
            # Follows the rows, so ticking the last run by hand ticks it too
            st.session_state[f'ida-sel-all-{v}'] = len(selected) == len(runs)
            head[0].checkbox('Select all', key=f'ida-sel-all-{v}',
                             label_visibility='collapsed', on_change=tick_all,
                             help='Tick or untick every run')
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
            # The favorites in the key: a new star redraws the labels (★ name)
            mkey = f'ida-m-{i}-{v}-{fav_stamp}'
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

            problem = fm.name_problem(run['name'], [r['name'] for r in runs[:i]])
            method = ff.run_method(run, folder)
            c[5].html(_seq_badges(run, current, method, problem, sequence, planned, present,
                                  usable and current not in methods, folder))
            if running:
                state = (status.get('runs') or {}).get(run['name'], {})
                if state.get('state') == 'running':
                    minutes = int((time.time() - state.get('since', time.time())) // 60)
                    d = ff.duration(method)
                    c[6].markdown(f'🟢 FLASHIda running · {minutes}' + (f' of {d:g} min' if d is not None else ' min'))
                else:
                    c[6].markdown(_STATE.get(state.get('state'), '⏳ Waiting for its raw file'))
            if c[-1].button('🗑️', key=f'ida-del-{i}', disabled=running, help='Remove this run'):
                del runs[i]
                fm.save_runs(workflow_dir, runs)
                _bump()
                st.rerun(scope='fragment')

        if not runs:
            st.caption('No runs yet. Add one per acquisition, named as the sequence names its raw file.')
        else:
            st.markdown(_duration_line(runs, folder, running, status))

        if not running:
            a, b, c3, d, e, g = st.columns([1.1, 2.0, 1.8, 1.4, 1.4, 1.9])
            if a.button('➕ Add run', use_container_width=True):
                runs.append({'name': '', 'method_file': _default_method(methods, favs, runs)})
                fm.save_runs(workflow_dir, runs)
                _bump()
                st.rerun(scope='fragment')
            seq_rows = _sequence_now(raw_dir).get('rows', []) if raw_dir else []
            have = {r['name'].casefold() for r in runs}
            new = [r['name'] for r in seq_rows if r['name'].casefold() not in have]
            seq_name = (_sequence_now(raw_dir) or {}).get('name') if raw_dir else None
            if b.button('📋 Add from sequence', use_container_width=True, disabled=not new,
                        help=(f"**Adds the {len(new)} raw file{'s' if len(new) != 1 else ''} of the Xcalibur "
                              f"sequence '{seq_name}' that the table does not have yet** "
                              f"({', '.join(new[:5])}{', …' if len(new) > 5 else ''}).\n\n"
                              'The sequence is the newest .sld in the raw file folder. The new runs get the '
                              'method of the last run (else the first favorite); pick another per run '
                              'afterwards. Runs already in the table stay as they are.') if new else
                             ('**Adds the raw files of the Xcalibur sequence the table does not have yet.**'
                              '\n\n' + (f"The table already has every raw file of '{seq_name}'." if seq_name
                                        else 'There is no sequence (.sld) in the raw file folder yet.'))):
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
            if d.button('📥 Import sequence', use_container_width=True, disabled=not usable,
                        help='Runs from a sequence exported here (.json or the .csv table), or the raw '
                             'file names of an Xcalibur sequence (.sld, or .csv exported from Xcalibur).'):
                st.session_state[DIALOG_KEY] = 'import'
                st.rerun()
            with e.popover('⬇️ Export sequence', use_container_width=True, disabled=not runs):
                st.download_button('FLASHApp sequence (.json)', fs.export_runs(runs, folder),
                                   use_container_width=True, file_name=f'flashida-runs-{fs.today()}.json',
                                   mime='application/json', key=f'ida-export-{v}',
                                   help='Every run with its method. Import it again here with **Import '
                                        'sequence**; methods the method folder lacks are saved there.')
                st.markdown('**For Xcalibur**')
                inst = st.text_input('Instrument method', key='ida-xcal-meth',
                                     placeholder=r'C:\Xcalibur\methods\TopDown',
                                     help='The Xcalibur instrument method every sample runs with, as Xcalibur '
                                          'writes it. Left empty, pick it in Xcalibur after the import.')
                st.download_button('Xcalibur sequence (.csv)',
                                   fs.export_xcalibur_csv(runs, raw_dir or '', inst.strip()),
                                   use_container_width=True, file_name=f'flashida-sequence-{fs.today()}.csv',
                                   mime='text/csv', key=f'ida-export-xcal-{v}',
                                   help='One sample per run, its raw file name as File Name, the raw file '
                                        'folder as Path. Import it in Xcalibur (Sequence Setup: File > Import); '
                                        'Import sequence here reads it too.')
                st.divider()
                st.download_button('Table (.csv)', fs.export_csv(runs), use_container_width=True,
                                   file_name=f'flashida-runs-{fs.today()}.csv', mime='text/csv',
                                   key=f'ida-export-csv-{v}',
                                   help='Row, raw file name and method of every run, e.g. for Excel. '
                                        'Import sequence reads it back with the methods.')
            if g.button('📝 New or edit method', use_container_width=True, disabled=not usable,
                        help='Opens the Method editor: make a new method from any method in the folder, '
                             'or change one (the ticked run\'s, if one is ticked). Saved as a file in the '
                             'method folder.'):
                open_in_editor(runs[selected[0]].get('method_file') if selected else None)

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


def etd_notice(workflow_dir, folder: str = '') -> None:
    """A warning that stands out when a run's method acquires with ETD: the
    reagent is set up in Tune, which no method does."""
    runs = []
    for i, run in enumerate(fm.load_runs(workflow_dir)):
        scans = sections.etd_scans(ff.run_method(run, folder))
        if scans:
            runs.append(f"Run {i + 1} `{run['name'] or 'unnamed'}` "
                        f"({run.get('method_file') or 'its method'}: {', '.join(scans)})")
    if not runs:
        return
    shown = runs[:5] + ([f'and {len(runs) - 5} more'] if len(runs) > 5 else [])
    st.warning('**ETD: set up ETD in Tune before the acquisition starts.** '
               f"{'This method uses' if len(runs) == 1 else 'These methods use'} ETD, "
               'which FLASHIda cannot switch on itself.\n\n'
               + '\n'.join(f'- {r}' for r in shown), icon='⚡')


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
