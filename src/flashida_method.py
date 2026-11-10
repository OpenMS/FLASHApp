"""FLASHIda methods: the keys a method file may hold, their defaults, checks
and the presets the acquisition page offers. Streamlit-free, so the workflow
process can check a method before FLASHIda gets it.

The keys come from FLASHIda's own JSON schema, ``assets/flashida/method.schema.json``
(``docs/site/method.schema.json`` of FlashIDA, the file its key reference is
generated from). Copy a newer one over it when FLASHIda gains keys.
"""

import copy
import json
import re
from dataclasses import dataclass, field
from pathlib import Path

ASSETS = Path(__file__).resolve().parents[1] / 'assets' / 'flashida'
SCHEMA_FILE = ASSETS / 'method.schema.json'
BUILTIN_PRESETS_FILE = ASSETS / 'presets.json'
# Presets the user saves, kept once per computer (the instrument PC) in the
# local workspaces folder, so every workspace offers them.
USER_PRESETS_FILE = '.flashida-presets.json'
DEFAULTS_PRESET = 'FLASHIda defaults'


@dataclass
class Key:
    """One leaf key of a method, as the page shows it."""
    path: str
    type: str
    description: str = ''
    default: object = None
    default_text: str = ''
    unit: str = ''
    enum: list = field(default_factory=list)
    enum_help: list = field(default_factory=list)
    min_items: int = 0
    max_items: int | None = None
    items_type: str = ''


_schema_cache = {}


def schema(path=SCHEMA_FILE) -> dict:
    path = Path(path)
    if path not in _schema_cache:
        _schema_cache[path] = json.loads(path.read_text(encoding='utf-8'))
    return _schema_cache[path]


def _resolve(node: dict, root: dict) -> dict:
    ref = node.get('$ref')
    if not ref:
        return node
    target = root
    for part in ref.lstrip('#/').split('/'):
        target = target[part]
    merged = dict(target)
    merged.update({k: v for k, v in node.items() if k != '$ref'})
    return merged


def _type(node: dict) -> str:
    t = node.get('type', '')
    if isinstance(t, list):
        t = next((x for x in t if x != 'null'), '')
    if not t and 'oneOf' in node:
        t = next((o.get('type') for o in node['oneOf'] if o.get('type') not in (None, 'null')), '')
    return t or 'string'


def keys(root: dict | None = None) -> dict:
    """{dotted path: Key} for every leaf key of a method (objects are walked,
    arrays of objects are one key)."""
    root = root or schema()
    out = {}

    def walk(props: dict, prefix: str, allowed=None):
        for name, node in props.items():
            if allowed is not None and name not in allowed:
                continue
            node = _resolve(node, root)
            path = prefix + name
            if node.get('properties') and _type(node) == 'object':
                names = (node.get('propertyNames') or {}).get('enum')
                walk(node['properties'], path + '.', set(names) if names else None)
                continue
            items = _resolve(node.get('items') or {}, root)
            out[path] = Key(
                path=path, type=_type(node), description=node.get('description', ''),
                default=copy.deepcopy(node.get('default')),
                default_text=node.get('x-default-text', ''), unit=node.get('x-unit', ''),
                enum=list(node.get('enum') or []), enum_help=list(node.get('enumDescriptions') or []),
                min_items=node.get('minItems', 0), max_items=node.get('maxItems'),
                items_type=_type(items) if items else '',
            )

    walk(root['properties'], '')
    return out


def get(method: dict, path: str, missing=None):
    node = method
    for part in path.split('.'):
        if not isinstance(node, dict) or part not in node:
            return missing
        node = node[part]
    return node


def put(method: dict, path: str, value) -> None:
    parts = path.split('.')
    node = method
    for part in parts[:-1]:
        if not isinstance(node.get(part), dict):
            node[part] = {}
        node = node[part]
    node[parts[-1]] = value


def drop(method: dict, path: str) -> None:
    """Remove a key (FLASHIda then uses its default), and parents left empty."""
    parts = path.split('.')
    chain = [method]
    for part in parts[:-1]:
        nxt = chain[-1].get(part)
        if not isinstance(nxt, dict):
            return
        chain.append(nxt)
    chain[-1].pop(parts[-1], None)
    for parent, part in zip(reversed(chain[:-1]), reversed(parts[:-1])):
        if parent.get(part) == {}:
            del parent[part]


def leaves(method: dict, prefix: str = '') -> dict:
    """{dotted path: value} of a method, lists and their objects kept whole."""
    out = {}
    for name, value in method.items():
        path = prefix + name
        if isinstance(value, dict) and value:
            out.update(leaves(value, path + '.'))
        else:
            out[path] = value
    return out


def defaults() -> dict:
    """A method with the default of every key that has one. The MS1/MS2/MS3
    tolerance gets three entries, which FLASHIda requires."""
    method = {}
    for key in keys().values():
        if key.default is None or key.default == {}:
            continue
        put(method, key.path, copy.deepcopy(key.default))
    tol = get(method, 'deconvolution.tol')
    if isinstance(tol, list) and len(tol) < 3:
        put(method, 'deconvolution.tol', (tol + [tol[-1] if tol else 10] * 3)[:3])
    # The defaults of a block FLASHIda leaves out unless asked for (exploration,
    # characterization, ...) would switch it on; keep only the plain sections.
    for path in [k.path for k in keys().values() if '.exploration.' in k.path
                 or '.characterization.' in k.path or '.refine.' in k.path]:
        drop(method, path)
    return method


def _validate(method):
    """(problems, unknown keys) of a method against the bundled schema.

    A key the schema does not list is not a problem here: FLASHIda may be newer
    than the schema this FLASHApp ships, and it checks the method itself.
    """
    if not isinstance(method, dict):
        return ['A method is a JSON object ({...}).'], []
    import jsonschema

    validator = jsonschema.Draft202012Validator(schema())
    problems, unknown = [], []
    for error in sorted(validator.iter_errors(method), key=lambda e: list(e.path)):
        parts = [str(p) for p in error.path]
        if error.validator == 'additionalProperties' and isinstance(error.instance, dict):
            names = _unexpected(error.message)
            if names is not None:
                unknown += ['.'.join(parts + [n]) for n in names]
                continue
        if 'propertyNames' in error.schema_path and isinstance(error.instance, str):
            unknown.append('.'.join(parts + [error.instance]))
            continue
        problems.append(f"{'.'.join(parts) or 'method'}: {error.message}")
    tol = get(method, 'deconvolution.tol')
    if not isinstance(tol, list) or len(tol) < 3:
        if not any(p.startswith('deconvolution.tol') for p in problems):
            problems.append('deconvolution.tol: needs three tolerances (MS1, MS2, MS3).')
    return problems, sorted(set(unknown))


def _unexpected(message: str):
    """The key names of jsonschema's "Additional properties ... unexpected" message."""
    import ast
    import re

    found = re.search(r'\((.*) (?:was|were) unexpected\)\s*$', message, re.S)
    if not found:
        return None
    try:
        names = ast.literal_eval(f'({found.group(1)},)')
    except (ValueError, SyntaxError):
        return None
    return [str(n) for n in names]


def check(method) -> list:
    """Problems FLASHIda would refuse the method for ([] when it loads).

    Keys the bundled schema does not know are left out (see unknown_keys()).
    """
    return _validate(method)[0]


def unknown_keys(method) -> list:
    """Dotted paths of keys the bundled schema does not know. They are kept and
    passed to FLASHIda as they are; a newer FLASHIda may know them."""
    return _validate(method)[1]


def changes(base: dict, method: dict) -> list:
    """[(path, before, after)] where a method differs from the one it started from."""
    a, b = leaves(base), leaves(method)
    return [(p, a.get(p), b.get(p)) for p in sorted(set(a) | set(b)) if a.get(p) != b.get(p)]


# --- presets --------------------------------------------------------------------

def builtin_presets() -> dict:
    presets = {DEFAULTS_PRESET: defaults()}
    try:
        shipped = json.loads(BUILTIN_PRESETS_FILE.read_text(encoding='utf-8'))
    except (OSError, ValueError):
        shipped = {}
    for name, method in shipped.items():
        presets[name] = merge(defaults(), method)
    return presets


def merge(base: dict, over: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in over.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def user_presets(folder) -> dict:
    try:
        data = json.loads(Path(folder, USER_PRESETS_FILE).read_text(encoding='utf-8'))
    except (OSError, ValueError):
        return {}
    return {str(k): v for k, v in data.items() if isinstance(v, dict)} if isinstance(data, dict) else {}


def save_user_preset(folder, name: str, method: dict) -> None:
    presets = user_presets(folder)
    presets[name] = method
    _write_json(Path(folder, USER_PRESETS_FILE), presets)


def delete_user_preset(folder, name: str) -> None:
    presets = user_presets(folder)
    if presets.pop(name, None) is not None:
        _write_json(Path(folder, USER_PRESETS_FILE), presets)


def all_presets(folder) -> dict:
    """Built-in presets first, then the user's (a user preset never hides a
    built-in one of the same name)."""
    presets = builtin_presets()
    for name, method in user_presets(folder).items():
        presets.setdefault(name, method)
    return presets


# --- the runs table ---------------------------------------------------------------

RUNS_FILE = 'ida-runs.json'
_BAD_NAME = re.compile(r'[\\/:*?"<>|]')


def load_runs(workflow_dir) -> list:
    """[{name, preset, method}] of the runs table, in table order."""
    try:
        runs = json.loads(Path(workflow_dir, RUNS_FILE).read_text(encoding='utf-8'))
    except (OSError, ValueError):
        return []
    return [r for r in runs if isinstance(r, dict) and isinstance(r.get('method'), dict)] \
        if isinstance(runs, list) else []


def save_runs(workflow_dir, runs: list) -> None:
    _write_json(Path(workflow_dir, RUNS_FILE), runs)


def name_problem(name: str, others) -> str:
    """Why a raw file name cannot be a run ('' when it can)."""
    name = (name or '').strip()
    if not name:
        return 'Enter the raw file name.'
    if _BAD_NAME.search(name):
        return 'A raw file name cannot hold \\ / : * ? " < > |.'
    if name.lower().endswith('.raw'):
        return 'Leave out .raw.'
    if name.casefold() in {o.casefold() for o in others}:
        return 'Another run has this name.'
    return ''


def method_file_name(name: str) -> str:
    return _BAD_NAME.sub('_', name.strip()) + '.json'


def _write_json(path: Path, data) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(data, indent=2), encoding='utf-8')
    tmp.replace(path)
