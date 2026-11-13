"""Two FLASHIda methods compared, and settings of one taken over into the other.

Streamlit-free, so it is tested on its own. The Method editor's Compare mode
(src/render/ida_runs.py) draws it.

A comparison runs in three steps:

1. ``migrate(method)``: the method in the form of the current schema, and notes
   on what was changed to get there. Today it is the method as it is; a migrator
   for methods of older schema versions takes this place later, and its notes
   are shown above the comparison (and may name a row's key, ``Row.note``).
2. ``compare(a, b)``: the settings that make a difference to the acquisition.
   A key counts when FLASHIda reads it in at least one of the two methods and
   the values it works with differ. A key left out counts as its default, so a
   key left out in one and set to its default in the other is no difference.
3. ``adopt``: values of A taken over into B; a key A leaves out is left out
   of B too, and a part A does not have (identification, MS3) goes from B.
"""

import copy
from dataclasses import dataclass

from src import flashida_method as fm
from src.render import ida_sections as sections

UNKNOWN = 'Unknown keys'  # the tab of keys the bundled schema does not know
# Neither the host nor the engine acts on these (the schema says so): no difference to an acquisition
NO_EFFECT = {'global.method_name', 'global.method_description'}
_MISSING = object()


def migrate(method: dict) -> tuple:
    """(method in the form of the current schema, notes). The place of the
    method migrator: until there is one, the method is taken as it is."""
    return copy.deepcopy(method), []


@dataclass
class Row:
    """One setting in which two methods differ."""
    path: str
    tab: str
    part: str
    a: object  # the value as the method holds it (None: left out)
    b: object
    note: str = ''  # what the migrator did to this key, if anything


def _tab_of(part: str) -> str:
    for name, (_, parts) in sections.TAB_PARTS.items():
        if part in {p for p, _ in parts}:
            return name
    return 'Global'


def effective(method: dict, path: str):
    """The value FLASHIda works with: the method's, else the one the key takes
    from MS1, else FLASHIda's default."""
    value = fm.get(method, path, missing=_MISSING)
    if value is not _MISSING and value not in (None, ''):
        return value
    key = fm.keys().get(path)
    if key is not None and (key.default_text or '').startswith('inherits ms1'):
        ms1 = 'ms_settings.ms1.' + path.rsplit('.', 1)[-1]
        inherited = fm.get(method, ms1, missing=_MISSING)
        if inherited is not _MISSING and inherited not in (None, ''):
            return inherited
        return sections.default_value(ms1)
    return sections.default_value(path)


# The keys that switch a part on and off are no field of it, but belong to it
SWITCHES = {'quantification.enabled': 'quant'}


def part_of(path: str) -> str:
    return SWITCHES.get(path) or sections.owner(path)


def is_read(method: dict, path: str) -> bool:
    """Whether FLASHIda reads the key with the method as it is: its part is on
    (a switch is always read) and nothing else switches it off."""
    if path in SWITCHES:
        return True
    part = part_of(path)
    on = part == 'global' or part in sections.present(method)
    return on and sections.inactive(method, path) is None


# Parts that are on while their block is there; no key says so, so a row of its own
BLOCK_SWITCHES = {'identification': ('ident', 'Identify in real time')}


def switch_label(path: str):
    """The name of a row that switches a block, or None."""
    return BLOCK_SWITCHES[path][1] if path in BLOCK_SWITCHES else None


def compare(a: dict, b: dict) -> list:
    """[Row] of the settings in which ``b`` would acquire differently from
    ``a``, in the order of the editor's tabs, then the keys the schema does not
    know (their effect cannot be judged, so any difference counts)."""
    order = {name: i for i, name in enumerate(sections.TAB_PARTS)}
    rows = []
    for block, (part, _) in BLOCK_SWITCHES.items():
        on_a, on_b = part in sections.present(a), part in sections.present(b)
        if on_a != on_b:
            rows.append(Row(block, _tab_of(part), part, 'on' if on_a else 'off', 'on' if on_b else 'off'))
    for path in fm.keys():
        if path in NO_EFFECT or not (is_read(a, path) or is_read(b, path)):
            continue
        if effective(a, path) == effective(b, path):
            continue
        part = part_of(path)
        rows.append(Row(path, _tab_of(part), part,
                        fm.get(a, path), fm.get(b, path)))
    rows.sort(key=lambda r: order.get(r.tab, len(order)))
    # An unknown key may be a whole section: compare each of its values
    la, lb = fm.leaves(a), fm.leaves(b)
    unknown = set(fm.unknown_keys(a)) | set(fm.unknown_keys(b))
    for path in sorted(p for p in set(la) | set(lb) if any(p == u or p.startswith(u + '.') for u in unknown)):
        if la.get(path) != lb.get(path):
            rows.append(Row(path, UNKNOWN, '', la.get(path), lb.get(path)))
    return rows


def all_settings(a: dict, b: dict) -> list:
    """(tab, path, value A, value B, differs) of every key, with the values
    FLASHIda works with (a key left out shows its default), in tab order."""
    differs = {r.path for r in compare(a, b)}
    order = {name: i for i, name in enumerate(sections.TAB_PARTS)}
    out = [(_tab_of(part_of(p)), p, effective(a, p), effective(b, p), p in differs) for p in fm.keys()]
    return sorted(out, key=lambda r: order.get(r[0], len(order)))


def adopt(b: dict, a: dict, paths) -> dict:
    """``b`` with the values of ``a`` at ``paths``; a key ``a`` leaves out is
    removed from ``b``. ``b`` itself is left as it was."""
    out = copy.deepcopy(b)
    for path in paths:
        value = fm.get(a, path, missing=_MISSING)
        if path in BLOCK_SWITCHES:  # the whole block, or none
            if value is _MISSING:
                fm.drop(out, path)
            else:
                fm.put(out, path, copy.deepcopy(value))
        elif value is _MISSING and path in sections.CLEARS:  # MS3 off: its block goes
            fm.drop(out, sections.CLEARS[path])
        elif value is _MISSING:
            fm.drop(out, path)
        else:
            fm.put(out, path, copy.deepcopy(value))
    return out
