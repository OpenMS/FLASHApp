"""Helpers for the FLASHDeconv OpenMS-Insight SequenceView path.

The FLASHDeconv "submitted sequence" path lets the user pick a fixed
modification on cysteine and/or methionine (``fixed_mod_cysteine`` /
``fixed_mod_methionine`` in ``src/render/sequence.py``). The legacy renderer
applied those via ``setFixedModification`` *before* computing theoretical
fragment masses, so the masses reflected the mods.

The OpenMS-Insight ``SequenceView`` computes theoretical fragment masses from
the literal sequence string (``calculate_fragment_masses_pyopenms``), and its
``compute_fixed_mods`` flag only *marks* which residue types carry a mod (for
display) -- it does NOT shift the fragment masses. To get parity we therefore
BAKE the selected fixed mods into the sequence string (e.g.
``C(Carbamidomethyl)``) so pyOpenMS includes the mass shift in every fragment.

Mapping the FLASHApp option label (e.g. ``'Carbamidomethyl (+57)'``) to an
OpenMS modification name is done by mass, mirroring ``setFixedModification``'s
``ModificationsDB().getBestModificationByDiffMonoMass`` lookup, so the baked
name is one ``AASequence.fromString`` accepts.
"""

from __future__ import annotations

from typing import Optional

# Mass shifts for the selectable fixed modifications, mirroring
# ``src/render/sequence.py`` (``fixed_mod_cysteine`` / ``fixed_mod_methionine``).
# Duplicated here (rather than imported) so this helper does not pull in
# ``src/render/sequence.py``'s top-level ``pyopenms`` import at module load: that
# keeps the helper importable/testable when pyOpenMS is absent (the mass-based
# name resolution and theoretical-mass calc degrade gracefully below).
fixed_mod_cysteine = {
    "No modification": 0,
    "Carbamidomethyl (+57)": 57.021464,
    "Carboxymethyl (+58)": 58.005479,
    "Xlink:Disulfide (-1 per C)": -1.007825,
}
fixed_mod_methionine = {
    "No modification": 0,
    "L-methionine sulfoxide (+16)": 15.994915,
    "L-methionine sulfone (+32)": 31.989829,
}


def _resolve_mod_name(diff_mass: float, residue: str) -> Optional[str]:
    """Resolve an OpenMS modification id for a mass shift on ``residue``.

    Mirrors ``setFixedModification`` (``getBestModificationByDiffMonoMass``).
    Returns None if pyOpenMS is unavailable or no modification matches.
    """
    if diff_mass == 0:
        return None
    try:
        from pyopenms import ModificationsDB
    except Exception:
        return None
    try:
        mod = ModificationsDB().getBestModificationByDiffMonoMass(
            diff_mass, 0.001, residue, 0
        )
    except Exception:
        return None
    if mod is None:
        return None
    try:
        name = mod.getId()
    except Exception:
        return None
    return name or None


def bake_fixed_modifications(
    sequence: str, fix_c: Optional[str], fix_m: Optional[str]
) -> str:
    """Return ``sequence`` with the chosen C/M fixed mods baked in as OpenMS mods.

    ``fix_c`` / ``fix_m`` are FLASHApp option labels (keys of
    ``fixed_mod_cysteine`` / ``fixed_mod_methionine``); falsy / 'No modification'
    leave that residue untouched. Unknown labels or a missing pyOpenMS leave the
    sequence unchanged (graceful degradation; the static string still renders).
    """
    if not sequence:
        return sequence

    c_name = None
    if fix_c and fix_c in fixed_mod_cysteine:
        c_name = _resolve_mod_name(fixed_mod_cysteine[fix_c], "C")
    m_name = None
    if fix_m and fix_m in fixed_mod_methionine:
        m_name = _resolve_mod_name(fixed_mod_methionine[fix_m], "M")

    if c_name is None and m_name is None:
        return sequence

    out = []
    for aa in sequence:
        out.append(aa)
        if aa == "C" and c_name is not None:
            out.append(f"({c_name})")
        elif aa == "M" and m_name is not None:
            out.append(f"({m_name})")
    return "".join(out)


def theoretical_mass(sequence: str) -> Optional[float]:
    """Monoisotopic mass of the (possibly modified) sequence, or None.

    Used to populate the SequenceView mass header (``computed_mass``). Returns
    None when pyOpenMS is unavailable so the caller simply omits the column.
    """
    if not sequence:
        return None
    try:
        from pyopenms import AASequence

        return AASequence.fromString(sequence).getMonoWeight()
    except Exception:
        return None
