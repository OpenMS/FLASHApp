import pandas as pd


def _split_ints(s):
    return [int(x) for x in str(s).replace(';', ' ').split()
            if x.strip().lstrip('-').isdigit()]


def build_tagspace_to_proteoform_map(raw_tag_df, raw_protein_df):
    """Map FLASHTagger tag-space ProteoformIndex -> protein-space index.

    tags.tsv enumerates all candidate proteoforms (tag-space, incl. decoys);
    protein.tsv enumerates surviving proteoforms (protein-space). The bridge is
    protein.tsv.TagIndices crossed with tags.tsv.ProteoformIndex; the relation
    is a strictly monotonic bijection, resolved by greedy monotonic assignment
    over proteoforms in ascending protein-space order. Returns
    {tag_space_index: protein_space_index}; tag-space indices with no surviving
    proteoform are absent (callers map them to -1).
    """
    ti_to_qset = {
        int(ti): set(_split_ints(pis))
        for ti, pis in zip(raw_tag_df['TagIndex'], raw_tag_df['ProteoformIndex'])
    }
    q_to_p = {}
    prev_q = -1
    ordered = raw_protein_df.sort_values('ProteoformIndex')
    for p, tagidx in zip(ordered['ProteoformIndex'].astype(int), ordered['TagIndices']):
        cand = None
        for t in _split_ints(tagidx):
            s = ti_to_qset.get(t, set())
            cand = s if cand is None else (cand & s)
        if not cand:                          # empty intersection -> union fallback
            cand = set()
            for t in _split_ints(tagidx):
                cand |= ti_to_qset.get(t, set())
        nxt = [q for q in sorted(cand) if q > prev_q]
        if not nxt:
            continue
        q_to_p[nxt[0]] = int(p)
        prev_q = nxt[0]
    return q_to_p
