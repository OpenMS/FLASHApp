"""FLASHAnnotator output -> per-proteoform fragment annotations for the viewer.

FLASHAnnotator is run on the FLASHTnT proteoform file (-out_pro), so its
`-out` file is that file with the raw annotation columns appended, and the
first column of its `-out_annotation` file ("PrSMIndex") carries the
ProteoformIndex. Every accepted fragment isotope envelope becomes one row of
`annotated_fragments` in the proteoform's sequence-view entry.

Positions are converted to 0-based residue indices of the full database
sequence the Sequence View draws (`start`..`end`, inclusive), so the frontend
does not repeat FLASHAnnotator's coordinate rules.
"""
import numpy as np
import pandas as pd

PREFIX_IONS = {'a', 'b', 'c'}
SUFFIX_IONS = {'x', 'y', 'z', 'zp1', 'zp2'}
# Display names of FLASHTnT/FLASHAnnotator ion types (zp1: z+1, zp2: z+2)
ION_LABELS = {'zp1': 'z+1', 'zp2': 'z+2'}

ANNOTATION_COLUMNS = [
    'PrSMIndex', 'Scan', 'IonType', 'Length', 'Bond', 'Charge', 'MonoMass',
    'MostAbundantMz', 'MatchedPeaks', 'PredictedPeaks', 'Cosine', 'Pearson',
    'Intensity', 'PPMError', 'IntensityRank', 'LocalDensity', 'DecoyShift',
    'InternalFirst',
]

# Proteoform-level columns FLASHAnnotator appends to its -out file
SUMMARY_COLUMNS = {
    'FAScore': 'fa_score',
    'RawScore': 'raw_score',
    'RawMatchedFragments': 'matched_fragments',
    'RawFragmentHypotheses': 'fragment_hypotheses',
    'RawCoveredBonds': 'covered_bonds',
    'RawCoverage(%)': 'coverage',
    'RawCalibrationPPM': 'calibration_ppm',
    'RawInternalScore': 'internal_score',
    'RawInternalMatched': 'internal_matched',
    'RawInternalHypotheses': 'internal_hypotheses',
    'RawAuthenticCleavages': 'authentic_cleavages',
    'RawCoverageWithInternal(%)': 'coverage_with_internal',
}


def read_annotation(annotation_tsv):
    """Read the -out_annotation file; an absent or empty file gives an empty
    frame with the expected columns."""
    if annotation_tsv is None:
        return pd.DataFrame(columns=ANNOTATION_COLUMNS)
    try:
        df = pd.read_csv(annotation_tsv, sep='\t')
    except (pd.errors.EmptyDataError, FileNotFoundError):
        return pd.DataFrame(columns=ANNOTATION_COLUMNS)
    # Decoy envelopes are only written with report_decoy_envelopes
    if 'DecoyShift' in df.columns:
        df = df[df['DecoyShift'] == 0]
    return df


def read_summary(annotated_tsv):
    """{ProteoformIndex: {summary key: value}} from FLASHAnnotator's -out."""
    if annotated_tsv is None:
        return {}
    try:
        df = pd.read_csv(annotated_tsv, sep='\t')
    except (pd.errors.EmptyDataError, FileNotFoundError):
        return {}
    if 'ProteoformIndex' not in df.columns:
        return {}
    # FLASHTnT run with -in_raw already writes RawScore, RawCoverage(%), ...;
    # FLASHAnnotator appends its own after them, which pandas reads as "<name>.1"
    present = {
        c: (f'{c}.1' if f'{c}.1' in df.columns else c)
        for c in SUMMARY_COLUMNS if c in df.columns
    }
    summary = {}
    for _, row in df.iterrows():
        values = {}
        for c, column in present.items():
            v = row[column]
            # -1 marks values FLASHAnnotator could not compute, NaN the signed ones
            values[SUMMARY_COLUMNS[c]] = None if (pd.isna(v) or v == -1) else float(v)
        summary[int(row['ProteoformIndex'])] = values
    return summary


def proteoform_offset(database_sequence, proteoform_sequence, start, end):
    """0-based index of the proteoform's first residue in the database
    sequence, following FLASHAnnotator: StartPosition when determined, else
    EndPosition minus the proteoform length, else where the proteoform occurs."""
    if start > 0:
        return start - 1
    if end > 0:
        return max(0, end - len(proteoform_sequence))
    at = str(database_sequence).find(str(proteoform_sequence))
    return at if at >= 0 else 0


def _fragment_name(ion, kind, length, first, last):
    if kind == 'I':
        return f'{ion}[{first + 1}-{last + 1}]'
    if ion in ION_LABELS:
        return f'{ION_LABELS[ion]}({length})'
    return f'{ion}{length}'


def build_fragments(annotation_df, offset, proteoform_length):
    """Rows of one proteoform's envelopes -> list of fragment dicts.

    `first`/`last` are 0-based residues of the proteoform, `start`/`end` the
    same residues in database-sequence coordinates. `kind` is 'N' (a/b/c),
    'C' (x/y/z) or 'I' (internal, ion type 'N-end type/C-end type')."""
    fragments = []
    for row in annotation_df.itertuples(index=False):
        ion = str(row.IonType)
        length = int(row.Length)
        internal_first = int(row.InternalFirst) if not pd.isna(row.InternalFirst) else 0
        if ion.startswith('int(') or internal_first > 0:
            kind = 'I'
            ion = ion[4:-1] if ion.startswith('int(') else ion
            first = internal_first - 1
            last = first + length - 1
        elif ion in PREFIX_IONS:
            kind = 'N'
            first, last = 0, length - 1
        elif ion in SUFFIX_IONS:
            kind = 'C'
            first, last = proteoform_length - length, proteoform_length - 1
        else:
            continue
        ppm = float(row.PPMError)
        mono_mass = float(row.MonoMass)
        fragments.append({
            'name': _fragment_name(ion, kind, length, first, last),
            'ion_type': ION_LABELS.get(ion, ion),
            'kind': kind,
            'length': length,
            'first': first,
            'last': last,
            'start': first + offset,
            'end': last + offset,
            'charge': int(row.Charge),
            'mono_mass': mono_mass,
            'observed_mass': mono_mass * (1 + ppm * 1e-6),
            'mz': float(row.MostAbundantMz),
            'ppm_error': ppm,
            'matched_peaks': int(row.MatchedPeaks),
            'predicted_peaks': int(row.PredictedPeaks),
            'cosine': float(row.Cosine),
            'intensity': float(row.Intensity),
        })
    fragments.sort(key=lambda f: ({'N': 0, 'C': 1, 'I': 2}[f['kind']], f['start'], f['end'], f['charge']))
    return fragments


def group_annotation(annotation_df):
    """{ProteoformIndex: frame of its envelopes}"""
    if annotation_df.empty:
        return {}
    return {
        int(pid): g for pid, g in annotation_df.groupby(annotation_df['PrSMIndex'].astype(np.int64))
    }
