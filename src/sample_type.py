"""Native or lysate: which FLASHTnT settings suit a set of MS runs.

The two presets in presets.json ("Native", "Lysate") differ in how FLASHTnT
searches (-native) and in the decoy database size. Native or targeted
top-down runs hold one or a few MS2 spectra and often no MS1 survey; an LC-MS/MS
run of a cell lysate has MS1 surveys and many MS2 spectra. That is all the
recommendation looks at, so it reads only the spectra's MS levels (the
`ms level` cvParam), stops as soon as a file is clearly a lysate run, and never
loads peaks. Kept free of Streamlit.
"""

import re
from dataclasses import dataclass
from pathlib import Path

NATIVE = 'Native'
LYSATE = 'Lysate'

# A run with more MS2 spectra than this and MS1 surveys is an LC-MS/MS run
NATIVE_MAX_MS2 = 10

_MS_LEVEL = re.compile(rb'<cvParam[^>]*?MS:1000511[^>]*?value="(\d+)"')
_CHUNK = 1 << 20


@dataclass
class RunSpectra:
    ms1: int
    ms2: int
    complete: bool = True  # False when counting stopped early

    @property
    def sample_type(self):
        return classify(self.ms1, self.ms2)


def classify(ms1, ms2):
    """NATIVE for runs without MS1 surveys or with at most NATIVE_MAX_MS2
    MS2 spectra, LYSATE otherwise."""
    if ms1 == 0 or ms2 <= NATIVE_MAX_MS2:
        return NATIVE
    return LYSATE


def count_spectra(mzml_path):
    """MS1 and MS2 spectra of an mzML file, counted from the `ms level`
    cvParams. Stops once the file is clearly a lysate run (`complete` is
    then False)."""
    ms1 = ms2 = 0
    tail = b''
    with open(mzml_path, 'rb') as f:
        while True:
            chunk = f.read(_CHUNK)
            if not chunk:
                break
            data = tail + chunk
            # Keep an unfinished tag for the next chunk
            cut = data.rfind(b'<')
            if cut < 0:
                cut = len(data)
            for match in _MS_LEVEL.finditer(data, 0, cut):
                if match.group(1) == b'1':
                    ms1 += 1
                else:
                    ms2 += 1
            tail = data[cut:]
            if ms1 > 0 and ms2 > NATIVE_MAX_MS2:
                return RunSpectra(ms1, ms2, complete=False)
    for match in _MS_LEVEL.finditer(tail):
        if match.group(1) == b'1':
            ms1 += 1
        else:
            ms2 += 1
    return RunSpectra(ms1, ms2)


@dataclass
class Recommendation:
    sample_type: str | None  # NATIVE, LYSATE, or None (mixed or unknown)
    runs: dict  # path -> RunSpectra
    unreadable: list  # paths that are not (yet) mzML

    def reason(self):
        """One sentence on what the recommendation is based on."""
        native = [p for p, r in self.runs.items() if r.sample_type == NATIVE]
        lysate = [p for p, r in self.runs.items() if r.sample_type == LYSATE]
        if self.sample_type == NATIVE:
            ms2 = sum(r.ms2 for r in self.runs.values())
            no_ms1 = all(r.ms1 == 0 for r in self.runs.values())
            survey = 'no MS1 survey scans' if no_ms1 else f'at most {NATIVE_MAX_MS2} MS2 spectra per file'
            subject = 'The selected file holds' if len(native) == 1 else f'The {len(native)} selected files hold'
            return f'{subject} {ms2} MS2 {"spectrum" if ms2 == 1 else "spectra"} with {survey}.'
        if self.sample_type == LYSATE:
            subject = 'The selected file is an LC-MS/MS run' if len(lysate) == 1 else f'The {len(lysate)} selected files are LC-MS/MS runs'
            return (f'{subject}: MS1 survey scans and more than {NATIVE_MAX_MS2} '
                    'MS2 spectra per file.')
        if native and lysate:
            return (f'{_files(len(native))} look native and {_files(len(lysate))} '
                    'look like lysate runs; search them separately.')
        return ''


def _files(n):
    return f'{n} file' if n == 1 else f'{n} files'


def mzml_for(path, converted_dir=None):
    """The mzML to inspect for an input: the file itself, or for Thermo .raw /
    Bruker .d the converted mzML when it already exists."""
    path = Path(path)
    if path.suffix.lower() == '.mzml':
        return path
    if converted_dir is not None:
        candidate = Path(converted_dir, path.stem + '.mzML')
        if candidate.exists():
            return candidate
    return None


def recommend(paths, converted_dir=None, counter=count_spectra):
    """Recommendation for the selected input files."""
    runs, unreadable = {}, []
    for path in paths or []:
        mzml = mzml_for(path, converted_dir)
        if mzml is None or not mzml.exists():
            unreadable.append(str(path))
            continue
        runs[str(path)] = counter(mzml)
    types = {r.sample_type for r in runs.values()}
    sample_type = types.pop() if len(types) == 1 else None
    return Recommendation(sample_type, runs, unreadable)


def of_result(file_manager, dataset_id):
    """Sample type of a FLASHTnT result: NATIVE when it was searched with
    -native, else from the MS levels of its spectra (classify)."""
    if file_manager.result_exists(dataset_id, 'FTnT_parameters'):
        settings = file_manager.get_results(dataset_id, ['FTnT_parameters'])['FTnT_parameters']
        # Only settings that differ from FLASHTnT's defaults are recorded
        if isinstance(settings, dict) and settings.get('native'):
            return NATIVE
    path = file_manager.get_result_path(dataset_id, 'scan_table')
    if path is None:
        return None
    import polars as pl  # pylint: disable=import-outside-toplevel
    levels = pl.scan_parquet(path).select(
        (pl.col('MSLevel') == 1).sum().alias('ms1'), (pl.col('MSLevel') > 1).sum().alias('ms2')
    ).collect().row(0, named=True)
    return classify(levels['ms1'], levels['ms2'])
