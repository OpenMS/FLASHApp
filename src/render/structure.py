"""3D structures for the FLASHTnT Structure View.

A structure is looked up by the protein's UniProt accession. Files live in one
folder per workspace (``<workspace>/structures/``), so every dataset of the
workspace shares them and a structure that was downloaded once stays available
offline. A file there is either a model fetched from the AlphaFold Protein
Structure Database or a PDB/mmCIF file the user uploaded; an upload wins.

Streamlit-free, so it can be tested and reused outside the viewer.
"""
from __future__ import annotations

import json
import re
import time
from pathlib import Path
from urllib.parse import urlsplit

import requests

ALPHAFOLD_API = 'https://alphafold.ebi.ac.uk/api/prediction/{accession}'
ALPHAFOLD_ENTRY = 'https://alphafold.ebi.ac.uk/entry/{accession}'
ALPHAFOLD_HOST = 'alphafold.ebi.ac.uk'
USER_AGENT = 'FLASHApp-StructureView/1.0'
TIMEOUT = 15
MAX_BYTES = 50 * 1024 * 1024
# A failed download is retried after this many seconds, not on every rerun
RETRY_AFTER = 600

FORMATS = {'.pdb': 'pdb', '.ent': 'pdb', '.cif': 'cif', '.mmcif': 'cif'}

# https://www.uniprot.org/help/accession_numbers
_UNIPROT = re.compile(
    r'^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?$'
)

_failures: dict[str, tuple[float, str]] = {}


def uniprot_accession(accession) -> str | None:
    """UniProt accession of a FASTA identifier (``sp|P60844|AQPZ_ECOLI``,
    ``P60844``, ``P60844-2``), without isoform suffix; None if there is none."""
    if not isinstance(accession, str):
        return None
    accession = accession.strip()
    if accession.startswith('DECOY_'):
        return None
    parts = accession.split('|')
    candidates = [parts[1]] if len(parts) >= 3 and parts[0] in ('sp', 'tr') else parts
    for candidate in candidates:
        candidate = candidate.strip().split()[0] if candidate.strip() else ''
        if _UNIPROT.match(candidate):
            return candidate.split('-')[0]
    return None


def structure_key(accession) -> str | None:
    """File name stem a structure for this protein is stored under: the UniProt
    accession when there is one, else the identifier made file-name safe."""
    uniprot = uniprot_accession(accession)
    if uniprot is not None:
        return uniprot
    if not isinstance(accession, str) or not accession.strip():
        return None
    key = re.sub(r'[^A-Za-z0-9_.-]+', '_', accession.strip()).strip('._')
    return key[:100] or None


def structure_format(file_name: str, data: bytes) -> str | None:
    """'pdb' or 'cif' from the file name, else from the content; None if it is neither."""
    fmt = FORMATS.get(Path(file_name).suffix.lower())
    head = data[:4096].decode('utf-8', errors='replace')
    looks_cif = head.lstrip().startswith('data_') or '_atom_site.' in data[:200000].decode('utf-8', errors='replace')
    looks_pdb = re.search(r'^(ATOM  |HETATM|HEADER|MODEL |CRYST1|REMARK)', head, re.MULTILINE) is not None
    if fmt == 'cif' and looks_cif:
        return 'cif'
    if fmt == 'pdb' and looks_pdb:
        return 'pdb'
    if fmt is None:
        return 'cif' if looks_cif else ('pdb' if looks_pdb else None)
    return None


class StructureStore:
    """Structure files of one workspace, one per protein key."""

    def __init__(self, root):
        self.root = Path(root)

    def _meta_path(self, key):
        return Path(self.root, f'{key}.json')

    def find(self, key):
        """(path, metadata) of the stored structure, or None."""
        meta_path = self._meta_path(key)
        if not meta_path.exists():
            return None
        try:
            meta = json.loads(meta_path.read_text())
        except (OSError, ValueError):
            return None
        path = Path(self.root, meta.get('file', ''))
        if not meta.get('file') or not path.is_file():
            return None
        return path, meta

    def save(self, key, data: bytes, fmt: str, source: str, name: str, url: str | None = None):
        self.root.mkdir(parents=True, exist_ok=True)
        self.remove(key)
        file_name = f'{key}.{fmt}'
        Path(self.root, file_name).write_bytes(data)
        meta = {'file': file_name, 'format': fmt, 'source': source, 'name': name, 'url': url}
        self._meta_path(key).write_text(json.dumps(meta))
        return Path(self.root, file_name)

    def remove(self, key):
        found = self.find(key)
        if found is not None:
            found[0].unlink(missing_ok=True)
        self._meta_path(key).unlink(missing_ok=True)

    def list(self):
        """{key: metadata} of every stored structure."""
        out = {}
        if not self.root.is_dir():
            return out
        for meta_path in sorted(self.root.glob('*.json')):
            found = self.find(meta_path.stem)
            if found is not None:
                out[meta_path.stem] = found[1]
        return out


def _get(session, url, max_bytes=MAX_BYTES):
    parts = urlsplit(url)
    # Only follow links into AlphaFold DB itself
    if parts.scheme != 'https' or parts.hostname != ALPHAFOLD_HOST:
        raise ValueError(f'unexpected download location {url}')
    response = session.get(url, timeout=TIMEOUT, stream=True, allow_redirects=False,
                           headers={'User-Agent': USER_AGENT})
    response.raise_for_status()
    data = bytearray()
    for chunk in response.iter_content(65536):
        data.extend(chunk)
        if len(data) > max_bytes:
            raise ValueError('structure file too large')
    return bytes(data)


def fetch_alphafold(accession: str, session=None):
    """(data, format, model name, entry url) of the AlphaFold DB model."""
    session = session or requests.Session()
    entries = json.loads(_get(session, ALPHAFOLD_API.format(accession=accession), 1024 * 1024))
    if not entries:
        raise LookupError('no AlphaFold DB model')
    entry = entries[0]
    url = entry.get('pdbUrl') or entry.get('cifUrl')
    if not url:
        raise LookupError('no AlphaFold DB model file')
    data = _get(session, url)
    fmt = structure_format(url, data)
    if fmt is None:
        raise ValueError('AlphaFold DB returned no structure file')
    name = entry.get('entryId') or Path(urlsplit(url).path).stem
    return data, fmt, name, ALPHAFOLD_ENTRY.format(accession=accession)


def get_structure(store: StructureStore, accession, download=True, session=None, now=None):
    """Structure for the Structure View, as sent to the frontend.

    ``status`` is 'ok' (with ``format``, ``data``, ``source``, ``name``,
    ``url``), 'missing' (nothing stored and nothing to download) or
    'unavailable' (the download failed; ``message`` says why)."""
    key = structure_key(accession)
    result = {'accession': str(accession)}
    if key is None:
        return {**result, 'status': 'missing',
                'message': 'This protein has no accession to look up a structure for.'}
    # An upload wins; the AlphaFold model stays stored next to it
    found = store.find(upload_key(key)) or store.find(key)
    if found is None and download and uniprot_accession(accession) is not None:
        now = time.time() if now is None else now
        failed = _failures.get(key)
        if failed is None or now - failed[0] > RETRY_AFTER:
            try:
                data, fmt, name, url = fetch_alphafold(key, session)
                store.save(key, data, fmt, 'alphafold', name, url)
                _failures.pop(key, None)
                found = store.find(key)
            except (requests.RequestException, LookupError, ValueError, OSError) as error:
                _failures[key] = (now, _describe(error))
        if found is None:
            return {**result, 'status': 'unavailable', 'message': (
                f'No structure for {key} could be downloaded from AlphaFold DB '
                f'({_failures[key][1]}). Upload a PDB or mmCIF file for it here.'
            )}
    if found is None:
        hint = '' if uniprot_accession(accession) is not None or not download else (
            ' Its identifier is not a UniProt accession, so none can be downloaded.')
        return {**result, 'status': 'missing', 'message': (
            f'No structure stored for {key}.{hint} Upload a PDB or mmCIF file for it here.'
        )}
    path, meta = found
    return {
        **result, 'status': 'ok', 'format': meta.get('format', 'pdb'),
        'data': path.read_text(errors='replace'), 'source': meta.get('source', 'upload'),
        'name': meta.get('name'), 'url': meta.get('url'),
    }


def upload_key(key):
    """Store key of the structure uploaded for a protein."""
    return f'{key}.upload'


def apply_structure_action(store: StructureStore, action) -> tuple[bool, str | None]:
    """Apply an upload or removal sent by the Structure View.

    ``action`` is ``{'op': 'upload', 'accession', 'name', 'data'}`` (the file's
    text) or ``{'op': 'remove', 'accession'}``. Returns (changed, error)."""
    if not isinstance(action, dict):
        return False, None
    key = structure_key(action.get('accession'))
    if key is None:
        return False, 'This protein has no identifier to store a structure under.'
    if action.get('op') == 'remove':
        if store.find(upload_key(key)) is None:
            return False, None
        store.remove(upload_key(key))
        return True, None
    if action.get('op') != 'upload':
        return False, None
    name = str(action.get('name') or 'structure')
    text = action.get('data')
    data = text.encode() if isinstance(text, str) else b''
    if len(data) > MAX_BYTES:
        return False, f'{name} is larger than {MAX_BYTES // (1024 * 1024)} MB.'
    fmt = structure_format(name, data)
    if fmt is None:
        return False, f'{name} is not a PDB or mmCIF file.'
    store.save(upload_key(key), data, fmt, 'upload', Path(name).name)
    return True, None


def _describe(error):
    if isinstance(error, requests.HTTPError) and error.response is not None:
        if error.response.status_code == 404:
            return 'no model for this accession'
        return f'server answered {error.response.status_code}'
    if isinstance(error, (requests.ConnectionError, requests.Timeout)):
        return 'no connection'
    return str(error) or type(error).__name__


def structure_dir(file_manager):
    """Folder of the workspace's structure files, next to its result cache."""
    return Path(file_manager.workflow_dir, 'structures')


def download_enabled(settings) -> bool:
    """AlphaFold DB downloads are on unless settings.json turns them off."""
    return bool((settings or {}).get('structure_view', {}).get('alphafold_download', True))
