"""Structure lookup for the FLASHTnT Structure View (src/render/structure.py)."""
import json

import pytest
import requests

from src.render import structure
from src.render.layout_editor.catalog import TNT_COMPONENTS, validate_layout
from src.render.structure import (
    StructureStore, get_structure, structure_format, structure_key, uniprot_accession,
)

PDB = (
    b"ATOM      1  N   MET A   1      11.104   6.134  -6.504  1.00 91.00           N\n"
    b"ATOM      2  CA  MET A   1      11.639   6.071  -5.147  1.00 91.00           C\n"
    b"END\n"
)
CIF = b"data_AF-P60844-F1\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\nATOM 1\n"


@pytest.mark.parametrize('accession, expected', [
    ('sp|P60844|AQPZ_ECOLI', 'P60844'),
    ('tr|A0A023GPI8|A0A023GPI8_CANAL', 'A0A023GPI8'),
    ('P0DOX5', 'P0DOX5'),
    ('P0DOX5-2', 'P0DOX5'),
    ('DECOY_sp|P60844|AQPZ_ECOLI', None),
    ('AQPZ', None),
    ('SigmaMabH', None),
    (None, None),
])
def test_uniprot_accession(accession, expected):
    assert uniprot_accession(accession) == expected


def test_structure_key_falls_back_to_a_safe_file_name():
    assert structure_key('sp|P60844|AQPZ_ECOLI') == 'P60844'
    assert structure_key('NIST mAb/Heavy Chain') == 'NIST_mAb_Heavy_Chain'
    assert structure_key('../..') is None


def test_structure_format():
    assert structure_format('model.pdb', PDB) == 'pdb'
    assert structure_format('model.cif', CIF) == 'cif'
    assert structure_format('model', CIF) == 'cif'
    assert structure_format('model.pdb', b'<html>not found</html>') is None
    assert structure_format('model.cif', PDB) is None


def test_store_roundtrip(tmp_path):
    store = StructureStore(tmp_path / 'structures')
    assert store.find('P60844') is None
    store.save('P60844', PDB, 'pdb', 'upload', 'aqpz.pdb')
    path, meta = store.find('P60844')
    assert path.read_bytes() == PDB
    assert meta['source'] == 'upload'
    # Replacing a structure removes the old file
    store.save('P60844', CIF, 'cif', 'upload', 'aqpz.cif')
    assert sorted(p.name for p in (tmp_path / 'structures').iterdir()) == ['P60844.cif', 'P60844.json']
    assert list(store.list()) == ['P60844']
    store.remove('P60844')
    assert store.list() == {}


class FakeResponse:
    def __init__(self, status, body):
        self.status_code = status
        self.body = body

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(response=self)

    def iter_content(self, size):
        yield self.body


class FakeSession:
    def __init__(self, routes):
        self.routes = routes
        self.calls = []

    def get(self, url, **kwargs):
        self.calls.append(url)
        if url not in self.routes:
            raise requests.ConnectionError('offline')
        return FakeResponse(*self.routes[url])


API = 'https://alphafold.ebi.ac.uk/api/prediction/P60844'
FILE = 'https://alphafold.ebi.ac.uk/files/AF-P60844-F1-model_v4.pdb'


@pytest.fixture(autouse=True)
def clear_failures():
    structure._failures.clear()


def test_downloads_alphafold_model_once(tmp_path):
    session = FakeSession({
        API: (200, json.dumps([{'entryId': 'AF-P60844-F1', 'pdbUrl': FILE}]).encode()),
        FILE: (200, PDB),
    })
    store = StructureStore(tmp_path)
    result = get_structure(store, 'sp|P60844|AQPZ_ECOLI', session=session)
    assert result['status'] == 'ok'
    assert result['source'] == 'alphafold'
    assert result['format'] == 'pdb'
    assert result['name'] == 'AF-P60844-F1'
    assert result['data'] == PDB.decode()
    assert result['url'] == 'https://alphafold.ebi.ac.uk/entry/P60844'
    # Stored: no second download
    assert get_structure(store, 'sp|P60844|AQPZ_ECOLI', session=session)['status'] == 'ok'
    assert len(session.calls) == 2


def test_upload_wins_and_needs_no_network(tmp_path):
    store = StructureStore(tmp_path)
    store.save('P60844', CIF, 'cif', 'upload', 'mine.cif')
    session = FakeSession({})
    result = get_structure(store, 'sp|P60844|AQPZ_ECOLI', session=session)
    assert (result['status'], result['source'], result['format']) == ('ok', 'upload', 'cif')
    assert session.calls == []


def test_offline_failure_is_not_retried_on_every_rerun(tmp_path):
    session = FakeSession({})
    store = StructureStore(tmp_path)
    result = get_structure(store, 'P60844', session=session, now=1000)
    assert result['status'] == 'unavailable'
    assert 'no connection' in result['message']
    get_structure(store, 'P60844', session=session, now=1100)
    assert len(session.calls) == 1
    get_structure(store, 'P60844', session=session, now=1000 + structure.RETRY_AFTER + 1)
    assert len(session.calls) == 2


def test_missing_model(tmp_path):
    session = FakeSession({API: (404, b'')})
    result = get_structure(StructureStore(tmp_path), 'P60844', session=session)
    assert result['status'] == 'unavailable'
    assert 'no model for this accession' in result['message']


def test_refuses_downloads_outside_alphafold(tmp_path):
    session = FakeSession({
        API: (200, json.dumps([{'pdbUrl': 'http://169.254.169.254/latest'}]).encode()),
    })
    result = get_structure(StructureStore(tmp_path), 'P60844', session=session)
    assert result['status'] == 'unavailable'
    assert session.calls == [API]


def test_no_download_without_uniprot_accession_or_when_disabled(tmp_path):
    session = FakeSession({})
    store = StructureStore(tmp_path)
    assert get_structure(store, 'SigmaMabH', session=session)['status'] == 'missing'
    assert get_structure(store, 'P60844', download=False, session=session)['status'] == 'missing'
    assert session.calls == []


def test_structure_view_needs_protein_table():
    assert validate_layout([[['protein_table'], ['sequence_view', 'structure_view']]], TNT_COMPONENTS) == ''
    assert 'needs Protein table' in validate_layout([[['structure_view']]], TNT_COMPONENTS)
