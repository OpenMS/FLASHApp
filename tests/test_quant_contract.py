"""The comparison of runs: FLASHQuantConsensus output read into the shared
contract (src/quant/contract.py), stored and analysed as a comparison dataset
(src/quant/session.py), and the design of a workflow's runs (src/quant/run.py)."""

from pathlib import Path

import polars as pl
import pytest

pytest.importorskip('pyopenms')

from src.quant import contract, session  # noqa: E402
from src.quant.run import design_of, run_name  # noqa: E402
from src.workflow import ExperimentalDesign  # noqa: E402
from src.workflow.FileManager import FileManager  # noqa: E402

DATA = Path(__file__).parent / 'data' / 'quant'
CONSENSUS = DATA / 'consensus.tsv'
IDENTIFICATIONS = DATA / 'identifications.tsv'


def test_consensus_alone_gives_the_feature_level():
    levels, runs = contract.build(CONSENSUS)
    assert list(levels) == ['feature']
    assert runs == ['FLASHQuantConsensus_input1', 'FLASHQuantConsensus_input2']
    long, annotations = levels['feature']
    assert set(long.columns) >= {'id', 'sample', 'intensity', 'status'}
    # One row per feature and run
    assert long.height == annotations.height * len(runs)
    assert set(long['status'].unique()) <= set(contract.STATUS_ORDER)


def test_identifications_add_proteoforms_and_proteins():
    levels, _ = contract.build(CONSENSUS, IDENTIFICATIONS)
    assert list(levels) == ['feature', 'proteoform', 'protein']
    features = levels['feature'][1]
    # Feature 3 has PrSMs of two proteoforms: it is ambiguous
    ambiguous = features.filter(pl.col('ambiguous'))['id'].to_list()
    assert ambiguous == ['3']
    for level in ('proteoform', 'protein'):
        long, annotations = levels[level]
        assert set(long['id']) == set(annotations['id'])
        # A rolled-up value is the sum of its parts
        assert long['intensity'].drop_nulls().min() > 0


def test_missing_values_stay_missing():
    _, long, _ = contract.read_consensus(CONSENSUS)
    assert long.filter(pl.col('intensity') <= 0).height == 0
    missing = long.filter(pl.col('status') == 'missing')
    assert missing['intensity'].null_count() == missing.height


def test_rollup_status_is_the_best_status_of_parts_with_a_value():
    long = pl.DataFrame({
        'id': ['a', 'b', 'a', 'b'], 'sample': ['r1', 'r1', 'r2', 'r2'],
        'intensity': [1.0, 2.0, None, 5.0],
        'status': ['detected', 'transferred', 'missing', 'transferred'],
        'score': [None] * 4, 'qvalue': [None] * 4,
    }, schema_overrides={'score': pl.Float64, 'qvalue': pl.Float64})
    mapping = pl.DataFrame({'id': ['a', 'b'], 'parent_id': ['P', 'P']})
    up = contract._rollup(long, mapping).sort('sample')
    assert up['intensity'].to_list() == [3.0, 5.0]
    assert up['status'].to_list() == ['detected', 'transferred']


def test_write_produces_the_contract_files(tmp_path):
    levels, runs = contract.build(CONSENSUS, IDENTIFICATIONS)
    design = session.design_frame(runs, {runs[0]: ('A', 1), runs[1]: ('B', 1)})
    contract.write(tmp_path, levels, design, {'tool': 'FLASHQuantConsensus'})
    names = {p.name for p in tmp_path.iterdir()}
    assert {'quant.json', 'design.tsv'} <= names
    for level in levels:
        assert f'{level}.tsv' in names or any(n.startswith(level) for n in names)


@pytest.fixture
def stored(tmp_path):
    file_manager = FileManager(tmp_path, Path(tmp_path, 'cache'))
    levels, runs = contract.build(CONSENSUS, IDENTIFICATIONS)
    design = {runs[0]: ('control', 1), runs[1]: ('treated', 1)}
    session.store(file_manager, 'comparison_20260101-120000', levels, runs, design, {}, 'flashtnt')
    return file_manager


def test_store_lists_the_comparison_under_its_tool(stored):
    assert session.comparisons(stored) == ['comparison_20260101-120000']
    assert session.comparisons(stored, 'flashtnt') == ['comparison_20260101-120000']
    assert session.comparisons(stored, 'flashdeconv') == []
    assert session.is_comparison(stored, 'comparison_20260101-120000')


def test_session_reads_back_design_and_level(stored):
    quant = session.QuantSession(stored, 'comparison_20260101-120000')
    assert quant.levels == ['feature', 'proteoform', 'protein']
    assert quant.level == 'proteoform'
    assert quant.conditions() == {
        'FLASHQuantConsensus_input1': 'control', 'FLASHQuantConsensus_input2': 'treated',
    }


def test_hidden_runs_are_left_out_and_come_back(stored):
    quant = session.QuantSession(stored, 'comparison_20260101-120000')
    first, second = quant.all_runs
    stamp = quant._stamp(quant.level)
    quant.set_hidden([second, 'not a run'])
    assert quant.hidden == [second]
    assert quant.runs == [first]
    assert quant.all_runs == [first, second]
    assert quant.conditions() == {first: 'control'}
    # Leaving a run out changes what the comparison chain computes
    assert quant._stamp(quant.level) != stamp
    quant.set_hidden([])
    assert quant.runs == [first, second]
    assert quant.conditions() == {first: 'control', second: 'treated'}


def test_two_runs_cannot_be_tested(stored):
    quant = session.QuantSession(stored, 'comparison_20260101-120000')
    assert 'three runs' in quant.problem()
    assert quant.analyse() == quant.problem()


def test_design_of_matches_runs_by_base_name(tmp_path):
    selected = ['/data/A_1.raw', '/data/B_1.mzML']
    design = ExperimentalDesign.default_design(selected)
    design['condition'] = ['control', 'treated']
    ExperimentalDesign.save_design(tmp_path, design, selected)
    mzmls = [tmp_path / 'converted-mzML' / 'A_1.mzML', '/data/B_1.mzML']
    assert run_name(mzmls[0]) == 'A_1'
    assert design_of(tmp_path, selected, mzmls) == {'A_1': ('control', 1), 'B_1': ('treated', 2)}


def test_a_proteoform_is_named_by_protein_residues_and_mass_shifts(tmp_path):
    fasta = tmp_path / 'db.fasta'
    fasta.write_text('>sp|P1|PROT_ECOLI Protein\nMKPEPTIDESEQ\nUENCE\n')
    sequences = contract.read_sequences(fasta)
    assert sequences == {'sp|P1|PROT_ECOLI': 'MKPEPTIDESEQUENCE'}

    form = contract.describe_proteoform('sp|P1|PROT_ECOLI', 'K[Acetyl]PEP(TIDE)[-17.02]SEQ', 1500.0, sequences)
    assert form['label'] == 'PROT_ECOLI 2–12 (Acetyl, −17.02)'
    assert (form['start'], form['end']) == (2, 12)
    assert form['shifts'] == 'Acetyl at 2; −17.02 in 6–9'

    # Not in the database: the start of its sequence stands for the residues
    assert contract.describe_proteoform('sp|P2|OTHER_HUMAN', 'ACDEFGHIKLMN[+15.99]', 0, {})['label'] == \
        'OTHER_HUMAN ACDEFGHIKL… (+15.99)'
    # No ProForma: the mass, or that it was not determined
    assert contract.describe_proteoform('sp|P1|PROT_ECOLI', None, 1234.56, sequences)['label'] == 'PROT_ECOLI 1234.6 Da'
    assert contract.describe_proteoform('sp|P1|PROT_ECOLI', None, -1.0, sequences)['label'] == \
        'PROT_ECOLI mass not determined'
