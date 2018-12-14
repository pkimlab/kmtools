import os.path as op

import kmbio.PDB
import pytest
from ruamel import yaml

from kmtools import structure_tools

TEST_FILES_DIR = op.abspath(op.splitext(__file__)[0])

with open(op.join(TEST_FILES_DIR, 'test_data.yml'), 'r') as fin:
    TEST_DATA = yaml.load(fin)


@pytest.mark.parametrize("structure_file, model_id, chain_id, aa_sequence_", [
    tuple(d[key] for key in 'structure_file, model_id, chain_id, aa_sequence_'.split(', '))
    for d in TEST_DATA['test_extract_sequence']
])
def test_extract_aa_sequence(structure_file, model_id, chain_id, aa_sequence_):
    structure_file = op.join(TEST_FILES_DIR, structure_file)
    structure = kmbio.PDB.load(structure_file)
    aa_sequence = structure_tools.extract_aa_sequence(structure, model_id, chain_id)
    assert aa_sequence == aa_sequence_


@pytest.mark.parametrize("structure_file, model_id, chain_id, residue_sequence_", [
    tuple(d[key] for key in 'structure_file, model_id, chain_id, residue_sequence_'.split(', '))
    for d in TEST_DATA['test_extract_sequence']
])
def test_extract_residue_sequence(structure_file, model_id, chain_id, residue_sequence_):
    structure_file = op.join(TEST_FILES_DIR, structure_file)
    structure = kmbio.PDB.load(structure_file)
    residue_sequence = structure_tools.extract_residue_sequence(structure, model_id, chain_id)
    assert residue_sequence == residue_sequence_


@pytest.mark.parametrize("residue_pair, residue_pair_hash_", [
    (('G', 'A'), '14965c5e8f0e81e3cfc839f6cfd73989'),
])
def test_hash_residue_pair(residue_pair, residue_pair_hash_):
    residue_pair_hash = structure_tools.hash_residue_pair(residue_pair)
    assert residue_pair_hash == residue_pair_hash_
    # Make sure it also works in reverse
    residue_pair = (residue_pair[1], residue_pair[0])
    residue_pair_hash = structure_tools.hash_residue_pair(residue_pair)
    assert residue_pair_hash == residue_pair_hash_
