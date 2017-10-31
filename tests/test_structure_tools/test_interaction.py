import logging
import os.path as op

import pandas as pd
import pytest
import yaml

import kmbio.PDB
from conftest import DIFFICULT, MISSING, PDB_IDS, random_subset
from kmtools import structure_tools

logger = logging.getLogger(__name__)

TEST_FILES_DIR = op.abspath(op.splitext(__file__)[0])

with open(op.join(TEST_FILES_DIR, 'test_data.yml'), 'r') as fin:
    TEST_DATA = yaml.load(fin)


@pytest.mark.parametrize(
    "structure_file, r_cutoff, interchain, interaction_file_",
    [tuple(d[key] for key in 'structure_file, r_cutoff, interchain, interaction_file_'.split(', '))
     for d in TEST_DATA['test_get_interactions']])
def test_get_interactions(structure_file, r_cutoff, interchain, interaction_file_):
    structure_file = op.join(TEST_FILES_DIR, structure_file)
    interaction_file_ = op.join(TEST_FILES_DIR, interaction_file_)
    structure = kmbio.PDB.load(structure_file)
    interaction_df = structure_tools.get_interactions(structure, r_cutoff, interchain)
    interaction_df_ = pd.read_csv(interaction_file_)
    assert interaction_df.equals(interaction_df_)


@pytest.mark.parametrize(
    "pdb_id, pdb_type, biounit, interchain",
    random_subset([
        (pdb_id, pdb_type, biounit, interchain)
        for pdb_id in PDB_IDS for pdb_type in ['pdb', 'cif']
        # TODO: enable testing for biounits
        for biounit in [False] for interchain in [False, True]
        if (pdb_id, pdb_type, biounit) not in MISSING
    ]))
def test_get_interactions_mock(pdb_id, pdb_type, biounit, interchain):
    structure = kmbio.PDB.load(
        'rcsb://{}.{}'.format(pdb_id, pdb_type), bioassembly_id=int(biounit))
    logger.debug('atoms: %s', list(structure.atoms))
    # Test 1
    df = structure_tools.get_interactions(structure, interchain=interchain)
    logger.info(df.head())
    # Test that get_interactions forks after process_structure
    structure = structure_tools.process_structure(structure)
    df = structure_tools.get_interactions(structure, interchain=interchain)
    logger.info(df.head())


@pytest.mark.parametrize("pdb_id, bioassembly_id",
                         random_subset([(pdb_id, bioassembly_id)
                                        for pdb_id in PDB_IDS for bioassembly_id in [0, 1]
                                        if (pdb_id, 'pdb', bool(bioassembly_id)) not in MISSING
                                        if pdb_id not in DIFFICULT]))
def test_get_interactions_symmetrical(pdb_id, bioassembly_id):
    """Make sure that the output of `get_interactions` is symmetrical.

    That is, for each (model_id_1, chain_id_1, residue_id_1)
    there is a corresponding (model_id_2, chain_id_2, residue_id_2).
    """

    def _assert_symmetrical(df):
        df['residue_pair_1'] = df[['residue_id_1', 'residue_id_2']].apply(tuple, axis=1)
        df['residue_pair_2'] = df[['residue_id_2', 'residue_id_1']].apply(tuple, axis=1)
        set1 = set(df[['model_id_1', 'chain_id_1', 'residue_pair_1']].apply(tuple, axis=1))
        set2 = set(df[['model_id_2', 'chain_id_2', 'residue_pair_2']].apply(tuple, axis=1))
        assert not set1 ^ set2

    s = kmbio.PDB.load('rcsb://{}.cif'.format(pdb_id), bioassembly_id=bioassembly_id)
    df = structure_tools.get_interactions(s)
    _assert_symmetrical(df)
