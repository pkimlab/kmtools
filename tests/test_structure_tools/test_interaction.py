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


@pytest.mark.parametrize("structure_file, r_cutoff, interchain, interaction_file_", [
    tuple(d[key] for key in 'structure_file, r_cutoff, interchain, interaction_file_'.split(', '))
    for d in TEST_DATA['test_get_interactions']
])
def test_get_interactions(structure_file, r_cutoff, interchain, interaction_file_):
    """Test `get_interactions` using contrived structure files."""
    structure_file = op.join(TEST_FILES_DIR, structure_file)
    interaction_file_ = op.join(TEST_FILES_DIR, interaction_file_)
    structure = kmbio.PDB.load(structure_file)
    interaction_df = structure_tools.get_interactions(structure, r_cutoff, interchain)
    interaction_df_ = pd.read_csv(interaction_file_)
    assert interaction_df.equals(interaction_df_)


@pytest.mark.parametrize("pdb_id", PDB_IDS)
@pytest.mark.parametrize("pdb_type", ['pdb', 'cif'])
@pytest.mark.parametrize("biounit", [False])
@pytest.mark.parametrize("interchain", [False, True])
def test_get_interactions_pdb(pdb_id, pdb_type, biounit, interchain):
    """Test `get_interactions` using real structure files."""
    if pdb_id in DIFFICULT or (pdb_id, pdb_type, biounit) in MISSING:
        pytest.skip("Difficult or missing from RCSB...")

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


GET_INTERACTIONS_SYMMETRICAL_DATA = [(pdb_id, bioassembly_id)
                                     for pdb_id in PDB_IDS for bioassembly_id in [0, 1]
                                     if (pdb_id, 'pdb', bool(bioassembly_id)) not in MISSING
                                     if pdb_id not in DIFFICULT]


@pytest.mark.parametrize("pdb_id, bioassembly_id",
                         random_subset(GET_INTERACTIONS_SYMMETRICAL_DATA))
def test_get_interactions_symmetrical(pdb_id, bioassembly_id):
    """Make sure that the output of `get_interactions` is symmetrical.

    That is, for each (model_id_1, chain_id_1, residue_id_1)
    there is a corresponding (model_id_2, chain_id_2, residue_id_2).
    """
    s = kmbio.PDB.load('rcsb://{}.cif'.format(pdb_id), bioassembly_id=bioassembly_id)
    df = structure_tools.get_interactions(s)
    _assert_symmetrical(df)


def _assert_symmetrical(df):
    df['residue_pair_1'] = df[['residue_id_1', 'residue_id_2']].apply(tuple, axis=1)
    df['residue_pair_2'] = df[['residue_id_2', 'residue_id_1']].apply(tuple, axis=1)
    set1 = set(df[['model_id_1', 'chain_id_1', 'residue_pair_1']].apply(tuple, axis=1))
    set2 = set(df[['model_id_2', 'chain_id_2', 'residue_pair_2']].apply(tuple, axis=1))
    assert not set1 ^ set2


@pytest.mark.parametrize(
    "interaction_file, interaction_file_core_, interaction_file_interface_",
    [(d['interaction_file'], d['interaction_file_core_'], d['interaction_file_interface_'])
     for d in TEST_DATA['test_process_interactions']])
def test_process_interactions(interaction_file, interaction_file_core_,
                              interaction_file_interface_):
    # Process interactions
    interaction_file = op.join(TEST_FILES_DIR, interaction_file)
    interactions = pd.read_csv(interaction_file)
    interactions_core, interactions_interface = structure_tools.process_interactions(interactions)
    assert len(interactions) == len(interactions_core) + len(interactions_interface)

    # Make sure the core interactions are correct
    interaction_file_core_ = op.join(TEST_FILES_DIR, interaction_file_core_)
    interactions_core_ = pd.read_csv(interaction_file_core_)
    assert interactions_core.reset_index(drop=True).equals(interactions_core_)

    # Make sure the interface interactions are correct
    if interactions_interface.empty:
        assert interaction_file_interface_ is None
    else:
        interaction_file_interface_ = op.join(TEST_FILES_DIR, interaction_file_interface_)
        interactions_interface_ = pd.read_csv(interaction_file_interface_)
        assert interactions_interface.reset_index(drop=True).equals(interactions_interface_)


def test_process_interactions_core():
    ...


def test_process_interactions_interface():
    ...
