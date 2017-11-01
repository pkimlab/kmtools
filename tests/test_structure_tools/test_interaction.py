import logging
import os.path as op

import numpy as np
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

# #############################################################################
# Test get interactions
# #############################################################################


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


def test_get_interaction_cutoffs():
    """Test ``interactions.py`` using structures that we generate programatically."""
    structure_file = op.join(TEST_FILES_DIR, 't001.pdb')
    structure = kmbio.PDB.load(structure_file)
    # Change atom coordinates such that the distance from residue 0 to residue 1 is 1.005,
    # the distance from residue 1 to residue 2 is 2.005, etc.
    for chain_id in ['A', 'B']:
        coord = 0
        for i, residue in enumerate(structure[0][chain_id].residues):
            for atom in residue:
                atom.coord = np.array([coord, coord, coord])
            r_cutoff = i + 1
            coord += (r_cutoff**2 / 3)**.5 + 0.005
    # Assert that there are no core interactions past this cutoff
    interaction_df = structure_tools.get_interactions(structure, 1.0, interchain=False)
    assert interaction_df.empty
    # Each residue interacts with exactly one residue on the opposite chain
    interaction_df = structure_tools.get_interactions(structure, 1.0, interchain=True)
    assert len(interaction_df) == len(list(structure[0].residues))
    # Make sure that those interactions are actually interface interactions
    interactions_core, interactions_interface = structure_tools.process_interactions(
        interaction_df)
    assert len(interactions_core) == 0
    assert len(interactions_interface) == len(interaction_df)

    # Make sure that we only interact with one successive residue
    for residue_idx, r_cutoff in enumerate([1.01, 2.01, 3.01, 4.01]):
        interaction_df = structure_tools.get_interactions(structure, r_cutoff, interchain=False)
        interaction_df = \
            interaction_df[
                (interaction_df['chain_id_1'] == 'A') &
                (interaction_df['residue_idx_1'] == residue_idx) &
                (interaction_df['residue_idx_2'] >= interaction_df['residue_idx_1'])
            ]
        assert len(interaction_df) == 1

    # Test `process_interactions_{core,interface}` and `drop_duplicates_{core,interface}`
    interaction_df = structure_tools.get_interactions(structure, 2, interchain=True)
    interactions_core, interactions_interface = structure_tools.process_interactions(
        interaction_df)
    # # Core
    interactions_core_aggbychain = structure_tools.process_interactions_core(
        structure, interactions_core)
    assert len(interactions_core_aggbychain) == 2
    interactions_core, interactions_core_aggbychain = structure_tools.drop_duplicates_core(
        interactions_core, interactions_core_aggbychain)
    assert len(interactions_core_aggbychain) == 1
    # # Interface
    interactions_interface_aggbychain = structure_tools.process_interactions_interface(
        structure, interactions_interface)
    assert len(interactions_interface_aggbychain) == 2
    interactions_interface, interactions_interface_aggbychain = \
        structure_tools.drop_duplicates_interface(
            interactions_interface, interactions_interface_aggbychain)
    assert len(interactions_interface_aggbychain) == 1


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
    df = structure_tools.get_interactions(structure, r_cutoff=5, interchain=interchain)
    logger.info(df.head())
    # Test that get_interactions forks after process_structure
    structure = structure_tools.process_structure(structure)
    df = structure_tools.get_interactions(structure, r_cutoff=5, interchain=interchain)
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
    structure = kmbio.PDB.load('rcsb://{}.cif'.format(pdb_id), bioassembly_id=bioassembly_id)
    df = structure_tools.get_interactions(structure, r_cutoff=5.0)
    _assert_symmetrical(df)


def _assert_symmetrical(df):
    df['residue_pair_1'] = df[['residue_id_1', 'residue_id_2']].apply(tuple, axis=1)
    df['residue_pair_2'] = df[['residue_id_2', 'residue_id_1']].apply(tuple, axis=1)
    set1 = set(df[['model_id_1', 'chain_id_1', 'residue_pair_1']].apply(tuple, axis=1))
    set2 = set(df[['model_id_2', 'chain_id_2', 'residue_pair_2']].apply(tuple, axis=1))
    assert not set1 ^ set2


# #############################################################################
# Test process interactions
# #############################################################################


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


@pytest.mark.parametrize(
    "structure_file, interaction_file_core, processed_interaction_file_core_",
    [(d['structure_file'], d['interaction_file_core'], d['processed_interaction_file_core_'])
     for d in TEST_DATA['test_process_interactions_core']])
def test_process_interactions_core(structure_file, interaction_file_core,
                                   processed_interaction_file_core_):
    structure_file = op.join(TEST_FILES_DIR, structure_file)
    interaction_file_core = op.join(TEST_FILES_DIR, interaction_file_core)
    processed_interaction_file_core_ = op.join(TEST_FILES_DIR, processed_interaction_file_core_)

    structure = kmbio.PDB.load(structure_file)
    interactions_core = pd.read_csv(interaction_file_core)
    processed_interactions_core = structure_tools.process_interactions_core(
        structure, interactions_core)
    processed_interactions_core.drop(pd.Index(['residue_pair']), axis=1, inplace=True)

    processed_interactions_core_ = pd.read_csv(processed_interaction_file_core_)
    processed_interactions_core_.drop(pd.Index(['residue_pair']), axis=1, inplace=True)

    assert processed_interactions_core.equals(processed_interactions_core_)


@pytest.mark.parametrize(
    "structure_file, interaction_file_interface, processed_interaction_file_interface_",
    [(d['structure_file'], d['interaction_file_interface'],
      d['processed_interaction_file_interface_'])
     for d in TEST_DATA['test_process_interactions_interface']])
def test_process_interactions_interface(structure_file, interaction_file_interface,
                                        processed_interaction_file_interface_):
    structure_file = op.join(TEST_FILES_DIR, structure_file)
    interaction_file_interface = op.join(TEST_FILES_DIR, interaction_file_interface)
    processed_interaction_file_interface_ = op.join(TEST_FILES_DIR,
                                                    processed_interaction_file_interface_)

    structure = kmbio.PDB.load(structure_file)
    interactions_interface = pd.read_csv(interaction_file_interface)
    processed_interactions_interface = structure_tools.process_interactions_interface(
        structure, interactions_interface)
    processed_interactions_interface.drop(pd.Index(['residue_pair']), axis=1, inplace=True)

    processed_interactions_interface_ = pd.read_csv(processed_interaction_file_interface_)
    processed_interactions_interface_.drop(pd.Index(['residue_pair']), axis=1, inplace=True)

    assert processed_interactions_interface.equals(processed_interactions_interface_)


# #############################################################################
# Test remove duplicates
# #############################################################################


@pytest.mark.parametrize(
    "interaction_file_core, interaction_file_core_aggbychain, kept_idxs_, kept_idxs_aggbychain_",
    [(d['interaction_file_core'], d['interaction_file_core_aggbychain'], d['kept_idxs_'],
      d['kept_idxs_aggbychain_']) for d in TEST_DATA['test_drop_duplicates_core']])
def test_drop_duplicates_core(interaction_file_core, interaction_file_core_aggbychain, kept_idxs_,
                              kept_idxs_aggbychain_):
    interaction_file_core = op.join(TEST_FILES_DIR, interaction_file_core)
    interaction_file_core_aggbychain = op.join(TEST_FILES_DIR, interaction_file_core_aggbychain)

    interaction_core = pd.read_csv(interaction_file_core)
    interaction_core_aggbychain = pd.read_csv(interaction_file_core_aggbychain)

    interaction_core_unique, interaction_core_aggbychain_unique = \
        structure_tools.drop_duplicates_core(interaction_core,
                                             interaction_core_aggbychain)

    # Make sure that remaining columns have correct indexes
    assert interaction_core_unique.index.tolist() == kept_idxs_
    assert interaction_core_aggbychain_unique.index.tolist() == kept_idxs_aggbychain_

    # TODO(AS): Add more checks here...


@pytest.mark.parametrize(
    "interaction_file_interface, interaction_file_interface_aggbychain, kept_idxs_, "
    "kept_idxs_aggbychain_",
    [(d['interaction_file_interface'], d['interaction_file_interface_aggbychain'], d['kept_idxs_'],
      d['kept_idxs_aggbychain_']) for d in TEST_DATA['test_drop_duplicates_interface']])
def test_drop_duplicates_interface(interaction_file_interface,
                                   interaction_file_interface_aggbychain, kept_idxs_,
                                   kept_idxs_aggbychain_):
    interaction_file_interface = op.join(TEST_FILES_DIR, interaction_file_interface)
    interaction_file_interface_aggbychain = op.join(TEST_FILES_DIR,
                                                    interaction_file_interface_aggbychain)

    interaction_interface = pd.read_csv(interaction_file_interface)
    interaction_interface_aggbychain = pd.read_csv(interaction_file_interface_aggbychain)

    interaction_interface_unique, interaction_interface_aggbychain_unique = \
        structure_tools.drop_duplicates_interface(interaction_interface,
                                                  interaction_interface_aggbychain)

    # Make sure that remaining columns have correct indexes
    assert interaction_interface_unique.index.tolist() == kept_idxs_
    assert interaction_interface_aggbychain_unique.index.tolist() == kept_idxs_aggbychain_

    # TODO(AS): Add more checks here...
