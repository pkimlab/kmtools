"""
More comprehensive tests using complete PDB structures.
"""
import os.path as op

import kmbio.PDB
import pytest
from .conftest import DIFFICULT, MISSING, PDB_IDS, random_subset

from kmtools import structure_tools

TEST_FILES_DIR = op.abspath(op.splitext(__file__)[0])


@pytest.mark.parametrize("pdb_id", PDB_IDS)
@pytest.mark.parametrize("pdb_type", ['pdb', 'cif'])
@pytest.mark.parametrize("bioassembly_id", [0, 1])
@pytest.mark.parametrize("interchain", [False, True])
def test_get_interactions_pdb(pdb_id, pdb_type, bioassembly_id, interchain):
    """Test the entire pipeline using real structure files as input."""
    if pdb_id in DIFFICULT or (pdb_id, pdb_type, bioassembly_id) in MISSING:
        pytest.skip("Difficult or missing from RCSB...")

    structure = kmbio.PDB.load(
        'rcsb://{}.{}'.format(pdb_id, pdb_type), bioassembly_id=bioassembly_id)
    interactions = structure_tools.get_interactions(structure, r_cutoff=5, interchain=interchain)
    interactions_core, interactions_interface = structure_tools.process_interactions(interactions)

    # Group interactions by chain / chan pair
    interactions_core_aggbychain = structure_tools.process_interactions_core(
        structure, interactions_core)
    interactions_interface_aggbychain = structure_tools.process_interactions_interface(
        structure, interactions_interface)

    # Drop duplicate rows
    interactions_core, interactions_core_aggbychain = \
        structure_tools.drop_duplicates_core(interactions_core, interactions_core_aggbychain)
    interactions_interface, interactions_interface_aggbychain = \
        structure_tools.drop_duplicates_interface(interactions_interface,
                                                  interactions_interface_aggbychain)


GET_INTERACTIONS_SYMMETRICAL_DATA = [(pdb_id, bioassembly_id)
                                     for pdb_id in PDB_IDS for bioassembly_id in [0, 1]
                                     if (pdb_id, 'pdb', bool(bioassembly_id)) not in MISSING
                                     if pdb_id not in DIFFICULT]


@pytest.mark.parametrize("pdb_id, bioassembly_id", random_subset(GET_INTERACTIONS_SYMMETRICAL_DATA))
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
