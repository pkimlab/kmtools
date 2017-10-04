import logging

import pandas as pd
import pytest

from conftest import (DIFFICULT, MISSING, PDB_IDS, NO_RESNAME_ATTRIBUTE_PDBS,
                      random_subset)
from kmbio.PDB import load
from kmtools import structure_tools

logger = logging.getLogger(__name__)


# @pytest.mark.xfail(reason="process_structure seems to be terribly broken...")
@pytest.mark.parametrize("pdb_id", random_subset(NO_RESNAME_ATTRIBUTE_PDBS))
def test_no_resname_attribute(pdb_id):
    """Test for the ``'NoneType' has no resname attribute`` error."""
    s = load('rcsb://{}.{}'.format(pdb_id, 'cif'), bioassembly_id=1)
    s = structure_tools.process_structure(s)
    assert s


@pytest.mark.parametrize("pdb_id, pdb_type, biounit", random_subset([
    (pdb_id, pdb_type, biounit)
    for pdb_id in PDB_IDS
    for pdb_type in ['pdb', 'cif']
    # TODO: enable testing for biounits
    for biounit in [False]
    if (pdb_id, pdb_type, biounit) not in MISSING
]))
def test_process_structure(pdb_id, pdb_type, biounit):
    structure = load('rcsb://{}.{}'.format(pdb_id, pdb_type), bioassembly_id=int(biounit))
    structure = structure_tools.process_structure(structure)
    # Test model ids
    for model_idx, model in enumerate(structure):
        assert model_idx == model.id
        assert structure[model_idx] == model
    # Test chain ids
    for model in structure:
        for chain_idx, chain in enumerate(model):
            assert chain_idx == structure_tools.CHAIN_IDS.index(chain.id)
    # Test residue ids
    for model in structure:
        for chain in model:
            for residue_idx, residue in enumerate(chain):
                assert residue_idx == residue.id[1]


@pytest.mark.parametrize("pdb_id, pdb_type, biounit, interchain", random_subset([
    (pdb_id, pdb_type, biounit, interchain)
    for pdb_id in PDB_IDS
    for pdb_type in ['pdb', 'cif']
    # TODO: enable testing for biounits
    for biounit in [False]
    for interchain in [False, True]
    if (pdb_id, pdb_type, biounit) not in MISSING
]))
def test_get_interactions(pdb_id, pdb_type, biounit, interchain):
    structure = load('rcsb://{}.{}'.format(pdb_id, pdb_type), bioassembly_id=int(biounit))
    logger.debug('atoms: %s', list(structure.atoms))
    # Test 1
    df = structure_tools.get_interactions(structure, interchain)
    logger.info(df.head())
    # Test that get_interactions forks after process_structure
    structure = structure_tools.process_structure(structure)
    df = structure_tools.get_interactions(structure, interchain)
    logger.info(df.head())


@pytest.mark.parametrize("pdb_id, bioassembly_id", random_subset([
    (pdb_id, bioassembly_id)
    for pdb_id in PDB_IDS
    for bioassembly_id in [0, 1]
    if (pdb_id, 'pdb', bool(bioassembly_id)) not in MISSING
    if pdb_id not in DIFFICULT
]))
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

    s = load('rcsb://{}.cif'.format(pdb_id), bioassembly_id=bioassembly_id)
    df = structure_tools.get_interactions(s)
    _assert_symmetrical(df)
