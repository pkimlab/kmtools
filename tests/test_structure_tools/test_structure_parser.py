import logging

import pandas as pd
import pytest

from conftest import DIFFICULT, MISSING, PDB_IDS, random_subset  # NO_RESNAME_ATTRIBUTE_PDBS
from kmtools import structure_tools

logger = logging.getLogger(__name__)


# @pytest.mark.parametrize("pdb_id", random_subset(NO_RESNAME_ATTRIBUTE_PDBS))
# def test_no_resname_attribute(pdb_id):
#     """Test for the ``'NoneType' has no resname attribute`` error."""
#     s = structure_tools.fetch_structure(pdb_id, 'cif', True)
#     structure_tools.process_structure(s)
#     assert s
#

@pytest.mark.parametrize("pdb_id, pdb_type, biounit", random_subset([
    (pdb_id, pdb_type, biounit)
    for pdb_id in PDB_IDS
    for pdb_type in ['pdb', 'cif']
    # TODO: enable testing for biounits
    for biounit in [False]
    if (pdb_id, pdb_type, biounit) not in MISSING
]))
def test_process_structure(pdb_id, pdb_type, biounit):
    structure = structure_tools.fetch_structure(pdb_id, pdb_type, biounit)
    structure_tools.process_structure(structure)
    # Test model ids
    for model_idx, model in enumerate(structure.values()):
        assert model_idx == model.id
        assert structure[model_idx] == model
    # Test chain ids
    for model in structure.values():
        for chain_idx, chain in enumerate(model.values()):
            assert chain_idx == structure_tools.CHAIN_IDS.index(chain.id)
            assert model.ix[chain_idx] == chain
    # Test residue ids
    for model in structure.values():
        for chain in model.values():
            for residue_idx, residue in enumerate(chain.values()):
                assert residue_idx == residue.id[1]
                assert chain.ix[residue_idx] == residue


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
    structure = structure_tools.fetch_structure(pdb_id, pdb_type, biounit)
    df = structure_tools.get_interactions(structure, interchain)
    logger.info(df.head())


@pytest.mark.parametrize("pdb_id, biounit, interchain", random_subset([
    (pdb_id, biounit, interchain)
    for pdb_id in PDB_IDS
    # TODO: enable testing for biounits
    for biounit in [False]
    for interchain in [False, True]
    if (pdb_id, 'pdb', biounit) not in MISSING
    if pdb_id not in DIFFICULT
]))
def test_get_interactions_2(pdb_id, biounit, interchain):
    s_1 = structure_tools.fetch_structure(pdb_id, 'pdb', biounit)
    s_2 = structure_tools.fetch_structure(pdb_id, 'cif', biounit)

    # df_1 = structure_tools.get_interactions(s_1, interchain)
    # df_2 = structure_tools.get_interactions(s_2, interchain)
    #
    # with pytest.raises(AssertionError):
    #     pd.util.testing.assert_frame_equal(df_1, df_2)
    #
    structure_tools.process_structure(s_1)
    structure_tools.process_structure(s_2)

    del s_1.ix[-1]
    del s_2.ix[-1]

    df_1 = structure_tools.get_interactions(s_1, interchain)
    df_2 = structure_tools.get_interactions(s_2, interchain)

    pd.util.testing.assert_frame_equal(df_1, df_2)


@pytest.mark.parametrize("kwargs", random_subset([
    {'pdb_id': pdb_id, 'biounit': biounit}
    for pdb_id in PDB_IDS
    # TODO: enable testing for biounits
    for biounit in [False]
    if (pdb_id, 'pdb', biounit) not in MISSING
    if pdb_id not in DIFFICULT
]))
def test_get_interactions_symmetrical(kwargs):
    """Make sure that the output of `get_interactions` is symmetrical.

    That is, for each (model_id_1, chain_id_1, residue_id_1)
    there is a corresponding (model_id_2, chain_id_2, residue_id_2).
    """
    logger.debug(kwargs)

    def _assert_symmetrical(df):
        df['residue_pair_1'] = df[['residue_id_1', 'residue_id_2']].apply(tuple, axis=1)
        df['residue_pair_2'] = df[['residue_id_2', 'residue_id_1']].apply(tuple, axis=1)
        set1 = set(df[['model_id_1', 'chain_id_1', 'residue_pair_1']].apply(tuple, axis=1))
        set2 = set(df[['model_id_2', 'chain_id_2', 'residue_pair_2']].apply(tuple, axis=1))
        assert not set1 ^ set2

    s = structure_tools.fetch_structure(**kwargs)
    df = structure_tools.get_interactions(s)
    _assert_symmetrical(df)

    structure_tools.process_structure(s)
    df = structure_tools.get_interactions(s)
    _assert_symmetrical(df)
