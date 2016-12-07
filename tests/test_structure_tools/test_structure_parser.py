import logging

import pytest
from pandas.util.testing import assert_frame_equal

from kmtools import structure_tools

logger = logging.getLogger(__name__)

PDB_IDS = [
    '4dkl', '1arr', '1dvf', '3mbp',
]

#: There are no precalculated biounit structures for PDBs which have no biounits...
MISSING_PDB_BIOUNITS = [
    # (pdb_id, pdb_type, biounit)
    ('1arr', 'pdb', True),
]


# @pytest.mark.parametrize("pdb_id, pdb_type, biounit", [
#     (pdb_id, pdb_type, biounit)
#     for pdb_id in PDB_IDS
#     for pdb_type in ['pdb', 'cif']
#     for biounit in [False, True]
#     if (pdb_id, pdb_type, biounit) not in MISSING_PDB_BIOUNITS
# ])
# def test_parse_structure(pdb_id, pdb_type, biounit):
#     structure = structure_tools.fetch_structure(pdb_id, pdb_type, biounit)
#     structure = structure_tools.parse_structure(structure)
#     # Test model ids
#     for model_idx, model in enumerate(structure):
#         assert model_idx == model.id
#         assert structure[model_idx] == model
#     # Test chain ids
#     for model in structure:
#         for chain_idx, chain in enumerate(model):
#             assert chain_idx == structure_tools.CHAIN_IDS.index(chain.id)
#             assert model.child_list[chain_idx] == chain
#     # Test residue ids
#     for model in structure:
#         for chain in model:
#             for residue_idx, residue in enumerate(chain):
#                 assert residue_idx == residue.id[1]
#                 assert chain.child_list[residue_idx] == residue
#

@pytest.mark.parametrize("pdb_id, pdb_type, biounit, interchain", [
    (pdb_id, pdb_type, biounit, interchain)
    for pdb_id in PDB_IDS
    for pdb_type in ['pdb', 'cif']
    for biounit in [False, True]
    for interchain in [False, True]
    if (pdb_id, pdb_type, biounit) not in MISSING_PDB_BIOUNITS
])
def test_get_interactions_between_chains(pdb_id, pdb_type, biounit, interchain):
    structure = structure_tools.fetch_structure(pdb_id, pdb_type, biounit)
    df = structure_tools.get_interactions_between_chains(structure, interchain)
    logger.info(df.head())


@pytest.mark.parametrize("pdb_id, biounit, interchain", [
    (pdb_id, biounit, interchain)
    for pdb_id in PDB_IDS
    for biounit in [False, True]
    for interchain in [False, True]
    if (pdb_id, 'pdb', biounit) not in MISSING_PDB_BIOUNITS
])
def test_get_interactions_between_chains_2(pdb_id, biounit, interchain):
    structure_1 = structure_tools.fetch_structure(pdb_id, 'pdb', biounit)
    structure_2 = structure_tools.fetch_structure(pdb_id, 'cif', biounit)
    df_1 = structure_tools.get_interactions_between_chains(structure_1, interchain)
    df_2 = structure_tools.get_interactions_between_chains(structure_2, interchain)
    assert assert_frame_equal(df_1, df_2)
