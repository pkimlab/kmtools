import logging

import pytest

from conftest import MISSING, NO_RESNAME_ATTRIBUTE_PDBS, PDB_IDS, random_subset
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


@pytest.mark.parametrize(
    "pdb_id, pdb_type, biounit",
    random_subset([
        (pdb_id, pdb_type, biounit)
        for pdb_id in PDB_IDS for pdb_type in ['pdb', 'cif']
        # TODO: enable testing for biounits
        for biounit in [False] if (pdb_id, pdb_type, biounit) not in MISSING
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
