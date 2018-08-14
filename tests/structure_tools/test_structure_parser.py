from pathlib import Path

import kmbio.PDB
import pytest
from kmbio.PDB import Model, Structure

from kmtools import structure_tools

from .conftest import MISSING, NO_RESNAME_ATTRIBUTE_PDBS, PDB_IDS, random_subset

TESTS_DIR = Path(__file__).absolute().parent


# @pytest.mark.xfail(reason="process_structure seems to be terribly broken...")
@pytest.mark.parametrize("pdb_id", random_subset(NO_RESNAME_ATTRIBUTE_PDBS))
def test_no_resname_attribute(pdb_id):
    """Test for the ``'NoneType' has no resname attribute`` error."""
    s = kmbio.PDB.load("rcsb://{}.{}".format(pdb_id, "cif"), bioassembly_id=1)
    s = structure_tools.process_structure(s)
    assert s


@pytest.mark.parametrize(
    "pdb_id, pdb_type, biounit",
    random_subset(
        [
            (pdb_id, pdb_type, biounit)
            for pdb_id in PDB_IDS
            for pdb_type in ["pdb", "cif"]
            # TODO: enable testing for biounits
            for biounit in [False]
            if (pdb_id, pdb_type, biounit) not in MISSING
        ]
    ),
)
def test_process_structure(pdb_id, pdb_type, biounit):
    structure = kmbio.PDB.load("rcsb://{}.{}".format(pdb_id, pdb_type), bioassembly_id=int(biounit))
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


@pytest.mark.parametrize(
    "structure_file, sequence_chains, hetatm_chains",
    [(TESTS_DIR.joinpath("structures", "1yf4.cif"), ["A", "B"], ["C", "D", "E"])],
)
def test_copy_hetatm_chain(structure_file, sequence_chains, hetatm_chains):
    structure = kmbio.PDB.load(structure_file)
    new_structure = Structure(structure.id, [Model(0)])
    for chain_id in sequence_chains:
        new_structure[0].add(structure[0][chain_id])
    for chain_id in hetatm_chains:
        old_chain = structure[0][chain_id]
        new_chain = structure_tools.copy_hetatm_chain(new_structure, old_chain, r_cutoff=5)
        assert len(list(old_chain.residues)) == len(list(new_chain.residues))


@pytest.mark.parametrize(
    "structure_file, hetatm_identity, hetatm_flags",
    [(TESTS_DIR.joinpath("structures", "1yf4.cif"), 0.9, [False, False, True, True, True])],
)
def test_chain_is_hetatm(structure_file, hetatm_identity, hetatm_flags):
    structure = kmbio.PDB.load(structure_file)
    hetatm_flags_ = [structure_tools.chain_is_hetatm(chain) for chain in structure[0].chains]
    assert hetatm_flags == hetatm_flags_
