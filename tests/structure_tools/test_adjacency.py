import pytest

import kmbio.PDB
from kmtools import structure_tools


@pytest.mark.parametrize("pdb_id", ["4dkl"])
def test_get_distances(pdb_id):
    structure = kmbio.PDB.load(f"rcsb://{pdb_id}.cif")
    dd = structure_tools.DomainDef(0, "A", 10, 40)
    domain = structure_tools.extract_domain(structure, [dd])
    # BioPython's KDTree method
    interactions = structure_tools.get_interactions(domain, 5, add_reverse=False)
    interactions_core, interactions_interface = structure_tools.process_interactions(interactions)
    assert (interactions_core["residue_idx_1"] <= interactions_core["residue_idx_2"]).all()
    interactions_set_1 = set(
        interactions_core[["residue_idx_1", "residue_idx_2"]].apply(tuple, axis=1)
    )
    # MDAnalysis method (with distances)
    pairs_df = structure_tools.get_distances(domain, 5, groupby="residue")
    interactions_set_2 = set(pairs_df[["residue_idx_1", "residue_idx_2"]].apply(tuple, axis=1))
    #
    assert not interactions_set_1 ^ interactions_set_2
