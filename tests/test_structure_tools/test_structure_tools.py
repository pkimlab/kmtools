import pytest

from kmtools import structure_tools


PDB_IDS = [
    '4dkl', '1arr', '1dvf', '3mbp',
    # '4p6f',
]


@pytest.mark.parametrize("pdb_id", PDB_IDS)
def test_allequal(pdb_id):
    for pdb_type in ['cif', 'pdb']:
        for biounit in [False, True]:
            s1 = structure_tools.fetch_structure(pdb_id)
            s2 = structure_tools.fetch_structure(pdb_id)
            assert structure_tools.allequal(s1, s2)


@pytest.mark.parametrize("pdb_id_1, pdb_id_2", [
    (pdb_id_1, pdb_id_2)
    for pdb_id_1 in PDB_IDS
    for pdb_id_2 in PDB_IDS
    if pdb_id_1 != pdb_id_2
])
def test_assert_notallequal(pdb_id_1, pdb_id_2):
    for pdb_type in ['cif', 'pdb']:
        for biounit in [False, True]:
            s1 = structure_tools.fetch_structure(pdb_id_1)
            s2 = structure_tools.fetch_structure(pdb_id_2)
            assert not structure_tools.allequal(s1, s2)
