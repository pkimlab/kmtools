import logging
import os

import pytest

from conftest import (
    LOCAL_REMOTE_MISMATCH, MISSING, PDB_IDS, ATOM_DEFINED_TWICE_PDBS, NO_RESNAME_ATTRIBUTE_PDBS)
from kmtools import structure_tools

logger = logging.getLogger(__name__)


# @pytest.mark.skipif(
#     'PDB_DATABASE_DIR' not in os.environ,
#     reason="set PDB_DATABASE_DIR environment variable to run this test!")
# @pytest.mark.parametrize("pdb_id, pdb_type, biounit", [
#     (pdb_id, pdb_type, biounit)
#     for pdb_id in PDB_IDS
#     for pdb_type in ['pdb', 'cif']
#     for biounit in [False, True]
#     if (pdb_id, pdb_type, biounit) not in MISSING
#     if (pdb_id, pdb_type, biounit) not in LOCAL_REMOTE_MISMATCH
#     if (pdb_type, biounit) != ('cif', True)
# ])
# def test_local_vs_remote(pdb_id, pdb_type, biounit):
#     """Make sure that loading local and remote files produces the same result."""
#     s1 = structure_tools.fetch_structure(pdb_id, pdb_type, biounit)
#     url = structure_tools.get_wwpdb_url(pdb_id, pdb_type, biounit, os.environ['PDB_DATABASE_DIR'])
#     s2 = structure_tools.load_structure(url, pdb_id, pdb_type)
#     assert structure_tools.allequal(s1, s2)
#
#
# @pytest.mark.parametrize("pdb_id", PDB_IDS)
# def test_allequal(pdb_id):
#     """Make sure that structures that should be equal are equal."""
#     for pdb_type in ['cif', 'pdb']:
#         for biounit in [False, True]:
#             s1 = structure_tools.fetch_structure(pdb_id)
#             s2 = structure_tools.fetch_structure(pdb_id)
#             assert structure_tools.allequal(s1, s2)
#
#
# @pytest.mark.parametrize("pdb_id_1, pdb_id_2", [
#     (pdb_id_1, pdb_id_2)
#     for pdb_id_1 in PDB_IDS
#     for pdb_id_2 in PDB_IDS
#     if pdb_id_1 != pdb_id_2
# ][:10])
# def test_allnotequal(pdb_id_1, pdb_id_2):
#     """Make sure that structures that should be different are different."""
#     for pdb_type in ['cif', 'pdb']:
#         for biounit in [False, True]:
#             s1 = structure_tools.fetch_structure(pdb_id_1)
#             s2 = structure_tools.fetch_structure(pdb_id_2)
#             assert not structure_tools.allequal(s1, s2)


# @pytest.mark.parametrize("pdb_id", ATOM_DEFINED_TWICE_PDBS)
# def test_atom_defined_twice(pdb_id):
#     """Tests for the ``Atom defined twice`` error."""
#     s = structure_tools.fetch_structure(pdb_id, 'cif', False)
#     assert s
#
#
@pytest.mark.parametrize("pdb_id", NO_RESNAME_ATTRIBUTE_PDBS)
def test_no_resname_attribute(pdb_id):
    """Test for the ``'NoneType' has no resname attribute`` error."""
    s = structure_tools.fetch_structure(pdb_id, 'cif', True)
    ps = structure_tools.process_structure(s)
    assert ps
