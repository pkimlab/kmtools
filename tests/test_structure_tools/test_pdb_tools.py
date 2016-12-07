import os
import logging

# from pandas.util.testing import assert_frame_equal

import pytest

from kmtools import structure_tools

logger = logging.getLogger(__name__)

PDB_IDS = [
    '4dkl', '1arr', '1dvf', '3mbp',
    '4p6f',
]

#: There are no precalculated mmCIF biounit structures for PDBs which have no biounits.
#: Conversely, PDBs which make up a large biounit often only have mmCIF structure.
MISSING = [
    # (pdb_id, pdb_type, biounit)
    ('1arr', 'pdb', True),
    ('4p6f', 'pdb', False),
    ('4p6f', 'pdb', True),
]


LOCAL_REMOTE_MISMATCH = [
    ('4dkl', 'pdb', False),
    ('4dkl', 'pdb', True),
]


@pytest.mark.skipif(
    'PDB_DATABASE_DIR' not in os.environ,
    reason="set PDB_DATABASE_DIR environment variable to run this test!")
@pytest.mark.parametrize("pdb_id, pdb_type, biounit", [
    (pdb_id, pdb_type, biounit)
    for pdb_id in PDB_IDS
    for pdb_type in ['pdb', 'cif']
    for biounit in [False, True]
    if (pdb_id, pdb_type, biounit) not in MISSING
    if (pdb_id, pdb_type, biounit) not in LOCAL_REMOTE_MISMATCH
    if (pdb_type, biounit) != ('cif', True)
])
def test_fetch_vs_load(pdb_id, pdb_type, biounit):
    s1 = structure_tools.fetch_structure(pdb_id, pdb_type, biounit)
    url = structure_tools.get_wwpdb_url(pdb_id, pdb_type, biounit, os.environ['PDB_DATABASE_DIR'])
    s2 = structure_tools.load_structure(url, pdb_id, pdb_type)
    assert structure_tools.allequal(s1, s2)
