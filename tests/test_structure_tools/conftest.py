import random

random.seed(42)


PDB_IDS = [
    '4dkl', '1arr', '1dvf', '3mbp',
    '4p6f', '1dcg'
]


#: There are no precalculated mmCIF biounit structures for PDBs which have no biounits.
#: Conversely, PDBs which make up a large biounit often only have mmCIF structure.
MISSING = [
    # (pdb_id, pdb_type, biounit)
    ('1arr', 'pdb', 1),
    ('4p6f', 'pdb', 0),
    ('4p6f', 'pdb', 1),
]


DIFFICULT = [
    '4p6f',
]


LOCAL_REMOTE_MISMATCH = [
    ('4dkl', 'pdb', False),
    ('4dkl', 'pdb', True),
]


# PDBs that cause errors
ATOM_DEFINED_TWICE_PDBS = [
    '2q3u', '2kax', '1wcn', '1wco', '2dii', '2eya',
]


NO_RESNAME_ATTRIBUTE_PDBS = [
    '1q3l', '4d1e', '1cty', '4pru', '1ctz', '2h9p',
]


NUM_TESTS_MAX = 4


def random_subset(lst, nmax=NUM_TESTS_MAX):
    random.shuffle(lst)
    return lst[:nmax]
