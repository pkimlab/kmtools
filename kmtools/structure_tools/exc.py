class ObjectsNotEqualError(Exception):
    pass


class PDBEmptySequenceError(Exception):
    """One of the sequences is missing from the alignment.

    The most likely cause is that the alignment domain definitions were incorrect.
    """
    pass


class PDBDomainDefsError(Exception):
    """PDB domain definitions not found in the pdb file."""
    pass


class MutationMismatchError(Exception):
    pass


class MSMSError(Exception):
    pass


class PDBError(Exception):
    pass


class PDBNotFoundError(Exception):
    pass


class PDBChainError(Exception):
    pass
