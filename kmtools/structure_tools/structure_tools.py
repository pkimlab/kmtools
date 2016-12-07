import Bio.PDB

__all__ = ['allequal']


def allequal(s1, s2):
    if type(s1) != type(s2):
        raise Exception
    if isinstance(s1, Bio.PDB.Atom.Atom):
        return s1 == s2
    return all(allequal(so1, so2) for (so1, so2) in zip(s1, s2))
