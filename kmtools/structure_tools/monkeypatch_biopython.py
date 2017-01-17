import numpy as np
import Bio.PDB
import Bio.Seq
import Bio.SeqRecord

from kmtools.structure_tools import CHAIN_IDS

# From . import * should not import anything (just monekypatch BioPython)
__all__ = []


# Make Bio objects hashable (hack!)
# TODO: Think more closely about which fields should be used to construct the hash.
# Bio.Seq.Seq.__hash__ = lambda self: hash(self.__repr__())
Bio.Seq.Seq.__eq__ = lambda self, other: hash(repr(self)) == hash(repr(other))
# Bio.SeqRecord.SeqRecord.__hash__ = lambda self: hash(self.__repr__())
Bio.SeqRecord.SeqRecord.__eq__ = lambda self, other: hash(repr(self)) == hash(repr(other))


# __hash__ and __eq__
# Bio.PDB.Atom.Atom.__hash__ = lambda self: hash(repr(self))
# Bio.PDB.Atom.Atom.__hash__ = lambda self: hash((self.name, self.coord.data.tobytes(), ))


# Set __repr__
def default_repr(self, defaults=None):
    """A `__repr__` method that makes more sense than what comes with BioPython.

    Inspiration comes from this `StackOverflow answer:
        <http://stackoverflow.com/a/7784214/2063031>`_

    Parameters
    ----------
    defaults : :cls:`dict`
        A dictionary of class arguments mapped to values.
        Useful if you want to override the defaults in `__init__`.

    Examples
    --------
    >>> Foo = type('Foo', (), {})
    >>> Foo.__repr__ = default_repr
    >>> Foo()
    Foo()
    """
    import inspect
    if defaults is None:
        defaults = {}
    cls = self.__class__
    name = cls.__name__
    params = list(inspect.signature(cls).parameters.keys())
    args = [(p, repr(defaults[p](self) if p in defaults else getattr(self, p))) for p in params]
    return "{}({})".format(name, ", ".join("{}={}".format(*a) for a in args))


def set_repr(cls, defaults=None):
    import inspect
    if defaults is None:
        defaults = {}
    name = cls.__name__
    params = list(inspect.signature(cls).parameters.keys())
    if defaults is None:
        get_attr = getattr
    else:
        get_attr = lambda self, p: defaults[p](self) if p in defaults else getattr(self, p)  # noqa
    cls.__repr__ = lambda self: (
        "{}({})".format(
            name,
            ", ".join("{}={}".format(p, repr(get_attr(self, p))) for p in params)))


atom_defaults = {
    'fullname': lambda self: getattr(self, 'name'),
    'serial_number': lambda self: None,
}
set_repr(Bio.PDB.Atom.Atom, atom_defaults)


for cls in [
    Bio.PDB.Residue.Residue,
    Bio.PDB.Chain.Chain,
    Bio.PDB.Model.Model,
    Bio.PDB.Structure.Structure,
]:
    set_repr(cls)


# Set __comp__
def set_atom_comp(cls):
    cls.__eq__ = lambda self, other: (
        self.name == other.name and np.allclose(self.coord, other.coord)
    )
    cls.__ne__ = lambda self, other: (
        self.name != other.name or not np.allclose(self.coord, other.coord)
    )


def set_residue_comp(cls):
    cls.__lt__ = lambda self, other: self.id[1] < other.id[1]
    cls.__le__ = lambda self, other: self.id[1] <= other.id[1]
    cls.__eq__ = lambda self, other: self.id == other.id and self.resname == other.resname
    cls.__ne__ = lambda self, other: self.id != other.id or self.resname != other.resname
    cls.__ge__ = lambda self, other: self.id[1] >= other.id[1]
    cls.__gt__ = lambda self, other: self.id[1] > other.id[1]


def set_chain_comp(cls):
    cls.__lt__ = lambda self, other: CHAIN_IDS.index(self.id) < CHAIN_IDS.index(other.id)
    cls.__le__ = lambda self, other: CHAIN_IDS.index(self.id) <= CHAIN_IDS.index(other.id)
    cls.__eq__ = lambda self, other: self.id == other.id
    cls.__ne__ = lambda self, other: self.id != other.id
    cls.__ge__ = lambda self, other: CHAIN_IDS.index(self.id) >= CHAIN_IDS.index(other.id)
    cls.__gt__ = lambda self, other: CHAIN_IDS.index(self.id) > CHAIN_IDS.index(other.id)


def set_model_comp(cls):
    cls.__lt__ = lambda self, other: self.id < other.id
    cls.__le__ = lambda self, other: self.id <= other.id
    cls.__eq__ = lambda self, other: self.id == other.id
    cls.__ne__ = lambda self, other: self.id != other.id
    cls.__ge__ = lambda self, other: self.id >= other.id
    cls.__gt__ = lambda self, other: self.id > other.id


def set_structure_comp(cls):
    cls.__lt__ = lambda self, other: self.id.lower() < other.id.lower()
    cls.__le__ = lambda self, other: self.id.lower() <= other.id.lower()
    cls.__eq__ = lambda self, other: self.id.lower() == other.id.lower()
    cls.__ne__ = lambda self, other: self.id.lower() != other.id.lower()
    cls.__ge__ = lambda self, other: self.id.lower() >= other.id.lower()
    cls.__gt__ = lambda self, other: self.id.lower() > other.id.lower()


set_atom_comp(Bio.PDB.Atom.Atom)
set_residue_comp(Bio.PDB.Residue.Residue)
set_chain_comp(Bio.PDB.Chain.Chain)
set_model_comp(Bio.PDB.Model.Model)
set_structure_comp(Bio.PDB.Structure.Structure)
