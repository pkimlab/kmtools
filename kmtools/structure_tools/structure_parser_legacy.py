# flake8: noqa
"""Structure Tools.

Tools for the analysis of PDB structures.

.. note::

    This whole module is deprecated and will be deleted when I am totally sure that I don't need anything else from it for ELASPIC.
"""
import logging
import os.path as op
import string
from collections import OrderedDict, defaultdict
from functools import wraps

import Bio
import numpy as np
import six
from Bio.Alphabet import IUPAC
from Bio.PDB import PDBIO, NeighborSearch, Select
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Seq import Seq

from kmtools import structure_tools
from kmtools.structure_tools import A_DICT, AAA_DICT, LYSINE_ATOMS, METHYLATED_LYSINES

logger = logging.getLogger(__name__)


def euclidean_distance(a, b):
    """Calculate the Euclidean distance between two lists or tuples of arbitrary length."""
    return np.sqrt(sum((a - b) ** 2 for a, b in zip(a, b)))


def calculate_distance(atom_1, atom_2, cutoff=None):
    """Calculate the distance between two points in 3D space.

    Parameters
    ----------
    cutoff : float, optional
        The maximum distance allowable between two points.
    """
    if (type(atom_1) == type(atom_2) == list) or (type(atom_1) == type(atom_2) == tuple):
        a = atom_1
        b = atom_2
    elif hasattr(atom_1, "coord") and hasattr(atom_2, "coord"):
        a = atom_1.coord
        b = atom_2.coord
    else:
        raise Exception("Unsupported format {} {}".format(type(atom_1), type(atom_2)))

    assert len(a) == 3 and len(b) == 3
    if cutoff is None or all(abs(p - q) <= cutoff for p, q in zip(a, b)):
        return euclidean_distance(a, b)


def get_chain_seqres_sequence(chain, aa_only=False):
    """Get the amino acid sequence for the construct coding for the given chain.

    Extracts a sequence from a PDB file. Usefull when interested in the
    sequence that was used for crystallization and not the ATOM sequence.

    Parameters
    ----------
    aa_only : bool
        If aa_only is set to `False`, selenomethionines will be included in the sequence.
        See: http://biopython.org/DIST/docs/api/Bio.PDB.Polypeptide-module.html.
    """
    sequence = Seq("", IUPAC.protein)
    for pb in PPBuilder().build_peptides(chain, aa_only=aa_only):
        sequence += sequence + pb.get_sequence()
    return sequence


def get_chain_sequence_and_numbering(chain, domain_def_tuple=None, include_hetatms=False):
    """Get the amino acid sequence and a list of residue ids for the given chain.

    Parameters
    ----------
    chain : Bio.PDB.Chain.Chain
        The chain for which to get the amino acid sequence and numbering.
    """
    if domain_def_tuple is not None:
        start_resid, end_resid = domain_def_tuple

    chain_numbering = []
    chain_numbering_extended = []
    chain_sequence = []
    inside_domain = False
    for res in chain:
        #
        resid = str(res.id[1]) + res.id[2].strip()

        if domain_def_tuple is None or resid == start_resid:
            inside_domain = True

        if inside_domain and (include_hetatms or res.resname in AAA_DICT):
            chain_numbering.append(res.id[1])
            chain_numbering_extended.append(resid)
            chain_sequence.append(AAA_DICT.get(res.resname, "."))

        if domain_def_tuple is not None and resid == end_resid:
            inside_domain = False

    chain_sequence = "".join(chain_sequence)
    return chain_sequence, chain_numbering_extended


def convert_position_to_resid(chain, positions, domain_def_tuple=None):
    """Convert mutation_domain to mutation_modeller.

    In mutation_modeller, the first amino acid in a chain may start
    with something other than 1.
    """
    __, chain_numbering = get_chain_sequence_and_numbering(chain, domain_def_tuple)
    logger.debug("chain_numbering: {}".format(chain_numbering))
    logger.debug("positions: {}".format(positions))
    return [chain_numbering[p - 1] for p in positions]


def get_structure_sequences(file_or_structure, seqres_sequence=False):
    """Return a dictionary of sequences for a given file or Structure.

    Parameters
    ----------
    file_or_structure : str | biopython.Structure | biopython.Model | biopython.Chain
        PDB filename or biopython object from which to extract the sequence.
    """
    if isinstance(file_or_structure, six.string_types):
        model = structure_tools.load_structure(file_or_structure)[0]
    elif isinstance(file_or_structure, Bio.PDB.Structure.Structure):
        model = file_or_structure[0]
    elif isinstance(file_or_structure, Bio.PDB.Model.Model):
        model = file_or_structure
    elif isinstance(file_or_structure, Bio.PDB.Chain.Chain):
        model = [file_or_structure]
    else:
        raise Exception

    chain_sequences = defaultdict(list)
    for chain in model:
        if seqres_sequence:
            chain_sequence = get_chain_seqres_sequence(chain)
        else:
            chain_sequence, __ = get_chain_sequence_and_numbering(chain)
        chain_sequences[chain.id] = chain_sequence
    return chain_sequences


def suppress_logger(fn):
    @wraps(fn)
    def fn_quiet(*args, **kwargs):
        level = logger.level
        logger.setLevel(logging.WARNING)
        try:
            return fn(*args, **kwargs)
        finally:
            logger.setLevel(level)

    return fn_quiet


def convert_aa(aa, quiet=False):
    """Convert amino acids from three letter code to one letter code or vice versa.

    .. note:: Deprecated!

       Use ``''.join(AAA_DICT[aaa] for aaa in aa)`` and ``''.join(A_DICT[a] for a in aa)``.
    """
    if quiet:
        return suppress_logger(convert_aa)(aa)

    if len(aa) == 3:
        try:
            return AAA_DICT[aa.upper()]
        except KeyError:
            if not quiet:
                logger.debug("Not a valid amino acid: {}".format(aa))
            return
    if len(aa) == 1:
        try:
            return A_DICT[aa.upper()]
        except KeyError:
            if not quiet:
                logger.debug("Not a valid amino acid: {}".format(aa))
            return
    if not quiet:
        logger.debug("Not a valid amino acid: {}".format(aa))


# STANDALONE FUNCTIONS
def get_interactions(model, chain_id, r_cutoff=6):
    """
    """
    interactions = {}
    for chain_id_2, chain_2 in model.child_dict.items():
        if chain_id == chain_id_2:
            continue
        interactions[chain_id_2] = get_interactions(model, chain_id, chain_id_2, r_cutoff)
    return {k: v for (k, v) in interactions.items() if v}


def get_interactions_between_chains(model, chain_id_1, chain_id_2, r_cutoff=6):
    """Calculate interactions between the residues of the two chains.

    An interaction is defines as a pair of residues where at least one pair of atom
    is closer than r_cutoff.

    .. deprecated:: 1.0
        Use python:fn:`get_interacting_residues` instead.
        It gives you both the residue index and the resnum.

    Returns
    -------
    OrderedDict
        Keys are (residue_number, residue_amino_acid) tuples
        (e.g. ('0', 'M'), ('1', 'Q'), ...).
        Values are lists of (residue_number, residue_amino_acid) tuples.
        (e.g. [('0', 'M'), ('1', 'Q'), ...]).
    """
    # Extract the chains of interest from the model
    chain_1 = None
    chain_2 = None
    for child in model.get_list():
        if child.id == chain_id_1:
            chain_1 = child
        if child.id == chain_id_2:
            chain_2 = child
    if chain_1 is None or chain_2 is None:
        raise Exception("Chains %s and %s were not found in the model" % (chain_id_1, chain_id_2))

    ns = NeighborSearch(list(chain_2.get_atoms()))
    interactions_between_chains = OrderedDict()
    for idx, residue_1 in enumerate(chain_1):
        if residue_1.resname in AAA_DICT and residue_1.id[0] == " ":
            resnum_1 = str(residue_1.id[1]) + residue_1.id[2].strip()
            resaa_1 = convert_aa(residue_1.get_resname(), quiet=True)
            interacting_residues = set()
            for atom_1 in residue_1:
                interacting_residues.update(ns.search(atom_1.get_coord(), r_cutoff, "R"))
            interacting_resids = []
            for residue_2 in interacting_residues:
                resnum_2 = str(residue_2.id[1]) + residue_2.id[2].strip()
                resaa_2 = convert_aa(residue_2.get_resname(), quiet=True)
                if residue_2.resname in AAA_DICT and residue_2.id[0] == " ":
                    interacting_resids.append((resnum_2, resaa_2))
            if interacting_resids:
                interacting_resids.sort(
                    key=lambda x: int("".join([c for c in x[0] if c.isdigit()]))
                )
                interactions_between_chains[(resnum_1, resaa_1)] = interacting_resids
    return interactions_between_chains


def get_interactions_between_chains_slow(model, pdb_chain_1, pdb_chain_2, r_cutoff=5):
    """Calculate interactions between residues in pdb_chain_1 and pdb_chain_2.

    An interaction is defines as a pair of residues where at least one pair of atom
    is closer than r_cutoff. The default value for r_cutoff is 5 Angstroms.

    .. deprecated:: 1.0
        Use :func:`get_interacting_residues` instead.
        It gives you both the residue index and the resnum.

    """
    # Extract the chains of interest from the model
    chain_1 = None
    chain_2 = None
    for child in model.get_list():
        if child.id == pdb_chain_1:
            chain_1 = child
        if child.id == pdb_chain_2:
            chain_2 = child
    if chain_1 is None or chain_2 is None:
        raise Exception("Chains %s and %s were not found in the model" % (pdb_chain_1, pdb_chain_2))

    interactions_between_chains = OrderedDict()
    for idx, residue_1 in enumerate(chain_1):
        if residue_1.resname in AAA_DICT and residue_1.id[0] == " ":
            resnum_1 = str(residue_1.id[1]) + residue_1.id[2].strip()
            resaa_1 = convert_aa(residue_1.get_resname())
            interacting_resids = []
            for residue_2 in chain_2:
                resnum_2 = str(residue_2.id[1]) + residue_2.id[2].strip()
                resaa_2 = convert_aa(residue_2.get_resname())
                r_min = None
                if residue_2.resname in AAA_DICT and residue_2.id[0] == " ":
                    for atom_1 in residue_1:
                        for atom_2 in residue_2:
                            r = calculate_distance(atom_1, atom_2, r_cutoff)
                            if r is not None:
                                if r_min and r < r_min:
                                    r_min = r
                                elif not r_min:
                                    r_min = r
                if r_min:
                    interacting_resids.append((resnum_2, resaa_2, r_min))
            if interacting_resids:
                interactions_between_chains[(resnum_1, resaa_1)] = interacting_resids
    return interactions_between_chains


def chain_is_hetatm(chain):
    """Return True if the chain is made up entirely of HETATMs."""
    hetatms = [None] * len(chain)
    for i in range(len(chain.child_list)):
        res = chain.child_list[i]
        hetatms[i] = res.resname not in AAA_DICT
    if all(hetatms):
        return True
    elif not any(hetatms):
        return False
    else:
        # Something went wrong.
        sequence, numbering = get_chain_sequence_and_numbering(chain)
        message = (
            "Some but not all residues in chain {} are hetatms!\n".format(chain.id)
            + "sequence: {}\n".format(sequence)
            + "numbering: {}\n".format(numbering)
        )
        logger.debug(message)
        False


def get_aa_residues(chain):
    aa_residues = [residue.id for residue in chain if residue.resname in AAA_DICT]
    return aa_residues


def get_interacting_residues(model, r_cutoff=5, skip_hetatm_chains=True):
    """Return residue-residue interactions between all chains in `model`.

    Parameters
    ----------
    model : biopython.Model
        Model to analyse.

    Returns
    -------
    dict
        A dictionary of interactions between chains i (0..n-1) and j (i+1..n).
        Keys are (chain_idx, chain_id, residue_idx, residue_resnum, residue_amino_acid) tuples.
        (e.g. (0, 'A', 0, '0', 'M'), (0, 1, '2', 'K'), ...)
        Values are a list of tuples having the same format as the keys.

    Examples
    --------
    You can reverse the order of keys and values like this::

        complement = dict()
        for key, values in get_interacting_chains(model):
            for value in values:
                complement.setdefault(value, set()).add(key)


    You can get a list of all interacting chains using this command::

        {(key[0], value[0])
         for (key, values) in get_interacting_chains(model).items()
         for value in values}

    """
    interactions_between_chains = dict()

    # Chain 1
    for chain_1_idx, chain_1 in enumerate(model):
        if skip_hetatm_chains and chain_is_hetatm(chain_1):
            message = "Skipping chain_1 with idx {} because it contains only hetatms.".format(
                chain_1_idx
            )
            logger.debug(message)
            continue
        chain_1_residue_ids = get_aa_residues(chain_1)

        # Chain 2
        for j, chain_2 in enumerate(model.child_list[chain_1_idx + 1 :]):
            chain_2_idx = chain_1_idx + 1 + j
            if skip_hetatm_chains and chain_is_hetatm(chain_2):
                message = "Skipping chain_2 with idx {} because it contains only hetatms.".format(
                    chain_2_idx
                )
                logger.debug(message)
                continue
            chain_2_residue_ids = get_aa_residues(chain_2)
            ns = NeighborSearch(list(chain_2.get_atoms()))

            # Residue 1
            for residue_1 in chain_1:
                try:
                    residue_1_idx = chain_1_residue_ids.index(residue_1.id)
                except ValueError:
                    continue
                residue_1_resnum = str(residue_1.id[1]) + residue_1.id[2].strip()
                residue_1_aa = convert_aa(residue_1.resname, quiet=True)
                residue_1_key = (
                    chain_1_idx,
                    chain_1.id,
                    residue_1_idx,
                    residue_1_resnum,
                    residue_1_aa,
                )
                interacting_residues = set()
                for atom_1 in residue_1:
                    interacting_residues.update(ns.search(atom_1.get_coord(), r_cutoff, "R"))

                # Residue 2
                interacting_residue_ids = []
                for residue_2 in interacting_residues:
                    try:
                        residue_2_idx = chain_2_residue_ids.index(residue_2.id)
                    except ValueError:
                        continue
                    residue_2_resnum = str(residue_2.id[1]) + residue_2.id[2].strip()
                    residue_2_aa = convert_aa(residue_2.get_resname(), quiet=True)
                    residue_2_key = (
                        chain_2_idx,
                        chain_2.id,
                        residue_2_idx,
                        residue_2_resnum,
                        residue_2_aa,
                    )
                    interacting_residue_ids.append(residue_2_key)
                if interacting_residue_ids:
                    interactions_between_chains.setdefault(residue_1_key, set()).update(
                        interacting_residue_ids
                    )

    return interactions_between_chains


class SelectChains(Select):
    """Select class which accepts only the specified chains."""

    def __init__(self, chain_letters, ns_chain_letters=None, ns=None, r_cutoff=None):
        self.chain_letters = chain_letters
        self.ns_chain_letters = ns_chain_letters
        self.ns = ns
        self.r_cutoff = r_cutoff

    def accept_residue(self, residue):
        chain_id = residue.parent.id
        if chain_id in self.chain_letters:
            return True
        elif (self.ns_chain_letters and self.ns) and (chain_id in self.ns_chain_letters):
            for atom in residue:
                if self.ns.search(atom.get_coord(), self.r_cutoff, "C"):
                    return True
        return False


class StructureParser:
    """
    """

    METHYLATED_LYSINES = ["MLZ", "MLY", "M3L"]
    LYSINE_ATOMS = ["N", "CA", "CB", "CG", "CD", "CE", "NZ", "C", "O"]

    def __init__(self, pdb_file, chain_ids=None, domain_defs=None):
        """.

        Parameters
        ----------
        pdb_file : str
            Full path and filename of the structure.
        output_dir : str
            Folder where to save extracted structures and sequences.
        chain_ids : list
            Chains of the structure that should be kept.
        """
        self.pdb_id = get_pdb_id(pdb_file)
        self.pdb_file = pdb_file
        self.input_structure = get_pdb_structure(self.pdb_file, self.pdb_id)

        if chain_ids is None:
            self.chain_ids = [chain.id for chain in self.input_structure[0].child_list]
        elif isinstance(chain_ids, str):
            self.chain_ids = chain_ids.split(",")
        elif isinstance(chain_ids, list) or isinstance(chain_ids, tuple):
            self.chain_ids = list(chain_ids)
        else:
            raise Exception
        # Rename empty chains (i.e. chain.id == ' ')
        model = structure[0]
        chain_ids = {chain.id for chain in model.child_list}
        for chain in model.child_list:
            if chain.id in [" ", "Z"]:
                chain_ids.remove(chain.id)
                chain.id = next(c for c in string.ascii_uppercase if c not in chain_ids)
                chain_ids.add(chain.id)
        model.child_dict = {chain.id: chain for chain in model.child_list}

        self.r_cutoff = 6  # remove hetatms more than x A away from the main chain(s)

        self.domain_boundaries = []
        for domain_def in domain_defs:
            self.domain_boundaries.append(
                decode_domain_def(domain_def, merge=False, return_string=True)
            )

        self.unique_id = "pdb_id: {}, chain_ids: {}".format(self.pdb_id, self.chain_ids)

    def __call__(self):
        # Contact distance
        contact_distance = None
        if chain_id_other is not None:
            try:
                contact_distance = self.get_interchain_distances(chain_id, mutation)
                contact_distance = contact_distance[chain_id][chain_id_other]
                logger.debug(
                    "The shortest interchain distance between chain {} and chain {} is {}".format(
                        chain_id, chain_id_other, contact_distance
                    )
                )
                if not contact_distance:
                    raise ValueError()
            except (IndexError, KeyError, ValueError) as e:
                logger.warning(
                    "Could not calculate the shortest contact distance between two chains!"
                )
                logger.warning(e)
                logger.warning(contact_distance)
                raise

    def extract(self):
        """Extract the wanted chains out of the PDB file.

        Remove water atoms and selects the domain regions (i.e. selects only those parts
        of the domain that are within the domain boundaries specified).
        """
        logger.debug("Extracting {}...".format(self.unique_id))

        model = self.input_structure[0]  # assuming that model 0 is always the desired one
        new_structure = Bio.PDB.Structure.Structure(self.pdb_id)
        new_model = Bio.PDB.Model.Model(0)

        # Always assigning hetatms to chain 'Z' may lead to undesirable performance
        # when the PDB stucture actually has a chain 'Z'.
        # As of 2015, there are ~1300 structures with chain 'Z' in the elaspic.domain table.
        # TODO: Convert `pdb_chain` tables in the database to use binary collation.
        # I think the Bio.PDB module may have to be upgraded too as it currently does not support
        # lowercase chain ids.
        hetatm_chain_id = "Z"
        hetatm_chain = Bio.PDB.Chain.Chain(hetatm_chain_id)

        # Loop over every chain and every residue and make sure that everything is ok
        chain_idx = 0
        while chain_idx < len(self.chain_ids):
            chain_id = self.chain_ids[chain_idx]
            chain = model[chain_id]

            (
                chain_numbering,
                domain_start_idxs,
                domain_end_idxs,
            ) = self._get_domain_def_idxs_for_chain(chain, chain_idx)
            logger.debug(
                "domain_def: %s, domain_start_idxs: %s, domain_end_idxs: %s",
                self.domain_boundaries,
                domain_start_idxs,
                domain_end_idxs,
            )

            res_idx = 0
            while res_idx < len(chain):
                res = chain.child_list[res_idx]
                original_res_id = res.id

                # Move water to the hetatm chain
                if res.id[0] == "W":
                    self._move_hetatm_to_hetatm_chain(chain, hetatm_chain, res, echo=False)
                    continue

                #                # Move heteroatoms to the hetatm chain
                #                if res.id[0] != ' ':
                #                    self._move_hetatm_to_hetatm_chain(chain, hetatm_chain, res, echo=True)
                #                    continue

                # Now treating all unusual amino acids as hetatms
                # Convert methylated lysines to regular lysines
                if res.resname in METHYLATED_LYSINES:
                    self._correct_methylated_lysines(res)

                # Move hetatms to the hetatm chain
                if res.resname not in AAA_DICT:
                    self._move_hetatm_to_hetatm_chain(chain, hetatm_chain, res)
                    continue

                # Cut each chain to domain boundaries
                residue_is_outside_domain = self._residue_outside_domain(
                    chain, chain_numbering, domain_start_idxs, domain_end_idxs, res
                )
                if residue_is_outside_domain:
                    chain.detach_child(original_res_id)
                    continue

                res_idx += 1

            if len(chain):
                new_model.add(chain)
                chain_idx += 1
            else:
                logger.debug("Chain {} is empty! Removing...".format(chain.id))
                self.chain_ids.remove(chain.id)

        # Make sure that the new model is not empty
        if not list(new_model.get_atoms()):
            raise structure_tools.exc.PDBEmptySequenceError(self.unique_id)

        # Remove hetatms if they are > 6A away from the chains of interest.
        self._remove_distant_hatatms(new_model, hetatm_chain)

        if hetatm_chain:
            logger.debug("Adding hetatm chain of length {}".format(len(hetatm_chain)))
            new_model.add(hetatm_chain)
            self.hetatm_chain_id = hetatm_chain_id
        else:
            self.hetatm_chain_id = None

        # If the hetatm chain is not empty, add it to the model
        new_structure.add(new_model)
        self.structure = new_structure
        logger.debug("PDB {} extracted successfully.".format(self.pdb_id))

        self.interactions_between_chains = get_interacting_residues(
            self.structure[0], self.r_cutoff, True
        )

        self.interacting_chain_ids = {
            (key[1], value[1])
            for (key, values) in self.interactions_between_chains.items()
            for value in values
        }
        self.interacting_chain_idxs = {
            (key[0], value[0])
            for (key, values) in self.interactions_between_chains.items()
            for value in values
        }

    def get_structure_file(self, chains, ext=".pdb"):
        return op.join(self.working_dir, self.sp.pdb_id + chains + ext)

    def _validate_mutation(self, resname, mutation):
        valid_aa = [mutation[0].upper(), mutation[-1].upper()]
        if AAA_DICT.get(resname, resname) not in valid_aa:
            logger.warning(resname)
            logger.warning(mutation)
            logger.warning(AAA_DICT[resname])
            raise structure_tools.exc.MutationMismatchError()

    def get_chain_sequence_and_numbering(self, chain_id, *args, **varargs):
        """Call ``get_chain_sequence_and_numbering`` using chain with id ``chain_id``."""
        chain = self.structure[0][chain_id]
        return get_chain_sequence_and_numbering(chain, *args, **varargs)

    def get_chain_seqres_sequence(self, chain_id, *args, **varargs):
        """Call ``get_chain_seqres_sequence`` using chain with id ``chain_id``."""
        chain = self.structure[0][chain_id]
        return get_chain_seqres_sequence(chain, *args, **varargs)

    def save_structure(self, output_dir="", remove_disordered=False):
        if remove_disordered:
            self._unset_disordered_flags()

        io = PDBIO()
        io.set_structure(self.structure)

        try:
            # Save all chains together
            outFile = op.join(output_dir, self.pdb_id + "".join(self.chain_ids) + ".pdb")
            io.save(outFile)
            if len(self.chain_ids) > 1:
                # Save each chain individually
                for chain_id in self.chain_ids:
                    chain = self.structure[0][chain_id]
                    if chain_is_hetatm(chain):
                        continue
                    outFile = op.join(output_dir, self.pdb_id + chain_id + ".pdb")
                    atom_list = [atom for atom in self.structure[0][chain_id].get_atoms()]
                    hetatm_chain_ns = NeighborSearch(atom_list)
                    select = SelectChains(
                        chain_id, self.hetatm_chain_id, hetatm_chain_ns, self.r_cutoff
                    )
                    io.save(outFile, select=select)
            if len(self.chain_ids) > 2:
                # Save each interacting chain pair.
                for chain_ids in self.interacting_chain_ids:
                    outFile = op.join(output_dir, self.pdb_id + "".join(chain_ids) + ".pdb")
                    atom_list = [
                        atom
                        for atom in self.structure[0][chain_id].get_atoms()
                        for chain_id in chain_ids
                    ]
                    hetatm_chain_ns = NeighborSearch(atom_list)
                    select = SelectChains(
                        chain_ids, self.hetatm_chain_id, hetatm_chain_ns, self.r_cutoff
                    )
                    io.save(outFile, select=select)

        except AttributeError as e:
            if remove_disordered:
                raise (e)
            self.save_structure(output_dir=output_dir, remove_disordered=True)

    def save_sequences(self, output_dir=""):
        self.chain_numbering_extended_dict = {}
        self.chain_sequence_dict = {}
        for chain_id in self.chain_ids:
            chain_sequence, chain_numbering_extended = self.get_chain_sequence_and_numbering(
                chain_id
            )
            self.chain_numbering_extended_dict[chain_id] = chain_numbering_extended
            self.chain_sequence_dict[chain_id] = chain_sequence
            with open(op.join(output_dir, self.pdb_id + chain_id + ".fasta"), "w") as f:
                f.write(">" + self.pdb_id + chain_id + "\n")
                f.write(chain_sequence + "\n")
                f.write("\n")

    def _get_domain_def_idxs_for_chain(self, chain, chain_idx):
        if not self.domain_boundaries or not self.domain_boundaries[chain_idx]:
            return None, None, None

        __, chain_numbering = get_chain_sequence_and_numbering(chain)
        try:
            domain_start_idxs, domain_end_idxs = [
                tuple(chain_numbering.index(resid) for resid in resids)
                for resids in zip(*self.domain_boundaries[chain_idx])
            ]
        except Exception as e:
            print(str(e))
            raise structure_tools.exc.PDBDomainDefsError(self.unique_id)

        return chain_numbering, domain_start_idxs, domain_end_idxs

    def _correct_methylated_lysines(self, res):
        new_resname = "LYS"
        new_resid = (" ", res.id[1], res.id[2])
        logger.debug(
            "Renaming residue {} {} to {} {}".format(res.resname, res.id, new_resname, new_resid)
        )
        res.resname = new_resname
        res.id = new_resid
        atom_idx = 0
        while atom_idx < len(res):
            atom_id = res.child_list[atom_idx].id
            if atom_id not in LYSINE_ATOMS:
                logger.debug(
                    "Removing atom {} from residue {} {}.".format(atom_id, res.resname, res.id)
                )
                res.detach_child(atom_id)
            else:
                atom_idx += 1

    def _move_hetatm_to_hetatm_chain(self, chain, hetatm_chain, res, echo=False):
        # logger.debug(
        #     'Moving hetatm residue {} {} to the hetatm chain'
        #     .format(res.resname, res.id))
        chain.detach_child(res.id)
        hetatm_res = res
        hetatm_res.id = (hetatm_res.id[0], len(hetatm_chain) + 1, hetatm_res.id[2])
        hetatm_chain.add(hetatm_res)

    def _residue_outside_domain(
        self, chain, chain_numbering, domain_start_idxs, domain_end_idxs, res
    ):
        """Return `True` if residue ``res`` is outside the domain."""
        if domain_start_idxs is None or domain_end_idxs is None:
            return False

        resid = str(res.id[1]) + res.id[2].strip()
        resid_idx = chain_numbering.index(resid)
        for domain_start_idx, domain_end_idx in zip(domain_start_idxs, domain_end_idxs):
            if resid_idx >= domain_start_idx and resid_idx <= domain_end_idx:
                # Residue is inside the domain
                return False

        # Residue is outside the domain
        return True

    def _remove_distant_hatatms(self, new_model, hetatm_chain):
        """Detach hetatms that are more than ``self.r_cutoff`` away from the main chain(s)."""
        ns = NeighborSearch(list(new_model.get_atoms()))
        hetatm_chain.id = [c for c in reversed(string.ascii_uppercase) if c not in self.chain_ids][
            0
        ]
        res_idx = 0
        while res_idx < len(hetatm_chain):
            res_1 = hetatm_chain.child_list[res_idx]
            in_contact = False
            for atom_1 in res_1:
                interacting_residues = ns.search(atom_1.get_coord(), self.r_cutoff, "R")
                if interacting_residues:
                    # logger.debug(res_1.id)
                    # logger.debug(interacting_residues)
                    in_contact = True
            if in_contact:
                res_idx += 1
                continue
            # logger.debug('Detaching child: {}'.format(res_1.id))
            hetatm_chain.detach_child(res_1.id)

    # def _unset_disordered_flags(self):
    #     """Change atom and residue ``disordered`` flag to `False`.
    #
    #     Otherwise, Biopython may crash when saving the PDB structure.
    #     """
    #     logger.debug('Setting all residues and atoms marked as disorded to non-disorded')
    #     for m in self.structure:
    #         for c in m:
    #             for r in c:
    #                 if r.is_disordered() or r.disordered:
    #                     logger.debug(
    #                         'Changing disordered_flag on residue {} from {} to 0'
    #                         .format(r, r.disordered))
    #                     r.disordered = 0
    #                 for a in r:
    #                     if a.is_disordered() or a.disordered_flag:
    #                         logger.debug(
    #                             'Changing disordered_flag on atom {} from {} to 0'
    #                             .format(a, a.disordered_flag))
    #                         a.disordered_flag = 0
    #
    def get_interchain_distances(self, pdb_chain=None, pdb_mutation=None, cutoff=None):
        """Calculate distance between two chains."""
        model = self.sp.structure[0]
        shortest_interchain_distances = {}
        # Chain 1
        for i, chain_1_id in enumerate(self.sp.chain_ids):
            shortest_interchain_distances[chain_1_id] = {}
            chain_1 = model[chain_1_id]
            if pdb_chain:
                if chain_1_id == pdb_chain:
                    continue
                chain_2_ids = [pdb_chain]
            else:
                chain_2_ids = self.sp.chain_ids[i + 1 :]
            # Chain 2
            for chain_2_id in chain_2_ids:
                chain_2 = model[chain_2_id]
                min_r = cutoff
                # Residue 1
                for residue_1 in chain_1:
                    if residue_1.resname not in AAA_DICT or residue_1.id[0] != " ":
                        continue
                    # Residue 2
                    for residue_2_idx, residue_2 in enumerate(chain_2):
                        if residue_2.resname not in AAA_DICT or residue_2.id[0] != " ":
                            continue
                        if pdb_mutation:
                            if str(residue_2.id[1]) != pdb_mutation[1:-1]:
                                continue
                            if (
                                structure_tools.convert_aa(residue_2.resname) != pdb_mutation[0]
                            ) and (
                                structure_tools.convert_aa(residue_2.resname) != pdb_mutation[-1]
                            ):
                                logger.debug(pdb_mutation)
                                logger.debug(structure_tools.convert_aa(residue_2.resname))
                                logger.debug(residue_2.id)
                                raise structure_tools.exc.MutationMismatchError()
                        # Atom 1
                        for atom_1 in residue_1:
                            # Atom 2
                            for atom_2 in residue_2:
                                r = structure_tools.calculate_distance(atom_1, atom_2, min_r)
                                if min_r is None or (r is not None and r < min_r):
                                    min_r = r

                shortest_interchain_distances[chain_1_id][chain_2_id] = min_r

        if not shortest_interchain_distances:
            logger.warning(
                "get_interchain_distances({pdb_chain}, {pdb_mutation}, {cutoff}) failed!".format(
                    pdb_chain=pdb_chain, pdb_mutation=pdb_mutation, cutoff=cutoff
                )
            )
            raise Exception()

        _shortest_interchain_distances_complement = {}
        for key in shortest_interchain_distances:
            for key_2, value in shortest_interchain_distances[key].items():
                _shortest_interchain_distances_complement.setdefault(key_2, dict())[key] = value
        shortest_interchain_distances.update(_shortest_interchain_distances_complement)

        all_chains = {key for key in shortest_interchain_distances}
        all_chains.update(
            {
                key_2
                for key in shortest_interchain_distances
                for key_2 in shortest_interchain_distances[key]
            }
        )

        if set(all_chains) != set(self.sp.chain_ids):
            logger.warning(
                "get_interchain_distances({pdb_chain}, {pdb_mutation}, {cutoff}) failed!".format(
                    pdb_chain=pdb_chain, pdb_mutation=pdb_mutation, cutoff=cutoff
                )
            )
            logger.warning("Did not calculate chain distances for all chain pairs!")
            logger.warning("all_chains: {}".format(all_chains))
            logger.warning("self.sp.chain_ids: {}".format(self.sp.chain_ids))
            raise Exception()

        for key_1, value_1 in shortest_interchain_distances.items():
            logger.debug(
                "Calculated interchain distances between chain {} and chains {}.".format(
                    key_1, ", ".join(list(value_1.keys()))
                )
            )

        return shortest_interchain_distances
