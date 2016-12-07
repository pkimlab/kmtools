import logging
from copy import copy

import Bio.PDB
import pandas as pd

import kmtools.structure_tools.monkeypatch_biopython  # noqa
from kmtools.structure_tools import AAA_DICT, AMINO_ACIDS, COMMON_HETATMS, CHAIN_IDS

logger = logging.getLogger(__name__)

#: Source: http://www.wwpdb.org/documentation/procedure#toc_4

R_CUTOFF = 5


class NotDisordered(Bio.PDB.Select):
    """Select only non-disordered residues and set their altloc flag to ' '.

    Source: http://biopython.org/wiki/Remove_PDB_disordered_atoms
    """
    def accept_residue(self, residue):
        if not residue.disordered:
            return True
        elif any(self.accept_atom(atom) for atom in residue):
            residue.disordered = False
            return True
        else:
            logger.debug("Ignoring residue %s.", residue)
            return False

    def accept_atom(self, atom):
        if not atom.disordered_flag:
            return True
        elif atom.altloc == 'A':
            atom.disordered_flag = False
            atom.altloc = ' '
            return True
        else:
            logger.debug("Ignoring atom %s.", atom)
            return False


def save_structure(structure, filename, include_hetatms=True):
    io = Bio.PDB.PDBIO()
    select = None if include_hetatms else NotDisordered()
    io.set_structure(structure)
    io.save(filename, select=select)


def _copy_residue(residue, residue_idx):
    new_residue = Bio.PDB.Residue.Residue(
        id=residue_idx,
        resname=residue.resname,
        segid=' ')
    # Note that this does not change the `.parent` attribute on the children.
    new_residue.child_list = copy(residue.child_list)
    new_residue.child_dict = copy(residue.child_dict)
    return new_residue


def _copy_structure(structure):
    raise NotImplementedError


def process_structure(structure):
    mapping = {}
    # Create new structure
    new_structure = Bio.PDB.Structure.Structure(structure.id)
    # Create model for storing common HETATMs
    hetatm_model = Bio.PDB.Model.Model(len(structure.child_list))
    hetatm_chain = Bio.PDB.Chain.Chain(CHAIN_IDS[0])
    hetatm_model.add(hetatm_chain)
    hetatm_residue_idx = 0
    # ---
    model_idx = 0
    for model in structure:
        # Create new model
        new_model = Bio.PDB.Model.Model(model_idx)
        chain_idx = 0
        for chain in model:
            # Create new chain
            new_chain = Bio.PDB.Chain.Chain(CHAIN_IDS[chain_idx])
            residue_idx = 0
            for residue in chain:
                if residue.resname in AMINO_ACIDS:
                    new_residue = _copy_residue(residue, (' ', residue_idx, ' '))
                    new_chain.add(new_residue)
                    mapping[(model.id, chain.id, residue.id)] = (
                        new_model.id, new_chain.id, new_residue.id)
                    residue_idx += 1
                else:
                    new_residue = _copy_residue(residue, (' ', hetatm_residue_idx, ' '))
                    hetatm_chain.add(new_residue)
                    mapping[(model.id, chain.id, residue.id)] = (
                        hetatm_model.id, hetatm_chain.id, new_residue.id)
                    hetatm_residue_idx += 1
            # Finallize chain
            new_model.add(new_chain)
            chain_idx += 1
        # Finallize model
        new_structure.add(new_model)
        model_idx += 1
    # Finallize HETATM model
    new_structure.add(hetatm_model)
    # Add residue mapping
    new_structure.xtra['mapping'] = mapping
    return new_structure


def _iter_interchain_ns(structure, interchain=True):
    from Bio.PDB import NeighborSearch
    for model_1_idx, model_1 in enumerate(structure):
        for chain_1_idx, chain_1 in enumerate(model_1):
            ns = NeighborSearch(list(chain_1.get_atoms()))
            for model_2_idx, model_2 in enumerate(structure):
                if model_1_idx > model_2_idx:
                    continue
                for chain_2_idx, chain_2 in enumerate(model_2):
                    if (interchain and
                            (model_1_idx == model_2_idx and chain_1_idx > chain_2_idx)):
                        continue
                    if (not interchain and
                            (model_1_idx != model_2_idx or chain_1_idx != chain_2_idx)):
                        continue
                    for residue_2_idx, residue_2 in enumerate(chain_2):
                        yield model_1.id, chain_1.id, model_2.id, chain_2.id, ns, residue_2


def get_interactions(structure, interchain=True):
    """Return residue-residue interactions within and between chains in `structure`.

    .. todo::

        This could probably be sped up by using the
        :py:meth:`Bio.PDB.NeighborSearch.search_all` method.

    Preconditions
    -------------
    - For each chain, `residue.id[1]` is unique.

    Parameters
    ----------
    structure : Bio.PDB.Structure.Structure
        Structure to analyse.
    interchain : :cls:`bool`
        Whether to include interactions between chains.
    """
    columns = [
        'model_1_id', 'chain_1_id', 'model_2_id', 'chain_2_id',
        'residue_1_id', 'residue_1_name', 'residue_1_aa',
        'residue_2_id', 'residue_2_name', 'residue_2_aa',
    ]
    results = []
    for model_1_id, chain_1_id, model_2_id, chain_2_id, ns, residue_2 in (
        _iter_interchain_ns(structure, interchain=interchain)
    ):
        if residue_2.resname in COMMON_HETATMS:
            continue
        seen = set()
        interacting_residues = [
            r
            for a in residue_2
            for r in ns.search(a.get_coord(), R_CUTOFF, 'R')
            if r.id not in seen and not seen.add(r.id)
        ]
        for residue_1 in interacting_residues:
            if residue_1.resname in COMMON_HETATMS:
                continue
            row = (
                model_1_id, chain_1_id, model_2_id, chain_2_id,
                residue_1.id[1], residue_1.resname, AAA_DICT.get(residue_1.resname),
                residue_2.id[1], residue_2.resname, AAA_DICT.get(residue_2.resname),
            )
            results.append(row)
    # df
    columns_to_sort = [
        'model_1_id', 'model_2_id', 'chain_1_id', 'chain_2_id', 'residue_1_id', 'residue_2_id'
    ]
    df = pd.DataFrame(results, columns=columns).sort_values(columns_to_sort)
    df.index = range(len(df))
    return df


def get_chain_sequence(chain):
    return ''.join(AAA_DICT[r.resname] for r in chain)
