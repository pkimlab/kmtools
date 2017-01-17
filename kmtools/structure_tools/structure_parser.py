import logging

import pandas as pd

import Bio.PDB
import kmtools.structure_tools.monkeypatch_biopython  # noqa
from kmtools import df_tools
from kmtools.structure_tools import (AAA_DICT, AMINO_ACIDS, CHAIN_IDS, COMMON_HETATMS,
                                     LYSINE_ATOMS, METHYLATED_LYSINES)

logger = logging.getLogger(__name__)

# Source: http://www.wwpdb.org/documentation/procedure#toc_4
R_CUTOFF = 5


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
                if residue.resname in METHYLATED_LYSINES:
                    residue = _correct_methylated_lysines(residue)
                if residue.resname in AMINO_ACIDS:
                    new_residue = _copy_residue(residue, (' ', residue_idx, ' '))
                    new_chain.add(new_residue)
                    mapping[(model.id, chain.id, residue.id[1])] = (
                        new_model.id, new_chain.id,
                        ''.join(str(x) for x in new_residue.id).strip())
                    residue_idx += 1
                else:
                    new_residue = _copy_residue(residue, (' ', hetatm_residue_idx, ' '))
                    hetatm_chain.add(new_residue)
                    mapping[(model.id, chain.id, residue.id[1])] = (
                        hetatm_model.id, hetatm_chain.id,
                        ''.join(str(x) for x in new_residue.id).strip())
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
    new_structure.xtra['reverse_mapping'] = {
        v: k for k, v in new_structure.xtra['mapping'].items()}
    assert all(' ' not in r.resname for r in new_structure.get_residues())
    return new_structure


def _copy_residue(residue, residue_idx):
    new_residue = residue.copy()
    new_residue.id = residue_idx
    new_residue.resname = residue.resname.strip(' ')
    new_residue.segin = ' '
    if residue.resname != new_residue.resname:
        logger.debug("Renamed residue %s to %s" % (residue.resname, new_residue.resname))
    return new_residue


def _correct_methylated_lysines(residue):
    """Remove methyl group from lysine residue."""
    logger.debug("Renaming residue %s, %s", residue.resname, residue.id)
    residue = _copy_residue(residue, (' ', residue.id[1], ' '))
    residue.resname = 'LYS'
    logger.debug("New name: %s %s", residue.resname, residue.id)
    atom_idx = 0
    while atom_idx < len(residue):
        atom_id = residue.child_list[atom_idx].id
        if atom_id not in LYSINE_ATOMS:
            logger.debug(
                'Removing atom %s from residue %s %s.', atom_id, residue.resname, residue.id)
            residue.detach_child(atom_id)
        else:
            atom_idx += 1
    return residue


def _iter_interchain_ns(structure, interchain=True):
    for model_1_idx, model_1 in enumerate(structure):
        for chain_1_idx, chain_1 in enumerate(model_1):
            atom_list = list(chain_1.get_atoms())
            if not atom_list:
                logger.debug("Skipping %s %s because it is empty.", model_1, chain_1)
                continue
            ns = Bio.PDB.NeighborSearch(atom_list)
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

    Parameters
    ----------
    structure : Bio.PDB.Structure.Structure
        Structure to analyse.
    interchain : :class:`bool`
        Whether to include interactions between chains.

    Notes
    -------------
    - For each chain, `residue.id[1]` is unique.

    """
    columns = [
        'structure_id', 'model_id_1', 'model_id_2', 'chain_id_1', 'chain_id_2',
        'residue_id_1', 'residue_id_2', 'residue_name_1', 'residue_name_2',
        'residue_aa_1', 'residue_aa_2',
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
                structure.id, model_1_id, model_2_id, chain_1_id, chain_2_id,
                residue_1.id[1], residue_2.id[1], residue_1.resname, residue_2.resname,
                AAA_DICT.get(residue_1.resname), AAA_DICT.get(residue_2.resname)
            )
            results.append(row)
    # df
    df = pd.DataFrame(results, columns=columns)
    df_reverse = df.rename(columns=df_tools.get_reverse_column).reindex(columns=pd.Index(columns))
    df = pd.concat([df, df_reverse]).drop_duplicates()
    df = df.sort_values(columns)  # will not work well for uppercase / lowercase columns
    df.index = range(len(df))
    return df


def extract(structure, chain_ids=None, domain_defs=None):
    """Select only the chains and residues of interest.

    .. todo::

        This has to be changed to work with structures that have several models
        (apart from the HETATM model).

    .. todo::

        Use a more unique chain id for HETATMS (e.g. `zz`).
    """
    model_id = 0
    # Validate input
    if chain_ids is None and domain_defs is None:
        logger.warning("Nothing to extract!")
        return structure
    if chain_ids is None:
        chain_ids = [chain.id for chain in structure[model_id]]
    elif isinstance(chain_ids, str):
        chain_ids = [chain_ids]
    elif not isinstance(chain_ids, list):
        raise ValueError
    if domain_defs is None:
        domain_defs = [
            slice(0, len(structure[model_id][chain_id]))
            for chain_id in chain_ids
        ]
    elif isinstance(domain_defs, slice):
        domain_defs = [slice]
    elif not isinstance(chain_ids, list):
        raise ValueError
    #
    new_structure = Bio.PDB.Structure.Structure(structure.id)
    new_model = Bio.PDB.Model.Model(model_id)
    new_structure.add(new_model)
    # Regular chains
    for chain_id, domain_def in zip(chain_ids, domain_defs):
        new_chain = Bio.PDB.Chain.Chain(chain_id)
        new_chain.add(structure[model_id][chain_id].child_list[domain_def])
        new_model.add(new_chain)
    # Hetatm chain
    assert len(structure.child_list[-1]) == 1
    hetatm_chain = structure.child_list[-1].child_list[-1]
    new_hetatm_chain = _get_hetatm_chain(new_structure, hetatm_chain)
    new_model.add(new_hetatm_chain)
    # Mapping
    new_structure_keys = {
        (model.id, chain.id, residue.id[1])
        for model in new_structure for chain in model for residue in chain
    }
    new_structure.xtra['mapping'] = {
        key: value
        for key, value in structure.items()
        if key in new_structure_keys
    }
    new_structure.xtra['reverse_mapping'] = {
        v: k for k, v in new_structure.xtra['mapping'].items()}
    return new_structure


def _get_hetatm_chain(structure, hetatm_chain):
    """Return `new_hetatm_chain` which contains only those residues from `hetatm_chain`
    which are less than `R_CUTOFF` away from a residue in `structure`.
    """
    # Old hetatm chain
    hetatm_residues = hetatm_chain.child_list
    hetatm_atoms = list(hetatm_chain.get_atoms())
    ns = Bio.PDB.NeighborSearch(hetatm_atoms)
    # New hetatm chain (keeping only proximal hetatms)
    new_hetatm_chain_id = 'Z'
    new_hetatm_chain = Bio.PDB.Chain.Chain(new_hetatm_chain_id)
    new_hetatm_residues = [
        residue
        for residue in hetatm_residues
        if any(ns.search(a.get_coord(), R_CUTOFF) for a in residue)
    ]
    new_hetatm_chain.add(new_hetatm_residues)
    return new_hetatm_chain


def get_chain_sequence(chain):
    return ''.join(AAA_DICT[r.resname] for r in chain)
