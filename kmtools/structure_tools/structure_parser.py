"""

Conventions
-----------

ID           | IDX                     | IDX starts with
-------------|-------------------------|-----------------
`model_id`   | `model_idx`             | 0
`chain_id`   | `CHAIN_IDS[chain_idx]`  | 0
`residue_id` | `('', residue_idx, '')` | 0

"""
import logging

import numpy as np
import pandas as pd

import kmbio.PDB
from kmtools import df_tools
from kmtools.structure_tools import (AAA_DICT, AMINO_ACIDS, CHAIN_IDS, COMMON_HETATMS,
                                     LYSINE_ATOMS, METHYLATED_LYSINES)

logger = logging.getLogger(__name__)

# Source: http://www.wwpdb.org/documentation/procedure#toc_4
R_CUTOFF = 5


def euclidean_distance(a, b):
    """Calculate the Euclidean distance between two lists or tuples of arbitrary length."""
    return np.sqrt(sum((a - b)**2 for a, b in zip(a, b)))


def calculate_distance(atom_1, atom_2, cutoff=None):
    """Calculate the distance between two points in 3D space.

    Parameters
    ----------
    cutoff : float, optional
        The maximum distance allowable between two points.
    """
    if ((type(atom_1) == type(atom_2) == list) or
            (type(atom_1) == type(atom_2) == tuple)):
        a = atom_1
        b = atom_2
    elif hasattr(atom_1, 'coord') and hasattr(atom_2, 'coord'):
        a = atom_1.coord
        b = atom_2.coord
    else:
        raise Exception('Unsupported format {} {}'.format(type(atom_1), type(atom_2)))

    assert(len(a) == 3 and len(b) == 3)
    if cutoff is None or all(abs(p - q) <= cutoff for p, q in zip(a, b)):
        return euclidean_distance(a, b)


def process_structure(structure):
    """Process structure inplace."""
    models = []
    model_idx = 0
    residue_mapping_fw = {}
    # Create model for storing common HETATMs
    hetatm_model = kmbio.PDB.Model.Model(len(structure))
    hetatm_chain = kmbio.PDB.Chain.Chain(CHAIN_IDS[0])
    hetatm_model.add(hetatm_chain)
    hetatm_residue_idx = 0
    for model in structure.values():
        # Create new model
        new_model = kmbio.PDB.Model.Model(model_idx)
        chain_idx = 0
        for chain in model.values():
            # Create new chain
            new_chain = kmbio.PDB.Chain.Chain(CHAIN_IDS[chain_idx])
            residue_idx = 0
            for residue in chain.values():
                if residue.resname in METHYLATED_LYSINES:
                    residue = _correct_methylated_lysines(residue)
                if residue.resname in AMINO_ACIDS:
                    new_residue = _copy_residue(residue, (' ', residue_idx, ' '))
                    new_chain.add(new_residue)
                    residue_mapping_fw[(model.id, chain.id, residue.id)] = (
                        new_model.id, new_chain.id, new_residue.id)
                    residue_idx += 1
                else:
                    new_residue = _copy_residue(residue, (' ', hetatm_residue_idx, ' '))
                    hetatm_chain.add(new_residue)
                    residue_mapping_fw[(model.id, chain.id, residue.id)] = (
                        hetatm_model.id, hetatm_chain.id, new_residue.id)
                    hetatm_residue_idx += 1
            # Finallize chain
            new_model.add(new_chain)
            chain_idx += 1
        # Finallize model
        models.append(new_model)
        model_idx += 1
    # Finallize HETATM model
    models.append(hetatm_model)
    # Add residue mapping
    structure.clear()
    structure.add(models)
    structure.xtra['residue_mapping_fw'] = residue_mapping_fw
    structure.xtra['residue_mapping_bw'] = {
        v: k for k, v in residue_mapping_fw.items()}
    assert len(structure.xtra['residue_mapping_fw']) == \
        len(structure.xtra['residue_mapping_bw'])
    assert all(' ' not in r.resname for r in structure.get_residues())


def _copy_residue(residue, residue_idx):
    new_residue = kmbio.PDB.Residue.Residue(
        id=residue_idx,
        resname=residue.resname.strip(' '),
        segid=residue.segid,  # ' '
        children=residue.values())
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
        atom_id = residue.ix[atom_idx].id
        if atom_id not in LYSINE_ATOMS:
            logger.debug(
                'Removing atom %s from residue %s %s.', atom_id, residue.resname, residue.id)
            del residue[atom_id]
        else:
            atom_idx += 1
    return residue


def _iter_interchain_ns(structure, interchain=True):
    for model_1_idx, model_1 in enumerate(structure.values()):
        for chain_1_idx, chain_1 in enumerate(model_1.values()):
            atom_list = list(chain_1.get_atoms())
            if not atom_list:
                logger.debug("Skipping %s %s because it is empty.", model_1, chain_1)
                continue
            ns = kmbio.PDB.NeighborSearch(atom_list)
            for model_2_idx, model_2 in enumerate(structure.values()):
                if model_1_idx > model_2_idx:
                    continue
                for chain_2_idx, chain_2 in enumerate(model_2.values()):
                    if (interchain and
                            (model_1_idx == model_2_idx and chain_1_idx > chain_2_idx)):
                        continue
                    if (not interchain and
                            (model_1_idx != model_2_idx or chain_1_idx != chain_2_idx)):
                        continue
                    for residue_2_idx, residue_2 in enumerate(chain_2.values()):
                        yield model_1.id, chain_1.id, model_2.id, chain_2.id, ns, residue_2


def get_interactions(structure, interchain=True):
    """Return residue-residue interactions within and between chains in `structure`.

    .. todo::

        This could probably be sped up by using the
        :py:meth:`kmbio.PDB.NeighborSearch.search_all` method.

    Parameters
    ----------
    structure : kmbio.PDB.Structure.Structure
        Structure to analyse.
    interchain : :class:`bool`
        Whether to include interactions between chains.

    Notes
    -----
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
            for a in residue_2.values()
            for r in ns.search(a.coord, R_CUTOFF, 'R')
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


def extract(structure, mcd=None, select_hetatms=None):
    """Update `structure` to contain only the chains and residues of interest."""
    if mcd is None:
        mcd = [(m.id, c.id, range(0, len(c))) for m in structure for c in m]
        select_hetatms = False  # will be copied as-is
    else:
        select_hetatms = True

    # Initialize
    new_model = kmbio.PDB.Model.Model(0)
    residue_mapping_fw = dict()

    # Regular chains
    for i, (model_id, chain_id, domain_def) in enumerate(mcd):
        new_chain = kmbio.PDB.Chain.Chain(
            CHAIN_IDS[i], children=structure[model_id][chain_id].ix[domain_def])
        new_model.add(new_chain)
        residue_mapping_fw.update({
            k: (new_model.id, new_chain.id, r_id)
            for k, (m_id, c_id, r_id) in structure.xtra['residue_mapping_fw'].items()
            if m_id == model_id and c_id == chain_id and r_id[1] in domain_def})

    # Hetatm chain
    if select_hetatms:
        assert len(structure) > 1 and len(structure.ix[-1]) == 1
        hetatm_model = structure.ix[-1]
        hetatm_chain = hetatm_model.ix[-1]
        new_hetatm_chain = copy_hetatm_chain(
            kmbio.PDB.Structure.Structure(id=structure.id, children=new_model),
            hetatm_chain, CHAIN_IDS[i + 1])
        new_model.add(new_hetatm_chain)
        hetatm_residue_ids = set(r.id for r in new_hetatm_chain)
        residue_mapping_fw.update({
            k: (new_model.id, new_chain.id, r_id)
            for k, (m_id, c_id, r_id) in structure.xtra['residue_mapping_fw'].items()
            if m_id == hetatm_model.id and c_id == hetatm_chain.id and r_id in hetatm_residue_ids})

    structure.clear()
    structure.add(new_model)
    structure.xtra['residue_mapping_fw'] = residue_mapping_fw
    structure.xtra['residue_mapping_bw'] = {v: k for (k, v) in residue_mapping_fw.items()}


def copy_hetatm_chain(structure, hetatm_chain, new_hetatm_chain_id):
    """Return `new_hetatm_chain` which contains only those residues from `hetatm_chain`
    which are less than `R_CUTOFF` away from some residue in `structure`.
    """
    # Old hetatm chain
    hetatm_residues = hetatm_chain.ix[:]
    hetatm_atoms = list(hetatm_chain.get_atoms())
    ns = kmbio.PDB.NeighborSearch(hetatm_atoms)
    # New hetatm chain (keeping only proximal hetatms)
    new_hetatm_chain = kmbio.PDB.Chain.Chain(new_hetatm_chain_id)
    new_hetatm_residues = [
        residue
        for residue in hetatm_residues
        if any(ns.search(a.coord, R_CUTOFF) for a in residue)
    ]
    for new_hetatm_residue in new_hetatm_residues:
        new_hetatm_chain.add(new_hetatm_residue)
    return new_hetatm_chain


def get_chain_sequence(chain):
    return ''.join(AAA_DICT[r.resname] for r in chain)
