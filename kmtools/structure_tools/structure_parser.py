import logging

import numpy as np

from kmbio.PDB import Chain, Model, NeighborSearch, Residue, Structure
from kmtools.structure_tools import AAA_DICT, CHAIN_IDS, LYSINE_ATOMS, METHYLATED_LYSINES

logger = logging.getLogger(__name__)

# Source: http://www.wwpdb.org/documentation/procedure#toc_4
R_CUTOFF = 4


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
    if ((type(atom_1) == type(atom_2) == list) or (type(atom_1) == type(atom_2) == tuple)):
        a = atom_1
        b = atom_2
    elif hasattr(atom_1, 'coord') and hasattr(atom_2, 'coord'):
        a = atom_1.coord
        b = atom_2.coord
    else:
        raise Exception('Unsupported format {} {}'.format(type(atom_1), type(atom_2)))

    assert (len(a) == 3 and len(b) == 3)
    if cutoff is None or all(abs(p - q) <= cutoff for p, q in zip(a, b)):
        return euclidean_distance(a, b)


def process_structure(structure: Structure) -> Structure:
    """Process structure inplace.

    Note
    ----
    Running this function should not be neccessary if you obtain your structure from an mmCIF file
    using ``use_auth_id=False``.

    Warnings
    --------
    - **Very** weird things happen if this function does not
      make a copy of the structure first.
    """
    structure = structure.copy()
    models = []
    model_idx = 0
    residue_mapping_fw = {}
    # Create model for storing common HETATMs
    hetatm_model = Model(len(structure))
    hetatm_chain = Chain(CHAIN_IDS[0])
    hetatm_model.add(hetatm_chain)
    hetatm_residue_idx = 0
    for model in structure:
        # Create new model
        new_model = Model(model_idx)
        chain_idx = 0
        for chain in model:
            # Create new chain
            new_chain = Chain(CHAIN_IDS[chain_idx])
            residue_idx = 0
            for residue in chain:
                if residue.resname in METHYLATED_LYSINES:
                    residue = _correct_methylated_lysines(residue)
                if residue.resname in AAA_DICT:
                    new_residue = _copy_residue(residue, (' ', residue_idx, ' '))
                    new_chain.add(new_residue)
                    residue_mapping_fw[(model.id, chain.id,
                                        residue.id)] = (new_model.id, new_chain.id, new_residue.id)
                    residue_idx += 1
                else:
                    new_residue = _copy_residue(residue, (' ', hetatm_residue_idx, ' '))
                    hetatm_chain.add(new_residue)
                    residue_mapping_fw[(model.id, chain.id,
                                        residue.id)] = (hetatm_model.id, hetatm_chain.id,
                                                        new_residue.id)
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
    structure.xtra['residue_mapping_bw'] = {v: k for k, v in residue_mapping_fw.items()}
    assert len(structure.xtra['residue_mapping_fw']) == \
        len(structure.xtra['residue_mapping_bw'])
    assert all(' ' not in r.resname for r in structure.residues)
    return structure


def _copy_residue(residue, residue_idx):
    new_residue = Residue(
        id=residue_idx,
        resname=residue.resname.strip(' '),
        segid=residue.segid,  # ' '
        children=list(residue.atoms))
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
        atom_id = list(residue.atoms)[atom_idx].id
        if atom_id not in LYSINE_ATOMS:
            logger.debug('Removing atom %s from residue %s %s.', atom_id, residue.resname,
                         residue.id)
            del residue[atom_id]
        else:
            atom_idx += 1
    return residue


def extract(structure, mcd=None, select_hetatms=None):
    """Update `structure` to contain only the chains and residues of interest."""
    if mcd is None:
        mcd = [(m.id, c.id, range(0, len(c))) for m in structure for c in m]
        select_hetatms = False  # will be copied as-is
    else:
        select_hetatms = True

    # Initialize
    new_model = Model(0)
    residue_mapping_fw = dict()

    # Regular chains
    for i, (model_id, chain_id, domain_def) in enumerate(mcd):
        new_chain = Chain(CHAIN_IDS[i], children=structure[model_id][chain_id].ix[domain_def])
        new_model.add(new_chain)
        residue_mapping_fw.update({
            k: (new_model.id, new_chain.id, r_id)
            for k, (m_id, c_id, r_id) in structure.xtra['residue_mapping_fw'].items()
            if m_id == model_id and c_id == chain_id and r_id[1] in domain_def
        })

    # Hetatm chain
    if select_hetatms:
        assert len(structure) > 1 and len(structure.ix[-1]) == 1
        hetatm_model = structure.ix[-1]
        hetatm_chain = hetatm_model.ix[-1]
        new_hetatm_chain = copy_hetatm_chain(
            Structure(id=structure.id, children=new_model), hetatm_chain, CHAIN_IDS[i + 1])
        new_model.add(new_hetatm_chain)
        hetatm_residue_ids = set(r.id for r in new_hetatm_chain)
        residue_mapping_fw.update({
            k: (new_model.id, new_chain.id, r_id)
            for k, (m_id, c_id, r_id) in structure.xtra['residue_mapping_fw'].items()
            if m_id == hetatm_model.id and c_id == hetatm_chain.id and r_id in hetatm_residue_ids
        })

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
    hetatm_atoms = list(hetatm_chain.atoms)
    ns = NeighborSearch(hetatm_atoms)
    # New hetatm chain (keeping only proximal hetatms)
    new_hetatm_chain = Chain(new_hetatm_chain_id)
    new_hetatm_residues = [
        residue for residue in hetatm_residues
        if any(ns.search(a.coord, R_CUTOFF) for a in residue)
    ]
    for new_hetatm_residue in new_hetatm_residues:
        new_hetatm_chain.add(new_hetatm_residue)
    return new_hetatm_chain


def get_chain_sequence(chain):
    return ''.join(AAA_DICT[r.resname] for r in chain)
