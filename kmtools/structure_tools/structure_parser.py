import logging
from typing import Optional

import numpy as np
from kmbio.PDB import Chain, Model, NeighborSearch, Residue, Structure

from kmtools.structure_tools import (
    AAA_DICT,
    CHAIN_IDS,
    DNA_MAPPING_TO_CANONICAL,
    LYSINE_ATOMS,
    METHYLATED_LYSINES,
    RESIDUE_MAPPING_TO_CANONICAL,
    RNA_MAPPING_TO_CANONICAL,
)

logger = logging.getLogger(__name__)

# Source: http://www.wwpdb.org/documentation/procedure#toc_4
R_CUTOFF = 5


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


def process_structure(structure: Structure) -> Structure:
    """Process structure inplace.

    Notes:
        Running this function should not be neccessary if you obtain your structure
        from an mmCIF file using ``use_auth_id=False``.

    Warnings:
        **Very** weird things happen if this function does not make a copy of the structure first.
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
                    new_residue = _copy_residue(residue, (" ", residue_idx, " "))
                    new_chain.add(new_residue)
                    residue_mapping_fw[(model.id, chain.id, residue.id)] = (
                        new_model.id,
                        new_chain.id,
                        new_residue.id,
                    )
                    residue_idx += 1
                else:
                    new_residue = _copy_residue(residue, (" ", hetatm_residue_idx, " "))
                    hetatm_chain.add(new_residue)
                    residue_mapping_fw[(model.id, chain.id, residue.id)] = (
                        hetatm_model.id,
                        hetatm_chain.id,
                        new_residue.id,
                    )
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
    structure.xtra["residue_mapping_fw"] = residue_mapping_fw
    structure.xtra["residue_mapping_bw"] = {v: k for k, v in residue_mapping_fw.items()}
    assert len(structure.xtra["residue_mapping_fw"]) == len(structure.xtra["residue_mapping_bw"])
    assert all(" " not in r.resname for r in structure.residues)
    return structure


def _copy_residue(residue, residue_idx):
    new_residue = Residue(
        id=residue_idx,
        resname=residue.resname.strip(" "),
        segid=residue.segid,  # ' '
        children=list(residue.atoms),
    )
    if residue.resname != new_residue.resname:
        logger.debug("Renamed residue %s to %s" % (residue.resname, new_residue.resname))
    return new_residue


def _correct_methylated_lysines(residue):
    """Remove methyl group from lysine residue."""
    logger.debug("Renaming residue %s, %s", residue.resname, residue.id)
    residue = _copy_residue(residue, (" ", residue.id[1], " "))
    residue.resname = "LYS"
    logger.debug("New name: %s %s", residue.resname, residue.id)
    atom_idx = 0
    while atom_idx < len(residue):
        atom_id = list(residue.atoms)[atom_idx].id
        if atom_id not in LYSINE_ATOMS:
            logger.debug(
                "Removing atom %s from residue %s %s.", atom_id, residue.resname, residue.id
            )
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
        residue_mapping_fw.update(
            {
                k: (new_model.id, new_chain.id, r_id)
                for k, (m_id, c_id, r_id) in structure.xtra["residue_mapping_fw"].items()
                if m_id == model_id and c_id == chain_id and r_id[1] in domain_def
            }
        )

    # Hetatm chain
    if select_hetatms:
        assert len(structure) > 1 and len(structure.ix[-1]) == 1
        hetatm_model = structure.ix[-1]
        hetatm_chain = hetatm_model.ix[-1]
        new_hetatm_chain = copy_hetatm_chain(
            Structure(id=structure.id, children=new_model), hetatm_chain
        )
        new_hetatm_chain.id = CHAIN_IDS[i + 1]
        new_model.add(new_hetatm_chain)
        hetatm_residue_ids = set(r.id for r in new_hetatm_chain)
        residue_mapping_fw.update(
            {
                k: (new_model.id, new_chain.id, r_id)
                for k, (m_id, c_id, r_id) in structure.xtra["residue_mapping_fw"].items()
                if m_id == hetatm_model.id
                and c_id == hetatm_chain.id
                and r_id in hetatm_residue_ids
            }
        )

    structure.clear()
    structure.add(new_model)
    structure.xtra["residue_mapping_fw"] = residue_mapping_fw
    structure.xtra["residue_mapping_bw"] = {v: k for (k, v) in residue_mapping_fw.items()}


def copy_hetatm_chain(
    structure: Structure, hetatm_chain: Chain, r_cutoff: float = R_CUTOFF
) -> Chain:
    """Return `new_hetatm_chain` which contains only those residues from `hetatm_chain`
    which are less than `r_cutoff` away from some residue in `structure`.

    TODO: Should probably do something smarter about DNA and RNA.
    """
    # Refernce structure
    structure_atoms = list(structure.atoms)
    ns = NeighborSearch(structure_atoms)
    # Old HETATM chain
    hetatm_residues = []
    for residue in hetatm_chain.residues:
        if residue.resname in RESIDUE_MAPPING_TO_CANONICAL and residue.id[0] == " ":
            pass
        elif residue.resname in RNA_MAPPING_TO_CANONICAL and residue.id[0] == " ":
            # Skipping DNA
            pass
        elif residue.resname in DNA_MAPPING_TO_CANONICAL and residue.id[0] == " ":
            # Skipping RNA
            pass
        else:
            polymer_units = (
                set(RESIDUE_MAPPING_TO_CANONICAL)
                | set(RNA_MAPPING_TO_CANONICAL)
                | set(DNA_MAPPING_TO_CANONICAL)
            )
            if residue.resname in polymer_units or residue.id[0] == " ":
                logger.warning(
                    f"Encountered a strange residue with resname '{residue.resname}' and "
                    f"id '{residue.id}'. Treating as HETATM."
                )
            hetatm_residues.append(residue)
    # New HETATM chain (keeping only proximal HETATMs)
    new_hetatm_chain = Chain(hetatm_chain.id)
    new_hetatm_residues = [
        residue for residue in hetatm_residues if any(ns.search(a.coord, r_cutoff) for a in residue)
    ]
    for new_hetatm_residue in new_hetatm_residues:
        new_hetatm_chain.add(new_hetatm_residue)
    return new_hetatm_chain


def get_chain_sequence(chain: Chain, unknown_residue_marker: Optional[str] = None) -> str:
    chain_aa_list = []
    for residue in chain.residues:
        aaa = RESIDUE_MAPPING_TO_CANONICAL.get(residue.resname)
        if aaa is not None:
            aa = AAA_DICT[aaa]
        elif unknown_residue_marker is not None:
            aa = unknown_residue_marker
        else:
            continue
        chain_aa_list.append(aa)
    return "".join(chain_aa_list)


def chain_is_hetatm(chain: Chain) -> bool:
    """Return `True` if `chain` contains predominantly heteroatoms."""
    fraction_hetatm = sum(bool(residue.id[0].strip()) for residue in chain) / len(chain)
    return fraction_hetatm > 0.80
