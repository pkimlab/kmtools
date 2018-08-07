import hashlib
import logging
from pathlib import Path
from typing import Optional

import kmbio.PDB
import numpy as np
from kmbio.PDB import Structure

from kmtools.structure_tools import A_DICT, AAA_DICT, PDB_RESIDUES, RESIDUE_MAPPING_TO_CANONICAL

logger = logging.getLogger(__name__)


def extract_aa_sequence(structure: Structure, model_id: int, chain_id: str) -> Optional[str]:
    """Return amino acid sequence of all residues in chain.

    Residues which cannot be assigned a single-character amino acid code represented as 'X'.

    Args:
        structure: Structure from which to extract the chain sequence.
        model_id: ID of the model which contains the chain of interest.
        chain_id: ID of the chain of interest.

    Returns:
        String of amino acids corresponding to the chain sequence.
        The length of the returned string must equal to the number of residues in the chain.
        Residues which cannot be mapped to a single-character amino acid are represented as 'X'.
    """
    aa_list = []
    for residue in structure[model_id][chain_id]:
        if residue.resname in AAA_DICT:
            residue_aa = AAA_DICT[residue.resname]
        elif residue.resname in RESIDUE_MAPPING_TO_CANONICAL:
            residue_aa = AAA_DICT[RESIDUE_MAPPING_TO_CANONICAL[residue.resname]]
        else:
            logger.debug(
                "Cannot map residue %s to a single-character amino acid string.", residue.resname
            )
            residue_aa = "X"
        aa_list.append(residue_aa)
    aa_string = "".join(aa_list)
    assert len(aa_string) == len(list(structure[model_id][chain_id].residues))
    return aa_string if aa_string else None


def extract_residue_sequence(structure: Structure, model_id: int, chain_id: str) -> Optional[str]:
    """Return comma-delimited residue sequence of all residues in chain."""
    aa_string = ",".join(r.resname for r in structure[model_id][chain_id])
    return aa_string if aa_string else None


def hash_residue_pair(residue_pair) -> str:
    """Create a hash of a pair of residues.

    Note:
        This function should give the same result if the order of the residues is switched.

    Examples:
        >>> hash_residue_pair(('G', 'A'))
        '14965c5e8f0e81e3cfc839f6cfd73989'
        >>> hash_residue_pair(('A', 'G'))
        '14965c5e8f0e81e3cfc839f6cfd73989'
    """
    myhash = hashlib.md5()
    a = np.array(sorted(residue_pair))
    myhash.update(a)
    return myhash.hexdigest()


def standardize_mutation(pdb_file: Path, pdb_mutation: str) -> str:
    """Convert a residue RESNUM to a residue ID."""
    structure = kmbio.PDB.load(pdb_file)

    if "-" in pdb_mutation:
        pdb_chain, mutation = pdb_mutation.split("-")
    elif "_" in pdb_mutation:
        pdb_chain, mutation = pdb_mutation.split("_")
    else:
        raise Exception(
            "`pdb_mutation` must contain the chain and the mutation separated by '-' or '_'!"
        )

    try:
        chain = structure[0][pdb_chain]
    except KeyError:
        raise Exception(f"Could not find chain '{pdb_chain}' in structure '{structure.id}'!")

    found = False
    residue_idx = 0
    for residue in chain:
        residue_id = f"{residue.id[1]}{residue.id[2]}"
        if residue_id == mutation[1:-1]:
            found = residue.resname == A_DICT[pdb_mutation[0]]
            break
        if residue.resname in PDB_RESIDUES:
            residue_idx += 1

    if found:
        return f"{pdb_chain}-{pdb_mutation[0]}{residue_idx + 1}{pdb_mutation[-1]}"
    else:
        return None
