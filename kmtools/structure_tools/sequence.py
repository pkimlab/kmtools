import hashlib
import logging
from pathlib import Path
from typing import Optional, Union

import kmbio.PDB
import numpy as np
from kmbio.PDB import Structure

from kmtools.structure_tools import AAA_DICT, RESIDUE_MAPPING_TO_CANONICAL

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


def standardize_mutation(structure: Union[str, Path, Structure], pdb_mutation: str) -> str:
    """Convert mutation from using residue resnum to using residue id.

    The residue id is a sequential number that starts at 1.

    Args:
        structure: Location of the PDB or mmCIF file describing the structure
            or the structure iteself.
        pdb_mutation: PDB chain and mutation for a single residue.

    Returns:
        Mutation in the same format as `pdb_mutation` but using the standardized residue id.

    Raises:
        ValueError: If `pdb_mutation` does not match the expected format.
        KeyError: If the chain specified in `pdb_mutation` is not found in the provided structure.
        IndexError: If we fail to map `pdb_mutation` to a residue in the structure,
            for whatever reason.

    Examples:
        >>> standardize_mutation("rscb://1ak4.pdb", "D-A488G")
        D_A88G
    """
    if isinstance(structure, (Path, str)):
        structure = kmbio.PDB.load(structure)

    seps = ["_", "-"]
    for sep in seps + [None]:
        if sep is None:
            raise ValueError(
                f"`pdb_mutation` must contain the chain and the mutation separated one of: {seps}"
            )
        if sep in pdb_mutation:
            pdb_chain, mutation = pdb_mutation.split(sep)
            break

    try:
        chain = structure[0][pdb_chain]
    except KeyError:
        raise KeyError(f"Could not find chain '{pdb_chain}' in structure '{structure.id}'!")

    residue_idx = 0
    for residue in chain:
        if residue.resname not in RESIDUE_MAPPING_TO_CANONICAL:
            continue
        if (str(residue.id[1]) + residue.id[2].strip()) == mutation[1:-1]:
            ref_wt = mutation[0]
            mapped_wt = AAA_DICT[RESIDUE_MAPPING_TO_CANONICAL[residue.resname]]
            if mapped_wt != ref_wt:
                raise IndexError(
                    f"Reference AA '{ref_wt}' does not match the mapped AA '{mapped_wt}' "
                    f"for residue '{residue.id}' and mutation '{pdb_mutation}'!"
                )
            else:
                return f"{pdb_chain}_{mutation[0]}{residue_idx + 1}{mutation[-1]}"
        residue_idx += 1
    raise IndexError(
        f"Could not map mutation {pdb_mutation} to a residue in the provided structure!"
    )
