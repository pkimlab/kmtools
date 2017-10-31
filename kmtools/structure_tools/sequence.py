import hashlib
import warnings
from typing import Optional

import numpy as np

from kmbio.PDB import Structure
from kmtools.structure_tools import AAA_DICT


def extract_aa_sequence(structure: Structure, model_id: int, chain_id: str) -> Optional[str]:
    """Return amino acid sequence of all residues in chain.

    Residues which cannot be assigned a single-character amino acid code are skipped.
    """
    aa_list = []
    skipped_resnames = set()
    for residue in structure[model_id][chain_id]:
        try:
            aa_list.append(AAA_DICT[residue.resname])
        except KeyError:
            skipped_resnames.add(residue.resname)
    if skipped_resnames:
        warning = ("Skipped the following residues when generating chain sequence: {}"
                   .format(skipped_resnames))
        warnings.warn(warning)
    aa_string = ''.join(aa_list)
    return aa_string if aa_string else None


def extract_residue_sequence(structure: Structure, model_id: int, chain_id: str) -> Optional[str]:
    """Return comma-delimited residue sequence of all residues in chain."""
    aa_string = ','.join(r.resname for r in structure[model_id][chain_id])
    return aa_string if aa_string else None


def hash_residue_pair(residue_pair):
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
