"""
Functions for creating a dataset of intra-chain and inter-chain residue-residue interactions
for all proteins in the Protein Data Bank.
"""
import logging
from typing import Tuple

import pandas as pd

import kmbio.PDB
from kmtools import structure_tools

logger = logging.getLogger(__name__)


def generate_interaction_dataset(structure_file: str, r_cutoff: float
                                ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Args:
        structure_file: Filename of the structure file (preferrably in mmCIF format).
        r_cutoff: Two residues are considered to be interacting if they are within
            this number of Ångströms from each other.

    Returns:
        A tuple of four dataframes containing the intra-chain and inter-chain interactions,
        grouped by residue and by chain.
    """
    try:
        structure = kmbio.PDB.load(structure_file, bioassembly_id=1, use_auth_id=False)
    except kmbio.PDB.exceptions.BioassemblyError as e:
        logger.error("Could not construct bioassembly for file '%s'", structure_file)
        structure = kmbio.PDB.load(structure_file, bioassembly_id=0, use_auth_id=False)

    interactions = structure_tools.get_interactions(structure, r_cutoff=r_cutoff, interchain=True)
    interactions_core, interactions_interface = structure_tools.process_interactions(interactions)

    # Group interactions by chain / chan pair
    interactions_core_aggbychain = structure_tools.process_interactions_core(
        structure, interactions_core)
    interactions_interface_aggbychain = structure_tools.process_interactions_interface(
        structure, interactions_interface)

    # Drop duplicate rows?
    interactions_core, interactions_core_aggbychain = \
        structure_tools.drop_duplicates_core(interactions_core, interactions_core_aggbychain)
    interactions_interface, interactions_interface_aggbychain = \
        structure_tools.drop_duplicates_interface(interactions_interface,
                                                  interactions_interface_aggbychain)

    return (interactions_core, interactions_core_aggbychain, interactions_interface,
            interactions_interface_aggbychain)
