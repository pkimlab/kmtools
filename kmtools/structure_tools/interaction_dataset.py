import hashlib
import logging
import os.path as op
import warnings
from typing import Tuple

import numpy as np
import pandas as pd

import kmbio.PDB
from kmbio.PDB import Structure
from kmtools import sequence_tools, structure_tools
from kmtools.structure_tools import AAA_DICT

logger = logging.getLogger(__name__)

CORE_ID_COLUMNS = ['structure_id', 'model_id', 'chain_id']
INTERFACE_ID_COLUMNS = ['structure_id', 'model_id_1', 'model_id_2', 'chain_id_1', 'chain_id_2']


def generate_interaction_dataset(filename: str, r_cutoff: float):
    pdb_id = op.basename(filename)[:4]
    pdb_path = op.dirname(filename)

    structure = kmbio.PDB.load(filename, type='cif', biounit_id=1)

    interactions = structure_tools.get_interactions(structure, r_cutoff=r_cutoff, interchain=True)
    interactions_core, interactions_interface = process_interactions(interactions)

    interactions_core_aggbychain = process_interactions_core(structure, interactions_core)
    interactions_interface_aggbychain = process_interactions_interface(
        structure, interactions_interface)

    interactions_core, interactions_core_aggbychain = drop_duplicates_core(
        *process_interactions_core(structure, interactions_core))

    interactions_interface, interactions_interface_aggbychain = drop_duplicates_interface(
        *process_interactions_interface(structure, interactions_interface))

    save_csv(interactions_core, op.join(pdb_path, '{}_core_chain.tsv.gz'.format(pdb_id)))
    save_csv(interactions_core_aggbychain,
             op.join(pdb_path, '{}_core_residue.tsv.gz'.format(pdb_id)))
    save_csv(interactions_interface, op.join(pdb_path, '{}_interface_chain.tsv.gz'.format(pdb_id)))
    save_csv(interactions_interface_aggbychain,
             op.join(pdb_path, '{}_interface_residue.tsv.gz'.format(pdb_id)))


def process_interactions(interactions: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Split the `interactions` DataFrame into intra- and iner-chain interactions."""
    # Intra-chain interactions
    interactions_core = \
        interactions[
            (interactions['model_id_1'] == interactions['model_id_2']) &
            (interactions['chain_id_1'] == interactions['chain_id_2'])
        ] \
        .drop(pd.Index(['model_id_2', 'chain_id_2']), axis=1) \
        .rename(columns={'model_id_1': 'model_id', 'chain_id_1': 'chain_id'}) \
        .assign(residue_pair=lambda df: (df.residue_id_1, df.residue_id_2))

    # Inter-chain interactions
    interactions_interface = \
        interactions[
            (interactions['model_id_1'] != interactions['model_id_2']) |
            (interactions['chain_id_1'] != interactions['chain_id_2'])
        ] \
        .assign(residue_pair=lambda df: (df.residue_id_1, df.residue_id_2))

    # Make sure everything adds up
    assert len(interactions) == len(interactions_core) + len(interactions_interface)
    return interactions_core, interactions_interface


def process_interactions_core(structure: Structure, interactions_core: pd.DataFrame
                             ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Calculate aggregate features for interactions within each chain."""
    interactions_core_aggbychain = \
        interactions_core \
        .groupby(CORE_ID_COLUMNS, as_index=False)['residue_pair'] \
        .agg(lambda x: tuple(x))

    interactions_core_aggbychain['protein_sequence'] = [
        extract_aa_sequence(structure, model_id, chain_id)
        for (model_id, chain_id) in interactions_core_aggbychain[['model_id', 'chain_id']].values
    ]

    interactions_core_aggbychain['protein_sequence_hash'] = \
        interactions_core_aggbychain['protein_sequence'] \
        .apply(sequence_tools.crc64)

    interactions_core_aggbychain['residue_sequence'] = [
        extract_residue_sequence(structure, model_id, chain_id)
        for model_id, chain_id in interactions_core_aggbychain[['model_id', 'chain_id']].values
    ]

    interactions_core_aggbychain['residue_pair_hash'] = \
        interactions_core_aggbychain['residue_pair'] \
        .apply(hash_residue_pair)

    return interactions_core, interactions_core_aggbychain


def process_interactions_interface(structure: Structure,
                                   interactions_interface: pd.DataFrame) -> pd.DataFrame:
    """Calculate aggregate features for interactions between each pair of chains."""
    if interactions_interface.empty:
        return interactions_interface[INTERFACE_ID_COLUMNS]

    assert 'residue_pair' in interactions_interface

    interactions_interface_aggbychain = \
        interactions_interface \
        .groupby(INTERFACE_ID_COLUMNS, as_index=False)['residue_pair'] \
        .agg(lambda x: tuple(x))

    # Interacting partner 1 properties
    interactions_interface_aggbychain['protein_sequence_1'] = [
        extract_aa_sequence(structure, model_id, chain_id)
        for model_id, chain_id in interactions_interface_aggbychain[['model_id_1', 'chain_id_1']]
        .values
    ]

    interactions_interface_aggbychain['protein_sequence_hash_1'] = \
        interactions_interface_aggbychain['protein_sequence_1'] \
        .apply(sequence_tools.crc64)

    # Interacting partner 2 properties
    interactions_interface_aggbychain['protein_sequence_2'] = [
        extract_aa_sequence(structure, model_id, chain_id)
        for model_id, chain_id in interactions_interface_aggbychain[['model_id_2', 'chain_id_2']]
        .values
    ]

    interactions_interface_aggbychain['protein_sequence_hash_2'] = \
        interactions_interface_aggbychain['protein_sequence_2'] \
        .apply(sequence_tools.crc64)

    # Interaction properties
    interactions_interface_aggbychain['residue_pair_hash'] = \
        interactions_interface_aggbychain['residue_pair'] \
        .apply(hash_residue_pair)

    return interactions_interface_aggbychain


def drop_duplicates_core(
        interactions_core: pd.DataFrame,
        interactions_core_aggbychain: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """

    Todo:
        * Figure out if we actually need this...
    """
    # Remove hetatm chain
    if interactions_core_aggbychain.iloc[-1]['protein_sequence'] is None:
        logger.info("Removing last HETATM chain.")
        interactions_core_aggbychain = interactions_core_aggbychain.iloc[:-1]

    # Remove duplicates
    _before = interactions_core_aggbychain.shape[0]
    interactions_core_aggbychain = interactions_core_aggbychain.drop_duplicates()
    logger.info("Lost %s core chains!", _before - len(interactions_core_aggbychain))

    _before = len(interactions_core)
    interactions_core = \
        interactions_core \
        .merge(interactions_core_aggbychain[CORE_ID_COLUMNS], on=CORE_ID_COLUMNS)
    logger.info("Lost %s core residues!", _before - len(interactions_core))

    return interactions_core, interactions_core_aggbychain


def drop_duplicates_interface(
        interactions_interface: pd.DataFrame,
        interactions_interface_aggbychain: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """

    Todo:
        * Figure out if we actually need this...
    """
    if interactions_interface.empty:
        return interactions_interface, interactions_interface_aggbychain

    _before = len(interactions_interface_aggbychain)
    interactions_interface_aggbychain = interactions_interface_aggbychain.drop_duplicates()
    logger.info("Lost %s chain pairs!", _before - len(interactions_interface_aggbychain))

    _before = len(interactions_interface)
    interactions_interface = \
        interactions_interface \
        .merge(interactions_interface_aggbychain[INTERFACE_ID_COLUMNS], on=INTERFACE_ID_COLUMNS)
    logger.info("Lost %s interacting residue pairs!", _before - len(interactions_interface))

    return interactions_interface, interactions_interface_aggbychain


def extract_aa_sequence(structure: Structure, model_id: int, chain_id: str) -> str:
    """Return amino acid sequence of all residues in chain.

    Residues which cannot be assigned a single-character amino acid code are skipped.
    """
    aa_list = []
    skipped_resname_list = set()
    for residue in structure[model_id][chain_id]:
        try:
            aa_list.append(AAA_DICT[residue.resname])
        except KeyError:
            skipped_resname_list.add(residue.resname)
    if skipped_resname_list:
        warnings.warn("Skipped the following residues when generating chain sequence:\n"
                      "{}".format(skipped_resname_list))
    aa_string = ''.join(aa_list)
    return aa_string


def extract_residue_sequence(structure: Structure, model_id: int, chain_id: str) -> str:
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


def save_csv(df: pd.DataFrame, filename: str) -> None:
    """Save residue-residue interaction data as a compressed CSV."""
    if df.empty:
        open(filename, 'w').close()
    else:
        columns_to_ignore = ['residue_pair']
        df = df.drop(pd.Index(columns_to_ignore), axis=1)
        df.to_csv(filename, sep='\t', na_rep='\\N', header=False, index=False, compression='gzip')
