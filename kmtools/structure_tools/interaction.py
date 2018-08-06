import logging
from collections import namedtuple
from typing import Iterable, List, Tuple, NamedTuple

import pandas as pd
from kmbio.PDB import NeighborSearch, Structure

from kmtools import df_tools, sequence_tools, structure_tools
from kmtools.structure_tools import AAA_DICT, COMMON_HETATMS

logger = logging.getLogger(__name__)

CORE_ID_COLUMNS = ["structure_id", "model_id", "chain_id"]
INTERFACE_ID_COLUMNS = ["structure_id", "model_id_1", "model_id_2", "chain_id_1", "chain_id_2"]


class Interaction(NamedTuple):
    structure_id: str
    model_id_1: str
    model_id_2: str
    chain_id_1: str
    chain_id_2: str
    residue_idx_1: str
    residue_idx_2: str
    residue_id_1: str
    residue_id_2: str
    residue_name_1: str
    residue_name_2: str
    residue_aa_1: str
    residue_aa_2: str


# #############################################################################
# Get interactions
# #############################################################################


def get_interactions(
    structure: Structure, r_cutoff: float, interchain: bool = True
) -> pd.DataFrame:
    """Return residue-residue interactions within and between chains in `structure`.

    Note:
        For each chain, ``residue.id[1]`` is unique.

    Todo:
        This could probably be sped up by using the
        :py:meth:`kmbio.PDB.NeighborSearch.search_all` method.

    Args:
        structure: Structure to analyse.
        interchain: Whether to include interactions between chains.

    Returns:
        A DataFrame with all interactions.
    """
    interactions = _get_interactions(structure, r_cutoff, interchain)
    # Specify `columns` in case the DataFrame is empty
    interaction_df = pd.DataFrame(interactions, columns=Interaction._fields)
    interaction_df = _add_reverse_interactions(interaction_df)
    return interaction_df


def _get_interactions(structure: Structure, r_cutoff: float, interchain: bool) -> List[Interaction]:
    """Compile a list of all interactions present in `structure`."""
    results = []
    for (
        model_1_idx,
        model_1,
        chain_1_idx,
        chain_1,
        chain_1_ns,
        model_2_idx,
        model_2,
        chain_2_idx,
        chain_2,
    ) in _iter_interchain_ns(structure, interchain=interchain):
        chain_1_residue_ids = [r.id for r in chain_1]
        for residue_2_idx, residue_2 in enumerate(chain_2):
            if residue_2.resname in COMMON_HETATMS:
                continue
            seen = set()
            chain_1_interacting_residues = [
                r
                for a in residue_2
                for r in chain_1_ns.search(a.coord, r_cutoff, "R")
                if r.id not in seen and not seen.add(r.id)
            ]
            for residue_1 in chain_1_interacting_residues:
                if residue_1.resname in COMMON_HETATMS:
                    continue
                residue_1_idx = chain_1_residue_ids.index(residue_1.id)
                if (
                    model_1.id == model_2.id
                    and chain_1.id == chain_2.id
                    and residue_1_idx >= residue_2_idx
                ):
                    continue
                row = Interaction(
                    structure.id,
                    model_1.id,
                    model_2.id,
                    chain_1.id,
                    chain_2.id,
                    residue_1_idx,
                    residue_2_idx,
                    residue_1.id[1],
                    residue_2.id[1],
                    residue_1.resname,
                    residue_2.resname,
                    AAA_DICT.get(residue_1.resname),
                    AAA_DICT.get(residue_2.resname),
                )
                results.append(row)
    return results


def _add_reverse_interactions(interaction_df: pd.DataFrame) -> pd.DataFrame:
    """Add reverse interactions to `interaction_df`.

    Args:
        interaction_df: DataFrame of residue-residue interactions.

    Returns:
        DataFrame of residue-residue interactions, where for each interactions
        between residue 1 and residue 2, there is a corresponding interaction between
        residue 2 and residue 1.
    """
    interaction_reverse_df = interaction_df.rename(columns=df_tools.reverse_column).reindex(
        columns=pd.Index(interaction_df.columns)
    )

    _length_before = len(interaction_df)
    df = pd.concat([interaction_df, interaction_reverse_df]).drop_duplicates()
    # Make sure that there were no duplicates
    assert len(df) == _length_before * 2
    # NB: Sorting will not work well for a mix of uppercase and lowercase columns
    df = df.sort_values(interaction_df.columns.tolist())
    df.index = range(len(df))
    return df


def _iter_interchain_ns(structure: Structure, interchain: bool = True) -> Iterable:
    """Iterate over interactions present in the `structure`."""
    for model_1_idx, model_1 in enumerate(structure):
        for chain_1_idx, chain_1 in enumerate(model_1):
            atom_list = list(chain_1.atoms)
            if not atom_list:
                logger.debug("Skipping %s %s because it is empty.", model_1, chain_1)
                continue
            chain_1_ns = NeighborSearch(atom_list)
            if not interchain:
                yield (
                    model_1_idx,
                    model_1,
                    chain_1_idx,
                    chain_1,
                    chain_1_ns,
                    model_1_idx,
                    model_1,
                    chain_1_idx,
                    chain_1,
                )
            else:
                for model_2_idx, model_2 in enumerate(structure):
                    if model_1_idx > model_2_idx:
                        continue
                    for chain_2_idx, chain_2 in enumerate(model_2):
                        if model_1_idx == model_2_idx and chain_1_idx > chain_2_idx:
                            continue
                        yield (
                            model_1_idx,
                            model_1,
                            chain_1_idx,
                            chain_1,
                            chain_1_ns,
                            model_2_idx,
                            model_2,
                            chain_2_idx,
                            chain_2,
                        )


# #############################################################################
# Process interactions
# #############################################################################


def process_interactions(interactions: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Split the `interactions` DataFrame into intra- and iner-chain interactions.

    Args:
        interactions: Interactions DataFrame, as returned by `get_interactions`.

    Returns:
        A DataFrame of core interactions and a DataFrame of interface interactions.
    """
    # Intra-chain interactions
    interactions_core = (
        interactions[
            (interactions["model_id_1"] == interactions["model_id_2"])
            & (interactions["chain_id_1"] == interactions["chain_id_2"])
        ]
        .drop(pd.Index(["model_id_2", "chain_id_2"]), axis=1)
        .rename(columns={"model_id_1": "model_id", "chain_id_1": "chain_id"})
        .copy()
    )

    # Inter-chain interactions
    interactions_interface = interactions[
        (interactions["model_id_1"] != interactions["model_id_2"])
        | (interactions["chain_id_1"] != interactions["chain_id_2"])
    ].copy()

    # Make sure everything adds up
    assert len(interactions) == len(interactions_core) + len(interactions_interface)
    return interactions_core, interactions_interface


def process_interactions_core(
    structure: Structure, interactions_core: pd.DataFrame
) -> pd.DataFrame:
    """Group core interactions by chain and add aggregate features.

    Args:
        structure: Structure for which the interactions were calculated.
        interactions_core: Core interactions between residues on the  same chain.

    Returns:
        A DataFrame of interface interactions grouped by chain.
    """
    if interactions_core.empty:
        return interactions_core[CORE_ID_COLUMNS]

    interactions_core = interactions_core.copy()
    interactions_core["residue_pair"] = interactions_core[["residue_idx_1", "residue_idx_2"]].apply(
        tuple, axis=1
    )

    interactions_core_aggbychain = (
        interactions_core.groupby(CORE_ID_COLUMNS)["residue_pair"]
        .agg(lambda x: tuple(x))
        .reset_index()
    )

    interactions_core_aggbychain["protein_sequence"] = [
        structure_tools.extract_aa_sequence(structure, model_id, chain_id)
        for model_id, chain_id in interactions_core_aggbychain[["model_id", "chain_id"]].values
    ]

    interactions_core_aggbychain["protein_sequence_hash"] = interactions_core_aggbychain[
        "protein_sequence"
    ].apply(sequence_tools.crc64)

    interactions_core_aggbychain["residue_sequence"] = [
        structure_tools.extract_residue_sequence(structure, model_id, chain_id)
        for model_id, chain_id in interactions_core_aggbychain[["model_id", "chain_id"]].values
    ]

    interactions_core_aggbychain["residue_pair_hash"] = interactions_core_aggbychain[
        "residue_pair"
    ].apply(structure_tools.hash_residue_pair)

    return interactions_core_aggbychain


def process_interactions_interface(
    structure: Structure, interactions_interface: pd.DataFrame
) -> pd.DataFrame:
    """Group interface interactions by chain pair and add aggregate features.

    Args:
        structure: Structure for which the interactions were calculated.
        interactions_interface: Interface interactions between residues on different chains.

    Returns:
        A DataFrame of interface interactions grouped by chain pair.
    """
    if interactions_interface.empty:
        return interactions_interface[INTERFACE_ID_COLUMNS]

    interactions_interface = interactions_interface.copy()
    interactions_interface["residue_pair"] = interactions_interface[
        ["residue_idx_1", "residue_idx_2"]
    ].apply(tuple, axis=1)

    interactions_interface_aggbychain = (
        interactions_interface.groupby(INTERFACE_ID_COLUMNS)["residue_pair"]
        .agg(lambda x: tuple(x))
        .reset_index()
    )

    def _extract_sequence(structure, model_id, chain_id):
        """Extract AA sequence if possible but fall back to residue sequence."""
        aa_sequence = structure_tools.extract_aa_sequence(structure, model_id, chain_id)
        if aa_sequence and aa_sequence.count("X") < (len(aa_sequence) // 3):
            return aa_sequence
        else:
            residue_sequence = structure_tools.extract_residue_sequence(
                structure, model_id, chain_id
            )
            return residue_sequence + ","

    # Interacting partner 1 properties
    interactions_interface_aggbychain["protein_sequence_1"] = [
        _extract_sequence(structure, model_id, chain_id)
        for model_id, chain_id in interactions_interface_aggbychain[
            ["model_id_1", "chain_id_1"]
        ].values
    ]

    interactions_interface_aggbychain[
        "protein_sequence_hash_1"
    ] = interactions_interface_aggbychain["protein_sequence_1"].apply(sequence_tools.crc64)

    # Interacting partner 2 properties
    interactions_interface_aggbychain["protein_sequence_2"] = [
        _extract_sequence(structure, model_id, chain_id)
        for model_id, chain_id in interactions_interface_aggbychain[
            ["model_id_2", "chain_id_2"]
        ].values
    ]

    interactions_interface_aggbychain[
        "protein_sequence_hash_2"
    ] = interactions_interface_aggbychain["protein_sequence_2"].apply(sequence_tools.crc64)

    # Interaction properties
    interactions_interface_aggbychain["residue_pair_hash"] = interactions_interface_aggbychain[
        "residue_pair"
    ].apply(structure_tools.hash_residue_pair)

    return interactions_interface_aggbychain


# #############################################################################
# Remove redundant rows from interaction dataframes
# #############################################################################


def drop_duplicates_core(
    interactions_core: pd.DataFrame, interactions_core_aggbychain: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Remove HETATM and duplicate chains and corresponding residues.

    Remove HETATM and duplicate chains from `interactions_core_aggbychain`,
    and remove corresponding residue pairs from `interactions_core`.

    Preconditions:
        - `interactions_core_aggbychain` should already be sorted from highest priority
          (i.e. removed last) to lowest priority (i.e. removed first).

    Args:
        interactions_core: DataFrame of interactions between residues on the same chains
            (as obtained using `process_interactions`).
        interactions_core_aggbychain: DataFrame of interactions within different chains
            (as obtained using `process_interactions_core`).

    Returns:
        Copy of `interactions_core` and `interactions_core_aggbychain`
        with unnecessary rows removed.
    """
    if interactions_core.empty:
        return interactions_core, interactions_core_aggbychain

    if "protein_sequence" not in interactions_core_aggbychain and (
        "protein_sequence_1" in interactions_core_aggbychain
        and "protein_sequence_2" in interactions_core_aggbychain
    ):
        raise Exception(
            "Looks like you are running `drop_duplicates_core` for **interface** interactions."
        )

    # Remove hetatm chains
    _before = len(interactions_core_aggbychain)
    interactions_core_aggbychain = interactions_core_aggbychain[
        interactions_core_aggbychain["protein_sequence"].notnull()
    ].copy()
    logger.debug("Removed %s hetatm chains!", _before - len(interactions_core_aggbychain))

    # Remove duplicates
    _before = interactions_core_aggbychain.shape[0]
    interactions_core_aggbychain = interactions_core_aggbychain.drop_duplicates(
        subset=["protein_sequence"]
    )
    logger.debug("Removed %s duplicated chains!", _before - len(interactions_core_aggbychain))

    # Remove corresponding residues from `interactions_core`
    _before = len(interactions_core)
    interactions_core = interactions_core.merge(
        interactions_core_aggbychain[structure_tools.CORE_ID_COLUMNS],
        on=structure_tools.CORE_ID_COLUMNS,
    )
    logger.debug(
        "Lost %s residues from hetatm or duplicated chains!", _before - len(interactions_core)
    )

    return interactions_core, interactions_core_aggbychain


def drop_duplicates_interface(
    interactions_interface: pd.DataFrame, interactions_interface_aggbychain: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Remove HETATM and duplicate chain pairs and corresponding residue pairs.

    Remove HETATM and duplicate chain pairs from `interactions_interface_aggbychain`,
    and remove corresponding residue pairs from `interactions_interface`.

    Preconditions:
        - `interactions_core_aggbychain` should already be sorted from highest priority
          (i.e. removed last) to lowest priority (i.e. removed first).

    Args:
        interactions_interface: DataFrame of interactions between residues on different chains
            (as obtained using `process_interactions`).
        interactions_interface_aggbychain: DataFrame of interactions between different chains
            (as obtained using `process_interactions_interface`).

    Returns:
        Copy of `interactions_interface` and `interactions_interface_aggbychain`
        with unnecessary rows removed.
    """
    if interactions_interface.empty:
        return interactions_interface, interactions_interface_aggbychain

    if "protein_sequence" in interactions_interface_aggbychain and not (
        "protein_sequence_1" in interactions_interface_aggbychain
        and "protein_sequence_2" in interactions_interface_aggbychain
    ):
        raise Exception(
            "Looks like you are running `drop_duplicates_interface` for **core** interactions."
        )

    # Make a copy that we will modify (prevents annothing 'indexing-view-versus-copy' errors)
    interactions_interface_aggbychain = interactions_interface_aggbychain.copy()

    # Remove rows that correpond to interactions between existing chain pairs
    _before = len(interactions_interface_aggbychain)
    interactions_interface_aggbychain = interactions_interface_aggbychain[
        ~(
            (interactions_interface_aggbychain["protein_sequence_1"].str.contains(","))
            & (interactions_interface_aggbychain["protein_sequence_2"].str.contains(","))
        )
    ].copy()
    logger.debug("Removed %s HETATM chain pairs!", _before - len(interactions_interface_aggbychain))

    if interactions_interface_aggbychain.empty:
        return (interactions_interface.iloc[0:0], interactions_interface_aggbychain.iloc[0:0])

    # Remove rows that correpond to interactions between existing chain pairs
    _before = len(interactions_interface_aggbychain)
    interactions_interface_aggbychain["unique_key"] = interactions_interface_aggbychain[
        ["protein_sequence_1", "protein_sequence_2"]
    ].apply(lambda x: x[0] + x[1] if x[0] <= x[1] else x[1] + x[0], axis=1)
    interactions_interface_aggbychain = interactions_interface_aggbychain.drop_duplicates(
        subset=["unique_key"]
    )
    interactions_interface_aggbychain.drop(pd.Index(["unique_key"]), axis=1, inplace=True)
    logger.debug("Removed %s chain pairs!", _before - len(interactions_interface_aggbychain))

    # Remove corresponding rows
    _before = len(interactions_interface)
    interactions_interface = interactions_interface.merge(
        interactions_interface_aggbychain[structure_tools.INTERFACE_ID_COLUMNS],
        on=structure_tools.INTERFACE_ID_COLUMNS,
    )
    logger.debug("Lost %s interacting residue pairs!", _before - len(interactions_interface))

    return interactions_interface, interactions_interface_aggbychain
