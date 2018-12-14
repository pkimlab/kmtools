from typing import List

import pandas as pd
from kmbio.PDB import Chain, Model, Structure
from MDAnalysis.lib.distances import self_capped_distance

from .types import DomainDef


def get_distances(structure: Structure, max_cutoff: int, min_cutoff=None, groupby="atom"):
    assert groupby in ["chain", "residue", "atom"]
    df = structure.to_dataframe()
    atom_to_residue_map = {
        row.atom_idx: (row.model_idx, row.chain_idx, row.residue_idx) for row in df.itertuples()
    }
    pairs_df = get_atom_distances(
        df[["atom_x", "atom_y", "atom_z"]].values, max_cutoff=max_cutoff, min_cutoff=min_cutoff
    )
    for suffix in ["_1", "_2"]:
        (
            pairs_df[f"model_idx{suffix}"],
            pairs_df[f"chain_idx{suffix}"],
            pairs_df[f"residue_idx{suffix}"],
        ) = list(zip(*pairs_df[f"atom_idx{suffix}"].map(atom_to_residue_map)))
    if groupby == "residue":
        pairs_df = _groupby_residue(pairs_df)
    elif groupby == "chain":
        pairs_df = _groupby_chain(pairs_df)
    return pairs_df


def _groupby_residue(pairs_df):
    residue_pairs_df = (
        pairs_df.groupby(["residue_idx_1", "residue_idx_2"]).agg({"distance": min}).reset_index()
    )
    residue_pairs_df = residue_pairs_df[
        (residue_pairs_df["residue_idx_1"] != residue_pairs_df["residue_idx_2"])
    ]
    return residue_pairs_df


def _groupby_chain(pairs_df):
    chain_pairs_df = (
        pairs_df.groupby(["chain_idx_1", "chain_idx_2"]).agg({"distance": min}).reset_index()
    )
    chain_pairs_df = chain_pairs_df[
        (chain_pairs_df["chain_idx_1"] != chain_pairs_df["chain_idx_2"])
    ]
    return chain_pairs_df


def get_atom_distances(xyz, max_cutoff, min_cutoff=None):
    pairs, distances = self_capped_distance(xyz, max_cutoff=max_cutoff, min_cutoff=min_cutoff)
    pairs.sort(axis=1)
    assert pairs.max() < xyz.shape[0]
    pairs_df = pd.DataFrame(
        {"atom_idx_1": pairs[:, 0], "atom_idx_2": pairs[:, 1], "distance": distances}
    )
    return pairs_df


def extract_domain(
    structure: Structure, dds: List[DomainDef], hetatm_residue_cutoff=None
) -> Structure:
    if hetatm_residue_cutoff is not None:
        # Keep residues from HETATM chain that are within `hetatm_residue_cutoff`
        # from any atom in the extracted chain.
        raise NotImplementedError
    assert len({(dd.model_id, dd.chain_id) for dd in dds}) == len(dds)
    new_structure = Structure(structure.id)
    new_model = Model(0)
    new_structure.add(new_model)
    for dd in dds:
        new_chain = Chain(dd.chain_id)
        new_model.add(new_chain)
        residues = list(structure[dd.model_id][dd.chain_id].residues)
        domain_residues = residues[dd.domain_start - 1 : dd.domain_end]
        new_chain.add(domain_residues)
    return new_structure
