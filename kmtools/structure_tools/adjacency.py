from typing import List, Optional

import numpy as np
import pandas as pd
from kmbio.PDB import Chain, Model, Structure
from scipy.spatial import cKDTree

from .types import DomainDef


def extract_domain(
    structure: Structure, dds: List[DomainDef], remove_hetatms=False, hetatm_residue_cutoff=None
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
        if remove_hetatms:
            domain_residues = [r for r in domain_residues if not r.id[0].strip()]
        new_chain.add(domain_residues)
    return new_structure


def get_distances(
    structure_df: pd.DataFrame, max_cutoff: Optional[float], groupby: str = "atom"
) -> pd.DataFrame:
    """Process structure dataframe to extract interacting chains, residues, or atoms.

    Args:
        structure_df: Structure dataframe, as returned by `kmbio.PDB.Structure.to_dataframe`.
        max_cutoff: Maximum distance for inclusion in results (Angstroms).
        groupby: Which pairs of objects to return?
            (Possible options are: `chain{,-backbone,-ca}`, `residue{,-backbone,-ca}`, `atom`).
    """
    assert groupby in [
        "chain",
        "chain-backbone",
        "chain-ca",
        "chain-cb",
        "residue",
        "residue-backbone",
        "residue-ca",
        "residue-cb",
        "atom",
    ]

    residue_idxs = set(structure_df["residue_idx"])
    if groupby.endswith("-backbone"):
        structure_df = structure_df[
            (structure_df["atom_name"] == "N")
            | (structure_df["atom_name"] == "CA")
            | (structure_df["atom_name"] == "C")
        ]
    elif groupby.endswith("-ca"):
        structure_df = structure_df[(structure_df["atom_name"] == "CA")]
    elif groupby.endswith("-cb"):
        structure_df = structure_df[(structure_df["atom_name"] == "CB")]
    assert not residue_idxs - set(structure_df["residue_idx"])

    pairs_df = get_atom_distances(structure_df, max_cutoff=max_cutoff)
    annotate_atom_distances(pairs_df, structure_df)

    if groupby.startswith("residue"):
        pairs_df = _groupby_residue(pairs_df)
    elif groupby.startswith("chain"):
        pairs_df = _groupby_chain(pairs_df)
    return pairs_df


def get_atom_distances(structure_df: pd.DataFrame, max_cutoff: Optional[float]) -> pd.DataFrame:
    if max_cutoff is None:
        return get_atom_distances_dense(structure_df)
    else:
        return get_atom_distances_sparse(structure_df, max_cutoff)


def get_atom_distances_dense(structure_df: pd.DataFrame) -> pd.DataFrame:
    num_atoms = len(structure_df)
    xyz = structure_df[["atom_x", "atom_y", "atom_z"]].values
    xyz_by_xyz = xyz[:, None, :] - xyz[None, :, :]
    xyz_by_xzy_dist = np.sqrt((xyz_by_xyz ** 2).sum(axis=2))
    assert xyz_by_xzy_dist.ndim == 2
    row = np.repeat(np.arange(num_atoms), num_atoms)
    col = np.tile(np.arange(num_atoms), num_atoms)
    atom_idx = structure_df["atom_idx"].values
    pairs_df = pd.DataFrame(
        {
            "atom_idx_1": atom_idx[row],
            "atom_idx_2": atom_idx[col],
            "distance": xyz_by_xzy_dist.reshape(-1),
        }
    )
    pairs_df = pairs_df[pairs_df["atom_idx_1"] < pairs_df["atom_idx_2"]]
    return pairs_df


def get_atom_distances_sparse(structure_df: pd.DataFrame, max_cutoff: float) -> pd.DataFrame:
    xyz = structure_df[["atom_x", "atom_y", "atom_z"]].values
    tree = cKDTree(xyz)
    coo_mat = tree.sparse_distance_matrix(tree, max_distance=max_cutoff, output_type="coo_matrix")
    assert coo_mat.row.max() < xyz.shape[0] and coo_mat.col.max() < xyz.shape[0]
    atom_idx = structure_df["atom_idx"].values
    pairs_df = pd.DataFrame(
        {
            "atom_idx_1": atom_idx[coo_mat.row],
            "atom_idx_2": atom_idx[coo_mat.col],
            "distance": coo_mat.data,
        }
    )
    pairs_df = pairs_df[pairs_df["atom_idx_1"] < pairs_df["atom_idx_2"]]
    return pairs_df


def annotate_atom_distances(pairs_df: pd.DataFrame, structure_df: pd.DataFrame) -> None:
    atom_to_residue_map = {
        row.atom_idx: (row.model_idx, row.chain_idx, row.residue_idx)
        for row in structure_df.itertuples()
    }
    for suffix in ["_1", "_2"]:
        (
            pairs_df[f"model_idx{suffix}"],
            pairs_df[f"chain_idx{suffix}"],
            pairs_df[f"residue_idx{suffix}"],
        ) = list(zip(*pairs_df[f"atom_idx{suffix}"].map(atom_to_residue_map)))


def _groupby_residue(pairs_df: pd.DataFrame) -> pd.DataFrame:
    residue_pairs_df = (
        pairs_df.groupby(["residue_idx_1", "residue_idx_2"]).agg({"distance": min}).reset_index()
    )
    residue_pairs_df = residue_pairs_df[
        (residue_pairs_df["residue_idx_1"] != residue_pairs_df["residue_idx_2"])
    ]
    return residue_pairs_df


def _groupby_chain(pairs_df: pd.DataFrame) -> pd.DataFrame:
    chain_pairs_df = (
        pairs_df.groupby(["chain_idx_1", "chain_idx_2"]).agg({"distance": min}).reset_index()
    )
    chain_pairs_df = chain_pairs_df[
        (chain_pairs_df["chain_idx_1"] != chain_pairs_df["chain_idx_2"])
    ]
    return chain_pairs_df


def complete_distances(distance_df: pd.DataFrame) -> pd.DataFrame:
    """Complete `distances_df` so that it corresponds to a symetric dataframe with a zero diagonal.

    Args:
        distance_df: A dataframe produced by `get_distances`, where
            'residue_idx_1' > 'residue_idx_2'.

    Returns:
        A dataframe which corresponds to a dense, symmetric adjacency matrix.
    """
    residues = sorted(
        {r for r in distance_df["residue_idx_1"]} | {r for r in distance_df["residue_idx_2"]}
    )
    complete_distance_df = pd.concat(
        [
            distance_df,
            distance_df.rename(
                columns={"residue_idx_1": "residue_idx_2", "residue_idx_2": "residue_idx_1"}
            ),
            pd.DataFrame({"residue_idx_1": residues, "residue_idx_2": residues, "distance": 0}),
        ],
        sort=False,
    ).sort_values(["residue_idx_1", "residue_idx_2"])
    assert len(complete_distance_df) == len(distance_df) * 2 + len(residues)
    return complete_distance_df
