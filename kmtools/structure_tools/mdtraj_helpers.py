import tempfile
from pathlib import Path
from typing import Union

import mdtraj
import numpy as np
import pandas as pd


def mdtraj_to_dataframe(traj: mdtraj.Trajectory) -> pd.DataFrame:
    """Generate a structure dataframe that is backwards-compatible with kmbio.

    Args:
        traj: An mdtraj trajectory.

    Returns:
        A dataframe with information describing the structure.
    """
    structure_df = traj.top.to_dataframe()[0]

    # Residue
    residue_map = {
        (chainID, resSeq): i
        for i, (chainID, resSeq) in enumerate(
            structure_df[["chainID", "resSeq"]].drop_duplicates().values
        )
    }
    structure_df["residue_idx"] = structure_df.apply(
        lambda s: residue_map[s.chainID, s.resSeq], axis=1
    )

    # Chain
    chain_map = {
        chainID: i for i, chainID in enumerate(structure_df["chainID"].drop_duplicates().values)
    }
    structure_df["chain_idx"] = structure_df["chainID"].map(chain_map)

    # Model
    structure_df["model_idx"] = 0

    structure_df["atom_serial_number"] = structure_df["serial"]
    structure_df["atom_name"] = structure_df["name"]
    structure_df["residue_id_1"] = structure_df["resSeq"]
    structure_df["residue_resname"] = structure_df["resName"]
    structure_df["atom_idx"] = np.arange(len(structure_df))
    # Convert nm back to Angstrom
    structure_df["atom_x"], structure_df["atom_y"], structure_df["atom_z"] = (
        traj.xyz[0].astype(np.double).T * 10
    )
    return structure_df


def mdtraj_to_pdb(traj: mdtraj.Trajectory, output_file: Union[str, Path]) -> None:
    with tempfile.NamedTemporaryFile(suffix=".pdb") as tmp_file:
        traj.save_pdb(tmp_file.name)
        with open(tmp_file.name, "rt") as fin:
            # For some reason, having 'MODEL' lines confuse reduce.
            lines = [
                r
                for r in fin.read().split("\n")
                if not r.startswith("MODEL") and not r.startswith("ENDMDL")
            ]
    with open(output_file, "wt") as fout:
        fout.write("\n".join(lines))
