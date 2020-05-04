import contextlib
import logging
import shlex
import subprocess
from pathlib import Path
from typing import IO, Optional, Union

import numpy as np
import pandas as pd
from scipy import stats

logger = logging.getLogger(__name__)


def protonate(
    input_file: Union[str, Path, IO], output_file: Union[str, Path, IO], method: str = "openmm"
) -> None:
    """Add hydrogen atoms to PDB structure.

    Args:
        input_file: Input PDB file or IO object.
        output_file: Output PDB file or IO object.
        method: Method to use for adding hydrogens. Supported methods are "openmm" and "reduce".
    """
    supported_methods = ["openmm", "reduce", "reduce-FLIP", "reduce-NOFLIP", "reduce-BUILD"]
    if method not in supported_methods:
        raise ValueError(f"The only supported methods are '{supported_methods}'.")
    with contextlib.ExitStack() as stack:
        if isinstance(input_file, (str, Path)):
            input_file = stack.enter_context(open(input_file, "rt"))
        if isinstance(output_file, (str, Path)):
            output_file = stack.enter_context(open(output_file, "wt"))
        if method.startswith("openmm"):
            _protonate_with_openmm(input_file, output_file)
        elif method.startswith("reduce"):
            _protonate_with_reduce(input_file, output_file, method="")


def _protonate_with_reduce(input_file: IO, output_file: IO, method: str = "") -> None:
    assert method in ["", "FLIP", "NOFLIP", "BUILD"]
    system_command = "reduce {} -".format("-" + method.upper() if method else "")
    proc = subprocess.run(
        shlex.split(system_command),
        input=input_file.read(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    if proc.returncode != 0:
        logger.warning(
            f"Command '{system_command}' returned non-zero exit status {proc.returncode}."
        )
    output_file.write(proc.stdout)


def _protonate_with_openmm(input_file: IO, output_file: IO) -> None:
    from simtk.openmm import app

    pdb = app.PDBFile(input_file)
    modeller = app.Modeller(pdb.getTopology(), pdb.getPositions())
    modeller.addHydrogens(pH=7.0)
    app.PDBFile.writeFile(modeller.topology, modeller.positions, output_file, keepIds=True)


def add_cbetas(
    structure_df: pd.DataFrame,
    internal_coords: np.ndarray,
    use_median: bool = False,
    random_state: Optional[int] = None,
) -> pd.DataFrame:
    """Add C-beta atoms to every residue (including glycine) in the provided dataframe.

    Useful in cases where we want to e.g. calculate the distance between the C-beta atom
    of every residue in the structure.
    """
    assert (
        structure_df["residue_idx"].sort_values().values == structure_df["residue_idx"].values
    ).all()
    assert (
        structure_df["atom_idx"].sort_values().unique() == structure_df["atom_idx"].values
    ).all()
    assert structure_df["residue_idx"].unique().shape[0] == internal_coords.shape[0]

    # The following notebook shows how these values were calculated:
    # https://gitlab.com/datapkg/pdb-analysis/-/blob/457cc6423396dca97b67e73c147ad01a85f55f3e/notebooks/01_pdb_cbeta_stats.validation.ipynb
    cb_vector_stats = [
        (0.9357068986630099, 0.027465371535687625),
        (-1.1999374252852002, 0.020877555051197268),
        (-0.0300870890760455, 0.0117160484381385),
    ]
    cb_vector_rvs = [stats.laplace(*stat) for stat in cb_vector_stats]

    cbeta_dfs = []
    for (residue_idx, residue_df), residue_internal_coords in zip(
        structure_df.groupby("residue_idx"), internal_coords
    ):
        if "CB" in residue_df["atom_name"].values:
            continue
        ca_row = residue_df[residue_df["atom_name"] == "CA"].iloc[:1].copy()
        ca_coords = ca_row.iloc[0][["atom_x", "atom_y", "atom_z"]].values
        if use_median:
            cb_vector_internal = np.hstack([stat[0] for stat in cb_vector_stats])
        else:
            cb_vector_internal = np.hstack(
                [rv.rvs(random_state=random_state) for rv in cb_vector_rvs]
            )
        cb_vector_pred = cb_vector_internal @ residue_internal_coords
        cb_coords_pred = ca_coords - cb_vector_pred
        ca_row["atom_name"] = "CB"
        ca_row["atom_fullname"] = " CB "
        ca_row["atom_x"], ca_row["atom_y"], ca_row["atom_z"] = cb_coords_pred
        cbeta_dfs.append(ca_row)

    if cbeta_dfs:
        structure_df = pd.concat([structure_df] + cbeta_dfs, ignore_index=True).sort_values(
            ["atom_idx", "atom_name"]
        )
        structure_df["atom_idx"] = np.arange(len(structure_df))
        structure_df.index = np.arange(len(structure_df))
        assert (
            structure_df["residue_idx"].sort_values().values == structure_df["residue_idx"].values
        ).all()

    return structure_df
