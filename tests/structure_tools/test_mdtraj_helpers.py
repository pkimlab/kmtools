import tempfile
from pathlib import Path

import mdtraj
import numpy as np

from kmtools.structure_tools.mdtraj_helpers import mdtraj_to_dataframe, mdtraj_to_pdb


def test_mdtraj_to_dataframe():
    structure_file = Path(__file__).parent.joinpath("structures", "1yf4b.pdb").resolve(strict=True)
    traj = mdtraj.load(structure_file.as_posix())
    structure_df = mdtraj_to_dataframe(traj)
    columns_dtypes_expected = [
        ("serial", np.dtype("int64")),
        ("name", np.dtype("O")),
        ("element", np.dtype("O")),
        ("resSeq", np.dtype("int64")),
        ("resName", np.dtype("O")),
        ("chainID", np.dtype("int64")),
        ("segmentID", np.dtype("O")),
        ("residue_idx", np.dtype("int64")),
        ("chain_idx", np.dtype("int64")),
        ("model_idx", np.dtype("int64")),
        ("atom_serial_number", np.dtype("int64")),
        ("atom_name", np.dtype("O")),
        ("residue_id_1", np.dtype("int64")),
        ("residue_resname", np.dtype("O")),
        ("atom_idx", np.dtype("int64")),
        ("atom_x", np.dtype("float64")),
        ("atom_y", np.dtype("float64")),
        ("atom_z", np.dtype("float64")),
    ]
    columns = set(structure_df.columns)
    assert all(c in columns for c, _ in columns_dtypes_expected)
    assert all(structure_df[c].dtype == dtype for c, dtype in columns_dtypes_expected)


def test_mdtraj_to_pdb():
    structure_file = Path(__file__).parent.joinpath("structures", "ADNA.pdb").resolve(strict=True)
    with structure_file.open("rt") as fin:
        structure_data = fin.read().split("\n")
    traj = mdtraj.load(structure_file.as_posix())
    with tempfile.NamedTemporaryFile(suffix=".pdb") as output_file:
        mdtraj_to_pdb(traj, output_file.name)
        with open(output_file.name, "rt") as fin:
            structure_file_out = fin.read().split("\n")
    print("\n".join(structure_file_out))
    assert structure_data == structure_file_out
