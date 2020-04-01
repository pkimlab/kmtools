from pathlib import Path

import mdtraj
import numpy as np
import pytest
from kmbio import PDB

from kmtools.structure_tools import protein_structure_analysis as psa
from kmtools.structure_tools.geometry import compute_angles, compute_dihedrals


def test_compute_angles_random():
    data = [
        ([1.0, 0.0], [0.0, 0.0], [0.0, 1.0], np.pi / 2),
        ([1.0, 0.0], [0.0, 0.0], [10.0, 10.0], np.pi / 4),
        ([1.0, 0.0], [0.0, 0.0], [0.0, -1.0], -np.pi / 2),
        ([1.0, 0.0], [0.0, 0.0], [-1.0, 0.000_000_000_000_1], np.pi),
        ([1.0, 0.0], [0.0, 0.0], [-1.0, -0.000_000_000_000_1], -np.pi),
    ]

    u1 = np.array([d[0] for d in data])
    u2 = np.array([d[1] for d in data])
    u3 = np.array([d[2] for d in data])
    angles_ = np.array([d[3] for d in data])

    angles = compute_angles(u1, u2, u3)
    assert np.allclose(angles, angles_)


def test_compute_angles_translations():
    data = [
        ([1.0, 0.0], [0.0, 0.0], [0.0, 1.0], np.pi / 2),
        ([1.0 * 100, 0.0], [0.0, 0.0], [0.0, 1.0], np.pi / 2),
        ([1.0, 0.0], [0.0, 0.0], [0.0, 1.0 * 100], np.pi / 2),
        ([1.0 - 100, 0.0 - 100], [0.0 - 100, 0.0 - 100], [0.0 - 100, 100.0 - 100], np.pi / 2),
    ]

    u1 = np.array([d[0] for d in data])
    u2 = np.array([d[1] for d in data])
    u3 = np.array([d[2] for d in data])
    angles_ = np.array([d[3] for d in data])

    angles = compute_angles(u1, u2, u3)
    assert np.allclose(angles, angles_)


def test_compute_angles_rotations():
    data = [
        ([1.0, 0.0], [0.0, 0.0], [0.0, 1.0], np.pi / 2),
        ([0.0, 1.0], [0.0, 0.0], [-1.0, 0.0], np.pi / 2),
        ([-1.0, 0.0], [0.0, 0.0], [0.0, -1.0], np.pi / 2),
        ([0.0, -1.0], [0.0, 0.0], [1.0, 0.0], np.pi / 2),
    ]

    u1 = np.array([d[0] for d in data])
    u2 = np.array([d[1] for d in data])
    u3 = np.array([d[2] for d in data])
    angles_ = np.array([d[3] for d in data])

    angles = compute_angles(u1, u2, u3)
    assert np.allclose(angles, angles_)


def test_compute_angles_3d():
    data = [
        ([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], np.pi / 2),
        ([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], np.pi / 2),
        ([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, -1.0, 0.0], np.pi / 2),
        ([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, -1.0], np.pi / 2),
        ([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [-1.0, 0.0, 0.0], np.pi),
        ([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.05, 0.0, 0.0], 0),
    ]

    u1 = np.array([d[0] for d in data])
    u2 = np.array([d[1] for d in data])
    u3 = np.array([d[2] for d in data])
    angles_ = np.array([d[3] for d in data])

    angles = compute_angles(u1, u2, u3)
    assert np.allclose(angles, angles_)


@pytest.mark.parametrize("structure_file", ["AAA.pdb", "AEA.pdb", "ADNA.pdb"])
def test_compute_angles_mdtraj(structure_file):
    structure_file = (
        Path(__file__).parent.joinpath("structures", structure_file).resolve(strict=True)
    )

    structure_df = PDB.load(structure_file).to_dataframe()
    traj = mdtraj.load(structure_file.as_posix())

    atoms_of_interest = ["N", "CA", "CB", "C"]
    residue_atom_coords_df = psa.get_residue_atom_coords(structure_df, atoms_of_interest)

    atom_coord_columns = ["atom_x", "atom_y", "atom_z"]

    # Validate angles between pairs of interacting residues
    residue_atom_coords_df["dummy"] = 1
    residue_pair_atom_coords_df = residue_atom_coords_df.merge(
        residue_atom_coords_df, on="dummy", suffixes=("_1", "_2")
    )
    residue_pair_atom_coords_df = residue_pair_atom_coords_df[
        residue_pair_atom_coords_df["residue_idx_1"] < residue_pair_atom_coords_df["residue_idx_2"]
    ]
    del residue_atom_coords_df["dummy"]
    del residue_pair_atom_coords_df["dummy"]

    for atom_name_suffix_lst in [
        #
        [("CA", 1), ("CB", 1), ("CB", 2)],
        [("CB", 1), ("CB", 2), ("CA", 2)],
    ]:
        dihedrals = compute_angles(
            *(
                residue_pair_atom_coords_df[
                    [f"{atom_name}_{col}_{suffix}" for col in atom_coord_columns]
                ].values
                for atom_name, suffix in atom_name_suffix_lst
            )
        )
        atom_indices = residue_pair_atom_coords_df[
            [f"{atom_name}_atom_idx_{suffix}" for atom_name, suffix in atom_name_suffix_lst]
        ].values
        dihedrals_ = mdtraj.compute_angles(traj, angle_indices=atom_indices)
        assert np.allclose(dihedrals, dihedrals_)


def test_compute_dihedrals():
    data = [
        ([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 1.0, 1.0], -np.pi / 2),
        ([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 1.0, 1000.0], -np.pi / 2),
        ([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 1.0], np.pi / 2),
        ([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 100.0, 1.0], np.pi / 2),
        ([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 1.0], -np.pi / 2),
    ]

    u1 = np.array([d[0] for d in data])
    u2 = np.array([d[1] for d in data])
    u3 = np.array([d[2] for d in data])
    u4 = np.array([d[3] for d in data])
    dihedrals_ = np.array([d[4] for d in data])

    dihedrals = compute_dihedrals(u1, u2, u3, u4)
    assert np.allclose(dihedrals, dihedrals_)


@pytest.mark.parametrize("structure_file", ["AAA.pdb", "AEA.pdb", "ADNA.pdb"])
def test_compute_dihedrals_mdtraj(structure_file):
    structure_file = (
        Path(__file__).parent.joinpath("structures", structure_file).resolve(strict=True)
    )

    structure_df = PDB.load(structure_file).to_dataframe()
    traj = mdtraj.load(structure_file.as_posix())

    atoms_of_interest = ["N", "CA", "CB", "C"]
    residue_atom_coords_df = psa.get_residue_atom_coords(structure_df, atoms_of_interest)

    atom_coord_columns = ["atom_x", "atom_y", "atom_z"]

    # Validate phi angles
    atom_coords = [
        residue_atom_coords_df.iloc[:-1][[f"C_{col}" for col in atom_coord_columns]].values,
        residue_atom_coords_df.iloc[1:][[f"N_{col}" for col in atom_coord_columns]].values,
        residue_atom_coords_df.iloc[1:][[f"CA_{col}" for col in atom_coord_columns]].values,
        residue_atom_coords_df.iloc[1:][[f"C_{col}" for col in atom_coord_columns]].values,
    ]
    dihedrals = compute_dihedrals(*atom_coords)
    _, dihedrals_ = mdtraj.compute_phi(traj, periodic=False)
    assert np.allclose(dihedrals, dihedrals_)

    # Validate psi angles
    atom_coords = [
        residue_atom_coords_df.iloc[:-1][[f"N_{col}" for col in atom_coord_columns]].values,
        residue_atom_coords_df.iloc[:-1][[f"CA_{col}" for col in atom_coord_columns]].values,
        residue_atom_coords_df.iloc[:-1][[f"C_{col}" for col in atom_coord_columns]].values,
        residue_atom_coords_df.iloc[1:][[f"N_{col}" for col in atom_coord_columns]].values,
    ]
    dihedrals = compute_dihedrals(*atom_coords)
    _, dihedrals_ = mdtraj.compute_psi(traj, periodic=False)
    assert np.allclose(dihedrals, dihedrals_)

    # Validate angles between pairs of interacting residues
    residue_atom_coords_df["dummy"] = 1
    residue_pair_atom_coords_df = residue_atom_coords_df.merge(
        residue_atom_coords_df, on="dummy", suffixes=("_1", "_2")
    )
    residue_pair_atom_coords_df = residue_pair_atom_coords_df[
        residue_pair_atom_coords_df["residue_idx_1"] < residue_pair_atom_coords_df["residue_idx_2"]
    ]
    del residue_atom_coords_df["dummy"]
    del residue_pair_atom_coords_df["dummy"]

    for atom_name_suffix_lst in [
        #
        [("N", 1), ("CA", 1), ("CB", 1), ("CB", 2)],
        [("CA", 1), ("CB", 1), ("CB", 2), ("CA", 2)],
        [("CB", 1), ("CB", 2), ("CA", 2), ("N", 2)],
    ]:
        dihedrals = compute_dihedrals(
            *(
                residue_pair_atom_coords_df[
                    [f"{atom_name}_{col}_{suffix}" for col in atom_coord_columns]
                ].values
                for atom_name, suffix in atom_name_suffix_lst
            )
        )
        atom_indices = residue_pair_atom_coords_df[
            [f"{atom_name}_atom_idx_{suffix}" for atom_name, suffix in atom_name_suffix_lst]
        ].values
        dihedrals_ = mdtraj.compute_dihedrals(traj, indices=atom_indices)
        assert np.allclose(dihedrals, dihedrals_)
