import tempfile
from pathlib import Path

import mdtraj
import numpy as np
import pandas as pd
import pytest

from kmtools import structure_tools
from kmtools.structure_tools.protein_structure_analysis import (
    _normalize_coords,
    calculate_backbone_angles,
    calculate_backbone_dihedrals,
    calculate_hydrogen_bonds,
    calculate_omega,
    calculate_phi,
    calculate_psi,
    calculate_sasa,
    get_internal_coords,
    get_rotations,
    get_translations,
    map_translations_to_internal_coords,
)

try:
    import torch
    import torch.nn.functional as F

    TORCH_IS_AVAILABLE = True
except ImportError:
    TORCH_IS_AVAILABLE = False

TEST_DATA_PATH = Path(__file__).parent.joinpath("structures").resolve(strict=True)


# === Features describing individual residues ===


def test_calculate_sasa():
    traj = mdtraj.load(TEST_DATA_PATH.joinpath("AEA.pdb").as_posix())
    sasa = calculate_sasa(traj)
    sasa_expected = [
        (1.594_771_6, 1.477_324_338_685_919_6),
        (1.826_687_1, 1.060_485_978_258_019_4),
        (1.476_028_2, 1.367_325_802_653_296_4),
    ]
    assert np.allclose(sasa, sasa_expected)
    assert sasa.dtype == np.double


def test_calculate_phi():
    traj = mdtraj.load(TEST_DATA_PATH.joinpath("AEA.pdb").as_posix())
    phi = calculate_phi(traj)
    phi_expected = [np.nan, -2.425_934_8, -2.426_028]
    assert np.allclose(phi, phi_expected, equal_nan=True)
    assert phi.dtype == np.double


def test_calculate_psi():
    traj = mdtraj.load(TEST_DATA_PATH.joinpath("AEA.pdb").as_posix())
    psi = calculate_psi(traj)
    psi_expected = [2.356_421, 2.356_195_2, np.nan]
    assert np.allclose(psi, psi_expected, equal_nan=True)
    assert psi.dtype == np.double


def test_calculate_omega():
    traj = mdtraj.load(TEST_DATA_PATH.joinpath("AEA.pdb").as_posix())
    omega = calculate_omega(traj)
    omega_expected = [(np.nan, -3.141_148_3), (-3.141_148_3, 3.141_450_4), (3.141_450_4, np.nan)]
    assert np.allclose(omega, omega_expected, equal_nan=True)
    assert omega.dtype == np.double


def test_calculate_backbone_angles():
    traj = mdtraj.load(TEST_DATA_PATH.joinpath("AEA.pdb").as_posix())
    backbone_angles = calculate_backbone_angles(traj)
    backbone_angles_expected = [(np.nan), (2.287_743_330_001_831), (np.nan)]
    assert np.allclose(backbone_angles, backbone_angles_expected, equal_nan=True)
    assert backbone_angles.dtype == np.double


def test_calculate_backbone_dihedrals_0():
    traj = mdtraj.load(TEST_DATA_PATH.joinpath("AEA.pdb").as_posix())
    backbone_dihedrals = calculate_backbone_dihedrals(traj)
    backbone_dihedrals_expected = [(np.nan, np.nan), (np.nan, np.nan), (np.nan, np.nan)]
    assert np.allclose(backbone_dihedrals, backbone_dihedrals_expected, equal_nan=True)
    assert backbone_dihedrals.dtype == np.double


def test_calculate_backbone_dihedrals_1():
    traj = mdtraj.load(TEST_DATA_PATH.joinpath("ADNA.pdb").as_posix())
    backbone_dihedrals = calculate_backbone_dihedrals(traj)
    backbone_dihedrals_expected = [
        (np.nan, np.nan),
        (np.nan, 3.078_418_970_108_032_2),
        (3.078_418_970_108_032_2, np.nan),
        (np.nan, np.nan),
    ]
    assert np.allclose(backbone_dihedrals, backbone_dihedrals_expected, equal_nan=True)
    assert backbone_dihedrals.dtype == np.double


# === Features describing interactions between residues ===


@pytest.mark.skipif(not TORCH_IS_AVAILABLE, reason="This test requires PyTorch.")
def test__normalize_coords_0():
    coords = np.array(
        [
            [[4.0, 1.0, -2.0], [7.0, -14.0, 7.0], [21.0, 42.0, 63.0]],
            [[0.0, 1.0, -9.0], [64.0, -63.0, -7.0], [574.0, 576.0, 64.0]],
        ]
    )
    coords_normed = _normalize_coords(coords)
    coords_normed_expected = F.normalize(torch.from_numpy(coords), p=2, dim=-1).data.numpy()
    assert np.allclose(coords_normed, coords_normed_expected)


@pytest.mark.skipif(not TORCH_IS_AVAILABLE, reason="This test requires PyTorch.")
def test__normalize_coords_1():
    coords = np.array(
        [
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]],
            [[0.8, 0.2, 0.0], [0.0, 0.0, 1.4], [-0.28, 1.12, 0.0]],
        ]
    )
    coords_normed = _normalize_coords(coords)
    coords_normed_expected = F.normalize(torch.from_numpy(coords), p=2, dim=-1).data.numpy()
    assert np.allclose(coords_normed, coords_normed_expected)


def test_get_internal_coords():
    input_data = [
        #
        (0, "N", 0.0, 0.0, 0.0),
        (0, "CA", 1.0, 2.0, 3.0),
        (0, "C", 6.0, 5.0, 4.0),
        (0, "H", 10.0, 20.0, 20.0),
        (0, "H", 40.0, 30.0, 30.0),
        (0, "H", 50.0, 60.0, 60.0),
        (1, "N", 7.0, 8.0, 10.0),
        (1, "CA", 14.0, 14.0, 20.0),
        (1, "C", 21.0, 21.0, 21.0),
    ]
    input_df = pd.DataFrame(
        input_data, columns=["residue_idx", "atom_name", "atom_x", "atom_y", "atom_z"]
    )
    internal_coords = get_internal_coords(input_df)
    internal_coords_normed_expected = np.array(
        [
            [
                [0.872_871_560_943_969_6, 0.218_217_890_235_992_4, -0.436_435_780_471_984_8],
                [0.408_248_290_463_862_96, -0.816_496_580_927_725_9, 0.408_248_290_463_862_96],
                [0.267_261_241_912_424_4, 0.534_522_483_824_848_8, 0.801_783_725_737_273_2],
            ],
            [
                [0.0, 0.110_431_526_074_846_53, -0.993_883_734_673_618_8],
                [0.710_497_366_125_577_3, -0.699_395_844_779_865_3, -0.077_710_649_419_985_03],
                [0.703_699_859_832_739_4, 0.706_151_775_720_658_3, 0.078_461_308_413_406_48],
            ],
        ]
    )
    assert np.allclose(internal_coords, internal_coords_normed_expected)


def test_get_translations():
    input_data = [
        (0, "CA", 0.0, 0.0, 0.0),
        (1, "CA", 10.0, 15.0, 10.0),
        (2, "CA", -10.0, 0.0, 10.0),
    ]
    input_df = pd.DataFrame(
        input_data, columns=["residue_idx", "atom_name", "atom_x", "atom_y", "atom_z"]
    )
    translations = get_translations(input_df)
    translations_expected = [
        [[0, 0, 0], [10, 15, 10], [-10, 0, 10]],
        [[-10, -15, -10], [0, 0, 0], [-20, -15, 0]],
        [[10, 0, -10], [20, 15, 0], [0, 0, 0]],
    ]
    assert np.allclose(translations, translations_expected)
    assert translations.dtype == np.double


def test_map_translations_to_internal_coords_0():
    # Test that translations do not change when mapping to Cartesian coordinate system.
    N = 100
    translations = np.random.uniform(-10, 10, size=3 * N ** 2).reshape(N, N, 3)
    translations[np.arange(N), np.arange(N), :] = 0
    internal_coords = np.array(
        [[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]] for _ in range(N)]
    )
    translations_internal = map_translations_to_internal_coords(translations, internal_coords)
    assert np.allclose(translations, translations_internal)
    assert translations_internal.dtype == np.double


def test_map_translations_to_internal_coords_1():
    translations = np.array(
        [
            [[0, 0, 0], [10, 15, 10], [-10, 0, 10]],
            [[-10, -15, -10], [0, 0, 0], [-20, -15, 0]],
            [[10, 0, -10], [20, 15, 0], [0, 0, 0]],
        ]
    )
    internal_coords = np.array(
        [
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]],
            [[0.8, 0.2, 0.0], [0.0, 0.0, 1.4], [-0.28, 1.12, 0.0]],
        ]
    )
    translations_internal = map_translations_to_internal_coords(translations, internal_coords)

    ti_0_0 = internal_coords[0] @ translations[0, 0].reshape(-1, 1)
    assert np.allclose(ti_0_0, translations_internal[0, 0])

    ti_1_2 = translations[1, 2] @ internal_coords[1]
    assert np.allclose(ti_1_2, translations_internal[1, 2])

    ti_2_1 = translations[2, 1] @ internal_coords[2]
    assert np.allclose(ti_2_1, translations_internal[2, 1])

    ti_2_0 = translations[2, 0] @ internal_coords[2]
    assert np.allclose(ti_2_0, translations_internal[2, 0])


def test_get_rotations_0():
    # Test to make sure we are getting the same results as with PyTorch.
    internal_coords = _normalize_coords(
        np.array(
            [
                [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]],
                [[0.8, 0.2, 0.0], [0.0, 0.0, 1.4], [-0.28, 1.12, 0.0]],
            ]
        )
    )
    rotations = get_rotations(internal_coords)
    rotations_expected = [
        [
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]],
            [
                [0.970_142_500_145_331_8, 0.242_535_625_036_332_94, 0.0],
                [0.0, 0.0, 1.0],
                [-0.242_535_625_036_333, 0.970_142_500_145_332, 0.0],
            ],
        ],
        [
            [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [
                [-0.242_535_625_036_333, 0.970_142_500_145_332, 0.0],
                [0.970_142_500_145_331_8, 0.242_535_625_036_332_94, 0.0],
                [0.0, 0.0, 1.0],
            ],
        ],
        [
            [
                [0.970_142_500_145_331_8, 0.0, -0.242_535_625_036_333],
                [0.242_535_625_036_332_94, 0.0, 0.970_142_500_145_332],
                [0.0, 1.0, 0.0],
            ],
            [
                [-0.242_535_625_036_333, 0.970_142_500_145_331_8, 0.0],
                [0.970_142_500_145_332, 0.242_535_625_036_332_94, 0.0],
                [0.0, 0.0, 1.0],
            ],
            [
                [0.999_999_999_999_999_8, -1.040_476_753_571_359e-16, 0.0],
                [-1.040_476_753_571_359e-16, 1.000_000_000_000_000_2, 0.0],
                [0.0, 0.0, 1.0],
            ],
        ],
    ]
    assert np.allclose(rotations, rotations_expected)
    assert rotations.dtype == np.double


def test_get_rotations_1():
    internal_coords = _normalize_coords(
        np.array(
            [
                [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]],
                [[0.8, 0.2, 0.0], [0.0, 0.0, 1.4], [-0.28, 1.12, 0.0]],
            ]
        )
    )
    rotations = get_rotations(internal_coords)

    # Make sure rotations are reversible
    test_array = np.array([1.0, 2.0, -3.456]).reshape(-1, 1)

    test_array_expected = rotations[2][0] @ (rotations[0][2] @ test_array)
    assert np.allclose(test_array, test_array_expected)

    test_array_expected = rotations[0][2] @ (rotations[2][0] @ test_array)
    assert np.allclose(test_array, test_array_expected)

    test_array_expected = rotations[2][1] @ (rotations[1][2] @ test_array)
    assert np.allclose(test_array, test_array_expected)

    # Make sure rotations make sense for demo
    test_array = np.array([1.0, 3.0, -5.0]).reshape(-1, 1)

    test_array_out = rotations[0][1] @ test_array
    test_array_out_expected = np.array([3.0, -5.0, 1.0]).reshape(-1, 1)
    assert np.allclose(test_array_out, test_array_out_expected)

    test_array_out = rotations[1][0] @ test_array
    test_array_out_expected = np.array([-5.0, 1.0, 3.0]).reshape(-1, 1)
    assert np.allclose(test_array_out, test_array_out_expected)


def test_calculate_hydrogen_bonds_0():
    input_file = TEST_DATA_PATH.joinpath("3fndA02.pdb")
    traj = mdtraj.load(input_file.as_posix())
    hbonds = calculate_hydrogen_bonds(traj)
    assert hbonds.columns.values.tolist() == ["residue_index_1", "residue_index_2"]
    # TODO: Need some consistency with dtypes of empty arrays?
    # assert hbonds.dtypes.values.tolist() == [np.int, np.int]
    assert hbonds.empty


def test_calculate_hydrogen_bonds_1():
    input_file = TEST_DATA_PATH.joinpath("3fndA02.pdb")
    with tempfile.NamedTemporaryFile(suffix=".pdb") as tmp_file:
        with open(tmp_file.name, "wt") as fout:
            structure_tools.protonate(input_file, fout, method="reduce")
        traj = mdtraj.load(tmp_file.name)
    hbonds = calculate_hydrogen_bonds(traj)
    hbonds_expected = [
        (0, 21),
        (2, 19),
        (3, 42),
        (4, 6),
        (4, 17),
        (4, 41),
        (5, 40),
        (6, 9),
        (6, 40),
        (7, 7),
        (7, 10),
        (9, 12),
        (14, 16),
        (14, 17),
        (14, 18),
        (20, 24),
        (21, 25),
        (22, 26),
        (23, 27),
        (24, 28),
        (25, 29),
        (29, 32),
        (30, 33),
        (32, 43),
        (36, 43),
    ]
    assert hbonds.columns.values.tolist() == ["residue_index_1", "residue_index_2"]
    assert hbonds.dtypes.values.tolist() == [np.int, np.int]
    assert np.allclose(hbonds.values, hbonds_expected)
