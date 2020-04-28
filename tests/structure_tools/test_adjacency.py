import inspect
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pytest
from kmbio import PDB
from ruamel import yaml

from kmtools import sequence_tools, structure_tools
from kmtools.structure_tools.adjacency import get_atom_distances, get_distances

TEST_DATA_DIR = Path(__file__).parent.joinpath("structures").resolve(strict=True)

with Path(__file__).with_suffix("").joinpath("data.yaml").open("rt") as fin:
    TEST_DATA = yaml.safe_load(fin)


def load_test_cases_from_file(fn):
    args = inspect.getfullargspec(fn).args
    data = TEST_DATA[fn.__name__]
    ids = list(data.keys())
    values = [tuple(d[a] for a in args) for d in data.values()]
    parametrize = pytest.mark.parametrize(", ".join(args), values, ids=ids)
    return parametrize(fn)


def test_get_atom_distances():
    # For some reason, this structure produces a segfault using MDAnalysis's `self_capped_distance`
    # with `method="nsgrid"`.
    structure_file = Path(__file__).with_suffix("").joinpath("5nleA.pdb").resolve(strict=True)
    atom_pairs_file = (
        Path(__file__).with_suffix("").joinpath("5nleA-atom-pairs.parquet").resolve(strict=True)
    )
    structure_df = PDB.load(structure_file).to_dataframe()
    pairs_df = (
        pq.read_table(atom_pairs_file)
        .to_pandas(integer_object_nulls=True)
        .sort_values(["atom_idx_1", "atom_idx_2", "distance"], ascending=True)
    )

    # Test with max_cutoff
    pairs_df_ = get_atom_distances(structure_df, max_cutoff=12).sort_values(
        ["atom_idx_1", "atom_idx_2", "distance"], ascending=True
    )
    assert np.allclose(pairs_df.values, pairs_df_.values, atol=1e-05, rtol=1e-05)

    # Test without max_cutoff
    pairs_df_ = get_atom_distances(structure_df, max_cutoff=None).sort_values(
        ["atom_idx_1", "atom_idx_2", "distance"], ascending=True
    )
    len(pairs_df) != len(pairs_df_)
    pairs_df_ = pairs_df_[pairs_df_["distance"] <= 12]
    assert np.allclose(pairs_df.values, pairs_df_.values, atol=1e-05, rtol=1e-05)


@pytest.mark.parametrize(
    "structure_name, groupby_method, distances_expected",
    [
        (
            "AE-AE.pdb",
            "residue",
            [
                [0.0, 1.28918, 6.08953, 1.44832],
                [1.28918, 0.0, 4.88483, 0.82686],
                [6.08953, 4.88483, 0.0, 1.28922],
                [1.44832, 0.82686, 1.28922, 0.0],
            ],
        ),
        (
            "AE-AE.pdb",
            "residue-backbone",
            [
                [0.0, 1.28917726, 6.08952855, 5.22890476],
                [1.28917726, 0.0, 6.46029798, 4.8154139],
                [6.08952855, 6.46029798, 0.0, 1.28921914],
                [5.22890476, 4.8154139, 1.28921914, 0.0],
            ],
        ),
        (
            "AE-AE.pdb",
            "residue-ca",
            [
                [0.0, 3.76792476, 7.40403653, 5.80728792],
                [3.76792476, 0.0, 7.94448551, 4.8340331],
                [7.40403653, 7.94448551, 0.0, 3.76698633],
                [5.80728792, 4.8340331, 3.76698633, 0.0],
            ],
        ),
        (
            "AE-AE.pdb",
            "residue-cb",
            [
                [0.0, 5.75226651, 10.37423886, 5.81188653],
                [5.75226651, 0.0, 7.99318272, 2.61213227],
                [10.37423886, 7.99318272, 0.0, 5.75206754],
                [5.81188653, 2.61213227, 5.75206754, 0.0],
            ],
        ),
    ],
)
def test_get_distances_residue(structure_name, groupby_method, distances_expected):
    structure_file = TEST_DATA_DIR.joinpath(structure_name)
    structure = PDB.load(structure_file)
    max_cutoff = np.max(distances_expected) + 0.1

    def distance_df_to_matrix(distances_df):
        distances_df = pd.concat(
            [
                distances_df,
                distances_df.rename(
                    columns={"residue_idx_1": "residue_idx_2", "residue_idx_2": "residue_idx_1"}
                ),
                pd.DataFrame(
                    {
                        "residue_idx_1": np.arange(len(distances_expected)),
                        "residue_idx_2": np.arange(len(distances_expected)),
                        "distance": 0.0,
                    }
                ),
            ],
            sort=False,
        ).sort_values(["residue_idx_1", "residue_idx_2"])
        distances = distances_df.pivot_table("distance", "residue_idx_1", "residue_idx_2").values
        return distances

    # Test with max_cutoff
    distances_df = get_distances(structure.to_dataframe(), max_cutoff, groupby=groupby_method)
    distances = distance_df_to_matrix(distances_df)
    assert np.allclose(distances, distances_expected, atol=0.01)

    # Test without max_cutoff
    distances_df = get_distances(structure.to_dataframe(), None, groupby=groupby_method)
    distances_df = distances_df[distances_df["distance"] <= max_cutoff]
    distances = distance_df_to_matrix(distances_df)
    assert np.allclose(distances, distances_expected, rtol=0.01)


@load_test_cases_from_file
def test_map_distances(
    structure_file, max_cutoff, b2a, residue_idx_1_corrected, residue_idx_2_corrected
):
    # Convert lists to arrays
    b2a = np.array(b2a)
    structure = PDB.load(Path(__file__).with_suffix("").joinpath(structure_file))
    # Calculate interactions
    distances = structure_tools.get_distances(
        structure.to_dataframe(), max_cutoff, groupby="residue"
    )
    # Map interactions to target sequence
    for i in [1, 2]:
        distances[f"residue_idx_{i}_corrected"] = distances[f"residue_idx_{i}"].apply(
            lambda idx: sequence_tools.convert_residue_index_a2b(idx, b2a)
        )
    interactions_1 = set(
        distances[[f"residue_idx_1_corrected", f"residue_idx_2_corrected"]]
        .dropna()
        .astype(int)
        .apply(tuple, axis=1)
    )
    # Get reference interactions
    interactions_2 = {
        (int(r1), int(r2)) if r1 <= r2 else (int(r2), int(r1))
        for r1, r2 in zip(residue_idx_1_corrected, residue_idx_2_corrected)
    }
    # Make sure interactions match
    assert not interactions_1 ^ interactions_2
