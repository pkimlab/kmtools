import inspect
from pathlib import Path

import numpy as np
import pytest
from kmbio import PDB
from ruamel import yaml

from kmtools import sequence_tools, structure_tools

with Path(__file__).with_suffix("").joinpath("data.yaml").open("rt") as fin:
    TEST_DATA = yaml.safe_load(fin)


def load_test_cases_from_file(fn):
    args = inspect.getfullargspec(fn).args
    data = TEST_DATA[fn.__name__]
    ids = list(data.keys())
    values = [tuple(d[a] for a in args) for d in data.values()]
    parametrize = pytest.mark.parametrize(", ".join(args), values, ids=ids)
    return parametrize(fn)


@pytest.mark.parametrize("pdb_id", ["4dkl"])
def test_get_distances(pdb_id):
    structure = PDB.load(f"rcsb://{pdb_id}.cif")
    dd = structure_tools.DomainDef(0, "A", 10, 40)
    domain = structure_tools.extract_domain(structure, [dd])
    # BioPython's KDTree method
    interactions = structure_tools.get_interactions(domain, 5, add_reverse=False)
    interactions_core, interactions_interface = structure_tools.process_interactions(interactions)
    assert (interactions_core["residue_idx_1"] <= interactions_core["residue_idx_2"]).all()
    interactions_set_1 = set(
        interactions_core[["residue_idx_1", "residue_idx_2"]].apply(tuple, axis=1)
    )
    # MDAnalysis method (with distances)
    pairs_df = structure_tools.get_distances(domain, 5, groupby="residue")
    interactions_set_2 = set(pairs_df[["residue_idx_1", "residue_idx_2"]].apply(tuple, axis=1))
    #
    assert not interactions_set_1 ^ interactions_set_2


@load_test_cases_from_file
def test_map_distances(
    structure_file, max_cutoff, min_cutoff, b2a, residue_idx_1_corrected, residue_idx_2_corrected
):
    # Convert lists to arrays
    b2a = np.array(b2a)
    structure = PDB.load(Path(__file__).with_suffix("").joinpath(structure_file))
    # Calculate interactions
    distances = structure_tools.get_distances(structure, max_cutoff, min_cutoff, groupby="residue")
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
