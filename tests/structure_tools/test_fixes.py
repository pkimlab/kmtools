import io
import tempfile
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest
from kmbio import PDB

from kmtools.structure_tools.fixes import (
    _protonate_with_openmm,
    _protonate_with_reduce,
    add_cbetas,
    protonate,
)
from kmtools.structure_tools.protein_structure_analysis import get_internal_coords

try:
    from simtk import openmm
except ImportError:
    openmm = None


@patch("kmtools.structure_tools.fixes._protonate_with_reduce")
@patch("kmtools.structure_tools.fixes._protonate_with_openmm")
def test_protonate(_mock_protonate_with_openmm, _mock_protonate_with_reduce):
    input_file = tempfile.NamedTemporaryFile(suffix=".pdb")
    output_file = tempfile.NamedTemporaryFile(suffix=".pdb")

    for mock_object, method in [
        (_mock_protonate_with_openmm, "openmm"),
        (_mock_protonate_with_reduce, "reduce"),
    ]:
        mock_object.reset_mock()
        mock_object.assert_not_called()
        protonate(input_file.name, output_file.name, method=method)
        mock_object.assert_called_once()
        assert mock_object.call_args[0][0] != input_file.name
        assert mock_object.call_args[0][1] != output_file.name

        mock_object.reset_mock()
        mock_object.assert_not_called()
        with open(input_file.name) as fin, open(output_file.name) as fout:
            protonate(fin, fout, method=method)
            mock_object.assert_called_once()
            assert mock_object.call_args[0][0] == fin
            assert mock_object.call_args[0][1] == fout


@pytest.mark.skipif(openmm is None, reason="openmm is not installed")
def test__protonate_with_openmm():
    # Test that adding hydrogens adds lines to the PDB file
    # (not a perfect test, but ¯\_(ツ)_/¯)
    output_file = io.StringIO()
    with Path(__file__).parent.joinpath("structures/1yf4b.pdb").open("rt") as input_file:
        _protonate_with_openmm(input_file, output_file)
    output_file.seek(0)
    lines = list(output_file.read().strip().split("\n"))
    assert len(lines) == 160


def test__protonate_with_reduce():
    output_file = io.StringIO()
    with Path(__file__).parent.joinpath("structures/1yf4b.pdb").open("rt") as input_file:
        _protonate_with_reduce(input_file, output_file)
    output_file.seek(0)
    lines = list(output_file.read().strip().split("\n"))
    assert len(lines) == 155


@pytest.mark.parametrize(
    "structure_name, use_median",
    [
        (structure_name, use_median)
        for structure_name in ["AAA.pdb", "ADNA.pdb", "AEA.pdb"]
        for use_median in [True, False]
    ],
)
def test_add_cbetas(structure_name, use_median):
    structure_file = Path(__file__).parent.joinpath("structures", structure_name)
    structure = PDB.load(structure_file)
    structure_df = structure.to_dataframe()
    internal_coords = get_internal_coords(structure_df)
    structure_df_nocbetas = structure_df[structure_df["atom_name"] != "CB"]
    structure_df_wcbetas = add_cbetas(
        structure_df_nocbetas, internal_coords, use_median=False, random_state=0
    )

    def reorder_cbetas(structure_df):
        structure_df_ = structure_df.copy()
        structure_df_.loc[structure_df_["atom_name"] == "CB", "atom_idx"] = structure_df_.loc[
            structure_df_["atom_name"] == "CA", "atom_idx"
        ].values
        structure_df_.sort_values(["atom_idx", "atom_name"], inplace=True)
        structure_df_["atom_idx"] = structure_df["atom_idx"].values
        structure_df_.index = structure_df.index.values
        return structure_df_

    structure_df_ = reorder_cbetas(structure_df)
    assert structure_df_.columns.equals(structure_df_wcbetas.columns)

    coord_columns = ["atom_x", "atom_y", "atom_z"]
    assert np.allclose(
        structure_df_[coord_columns].values,
        structure_df_wcbetas[coord_columns].values,
        rtol=0.1,
        atol=0.05,
    )

    columns_to_skip = ["atom_serial_number"]
    noncoord_columns = [
        c for c in structure_df_.columns if c not in coord_columns and c not in columns_to_skip
    ]
    assert structure_df_[noncoord_columns].equals(structure_df_wcbetas[noncoord_columns])
