import io
import tempfile
from pathlib import Path
from unittest.mock import patch

from kmtools.structure_tools.fixes import _protonate_with_openmm, _protonate_with_reduce, protonate


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
