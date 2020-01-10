import io
from pathlib import Path

import pytest

from kmtools.structure_tools.fixes import _protonate_with_openmm, _protonate_with_reduce, protonate


@pytest.mark.skip("To do.")
def test_protonate():
    assert protonate


def test__protonate_with_openmm():
    # Test that adding hydrogens adds lines to the PDB file
    # (not a perfect test, but ¯\_(ツ)_/¯)
    output_file = io.StringIO()
    with Path(__file__).parent.joinpath("structures/1yf4b.pdb").open("rt") as input_file:
        _protonate_with_openmm(input_file, output_file)
    output_file.seek(0)
    lines = list(output_file.read().strip().split("\n"))
    assert len(lines) == 160


@pytest.mark.skip("Not implemented.")
def test__protonate_with_reduce():
    output_file = io.StringIO()
    with Path(__file__).parent.joinpath("structures/1yf4b.pdb").open("rt") as input_file:
        _protonate_with_reduce(input_file, output_file)
    output_file.seek(0)
    lines = list(output_file.read().strip().split("\n"))
    assert len(lines) == 160
