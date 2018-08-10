from pathlib import Path

import kmbio.PDB
import pytest

from kmtools import structure_tools

TESTS_DIR = Path(__file__).absolute().parent


@pytest.mark.skip(reason="Not implemented yet!")
def test_prepare_for_modeling():
    # TODO: This is not a good test
    structure = kmbio.PDB.load(TESTS_DIR.joinpath("structures", "4dkl.cif"))
    assert structure

    structure_fm, alignment = structure_tools.prepare_for_modeling()
    assert len(alignment[0].seq) == len(alignment[1].seq)
