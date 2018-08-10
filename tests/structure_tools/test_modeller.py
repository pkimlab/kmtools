import tempfile
from pathlib import Path
from textwrap import dedent

import kmbio.PDB
import pytest
from Bio.AlignIO import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from kmtools import structure_tools
from kmtools.structure_tools import DomainTarget, run_modeller

TESTS_DIR = Path(__file__).absolute().parent


def test_write_pir_alignment():
    seqrec_1 = SeqRecord(Seq("MVTGITIM"), id="test")
    seqrec_2 = SeqRecord(Seq("MVTAITIM"), id="4dkl")
    alignment = MultipleSeqAlignment([seqrec_1, seqrec_2])
    with tempfile.NamedTemporaryFile() as file_:
        structure_tools.write_pir_alignment(alignment, file_.name)
        with open(file_.name, "rt") as fin:
            assert fin.read() == dedent(
                """\
                >P1;test
                sequence:test:.:.:.:.::::
                MVTGITIM*

                >P1;4dkl
                structure:4dkl:.:.:.:.::::
                MVTAITIM*
                """
            )


@pytest.mark.parametrize(
    "structure_file, t2t_align",
    [
        (
            TESTS_DIR.joinpath("structures", "1yf4b.pdb"),
            DomainTarget("1yf4b", 0, "B", "CAFQNCPRG", "CYFQNCPRG"),
        ),
        # TODO: Add test for multiple chains
    ],
)
def test_run_modeller(structure_file: Path, t2t_align: DomainTarget):
    structure = kmbio.PDB.load(structure_file)
    structure_fm, alignment = structure_tools.prepare_for_modeling(structure, [t2t_align])
    with tempfile.TemporaryDirectory() as temp_dir:
        results = run_modeller(structure_fm, alignment, temp_dir)
        structure_bm = kmbio.PDB.load(Path(temp_dir).joinpath(results["name"]))
    assert [m.id for m in structure_bm] == [0]
    assert [c.id for c in structure_bm[0]] == [" "]
    assert structure_tools.get_chain_sequence(structure_bm[0][" "]) == t2t_align.target_sequence
    assert "Normalized DOPE score" in results
    assert "GA341 score" in results
