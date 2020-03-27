import tempfile
from pathlib import Path
from textwrap import dedent

import kmbio.PDB
import pytest
from Bio.AlignIO import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from kmtools import structure_tools
from kmtools.structure_tools import DomainMutation, DomainTarget

TESTS_DIR = Path(__file__).absolute().parent


@pytest.mark.parametrize(
    "structure_file, targets, num_hetatms",
    [
        (
            TESTS_DIR.joinpath("structures", "1yf4b.pdb"),
            [DomainTarget(0, "B", "CYFQNCPRG", None, None, "CAFQNCPRG")],
            4,
        ),
        (
            TESTS_DIR.joinpath("structures", "1yf4.cif"),
            [
                DomainTarget(
                    0, "A", "LQGIVSWGYGCAQKNKPGVYT-----KV", 187, 209, "LGGGVSWGYGCAQKNKPGVYTKGGGGGV"
                ),
                DomainTarget(0, "B", "CYFQNCPRG", None, None, "CYFQNCPRG"),
            ],
            35,
        ),
    ],
)
def test_prepare_for_modeling(structure_file, targets, num_hetatms):
    structure = kmbio.PDB.load(structure_file)
    structure_fm, alignment = structure_tools.prepare_for_modeling(structure, targets)
    assert len(alignment[0].seq) == len(alignment[1].seq)
    # Validate structure
    num_chains = len(targets) + (1 if num_hetatms > 0 else 0)
    assert len(list(structure_fm[0].chains)) == num_chains
    for chain, target in zip(structure_fm.chains, targets):
        seq = target.template_sequence.replace("-", "")
        assert (
            structure_tools.get_chain_sequence(
                chain, if_unknown="replace", unknown_residue_marker=""
            )
            == seq
        )
    assert len(list(list(structure_fm[0].chains)[-1].residues)) == num_hetatms
    # Validate alignment
    num_seqs = len(targets) + (1 if num_hetatms > 0 else 0)
    assert len(alignment) == 2
    assert len(alignment[0].seq.split("/")) == num_seqs
    assert len(alignment[1].seq.split("/")) == num_seqs
    if num_hetatms:
        assert alignment[0].seq.split("/")[-1] == "." * num_hetatms
        assert alignment[1].seq.split("/")[-1] == "." * num_hetatms


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
    "mutation, target_in, target_out",
    [
        (
            DomainMutation(0, "B", "C", 1, "A"),
            DomainTarget(0, "B", "CYFQNCPRG", None, None, "CAFQNCPRG"),
            DomainTarget(0, "B", "CYFQNCPRG", None, None, "AAFQNCPRG"),
        )
    ],
)
def test_mutation_to_target(mutation, target_in, target_out):
    target_out_ = structure_tools.mutation_to_target(mutation, target_in)
    assert target_out_ == target_out
