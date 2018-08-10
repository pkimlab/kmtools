from typing import List, NamedTuple, Optional, Tuple

from Bio.AlignIO import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from kmbio.PDB import Chain, Model, Structure

from kmtools.structure_tools import CHAIN_IDS

from .structure_parser import copy_hetatm_chain, get_chain_sequence


class DomainTarget(NamedTuple):
    structure_id: str
    model_id: int
    chain_id: str
    #: Sequence of the target protein, with gaps in the alignmend indicated with '-'
    target_sequence: str
    #: Sequence of the template protein, with gaps in the alignment indicated with '-'
    template_sequence: str
    # Domain start and stop
    target_start: Optional[int] = None
    target_end: Optional[int] = None
    template_start: Optional[int] = None
    template_end: Optional[int] = None


def prepare_for_modeling(
    structure: Structure, targets: List[DomainTarget]
) -> Tuple[Structure, MultipleSeqAlignment]:
    """Return a structure and an alignment that can be provided as input to Modeller."""
    template_structure = Structure(structure.id, [Model(0)])
    # Add amino acid chains
    for chain_id, target in zip(CHAIN_IDS, targets):
        residues = list(structure[target.model_id][target.chain_id])
        if target.template_start is not None and target.template_end is not None:
            residues = residues[target.template_start - 1 : target.template_end]
        chain = Chain(chain_id, residues)
        chain_sequence = get_chain_sequence(chain)
        chain_sequence_expected = target.template_sequence.replace("-", "")
        assert chain_sequence == chain_sequence_expected, (chain_sequence, chain_sequence_expected)
        template_structure[0].add(chain)
    # Add hetatm chain
    hetatm_chain_final = Chain(CHAIN_IDS[CHAIN_IDS.index(chain_id) + 1])
    residue_idx = 0
    for chain in structure.chains:
        if not _chain_is_hetatm(chain):
            continue
        hetatm_chain = copy_hetatm_chain(template_structure, chain, r_cutoff=5)
        for residue in list(hetatm_chain.residues):
            residue.id = (residue.id[0], residue_idx + 1, residue.id[2])
            residue_idx += 1
            hetatm_chain_final.add(residue)
    if list(hetatm_chain_final.residues):
        template_structure[0].add(hetatm_chain_final)
    # Generate alignment
    target_seq = "/".join([target.target_sequence for target in targets])
    template_seq = "/".join([target.template_sequence for target in targets])
    num_hetatm_residues = len(list(hetatm_chain_final.residues))
    if num_hetatm_residues:
        target_seq += "/" + ("." * num_hetatm_residues)
        template_seq += "/" + ("." * num_hetatm_residues)
    alignment = MultipleSeqAlignment(
        [SeqRecord(Seq(target_seq), "target"), SeqRecord(Seq(template_seq), template_structure.id)]
    )
    return template_structure, alignment


def _chain_is_hetatm(chain: Chain) -> bool:
    """Return `True` if `chain` contains predominantly heteroatoms."""
    fraction_hetatm = sum(bool(residue.id[0].strip()) for residue in chain) / len(chain)
    return fraction_hetatm > 0.80
