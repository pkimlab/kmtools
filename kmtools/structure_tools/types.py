import re
from typing import NamedTuple, Optional

MUTATION_REGEX = re.compile(
    "([^-_\.]*[-_]{1})?([GVALICMFWPDESTYQNKRH])([0-9]+)([GVALICMFWPDESTYQNKRH]?)"
)


class DomainMutation(NamedTuple):
    #: ID of the model to mutate (usually 0)
    model_id: int
    #: ID of the chain to mutate
    chain_id: str
    #: Original residue (usually the wild-type residue)
    residue_wt: str
    #: Position of the residue to mutate (1-based)
    residue_id: int
    #: Mutant residue
    residue_mut: str

    def __str__(self) -> str:
        return f"{self.chain_id}_{self.residue_wt}{self.residue_id}{self.residue_mut}"

    @classmethod
    def from_string(cls, mutation: str, _mutation_regex=MUTATION_REGEX) -> "DomainMutation":
        matches = _mutation_regex.findall(mutation)
        assert len(matches) == 1
        chain_id, residue_wt, residue_id, residue_mut = matches[0]
        chain_id = chain_id.strip("-_")
        residue_id = int(residue_id)
        return cls(0, chain_id, residue_wt, residue_id, residue_mut)


class DomainTarget(NamedTuple):
    #: ID of the template model
    model_id: int
    #: ID of the template chain
    chain_id: str
    #: Sequence of the template protein, with gaps in the alignment indicated with '-'
    template_sequence: str
    #: Start position of the template chain (1-based)
    template_start: Optional[int]
    #: End position of the template chain (inclusive)
    template_end: Optional[int]
    #: Sequence of the target protein, with gaps in the alignmend indicated with '-'
    target_sequence: str
