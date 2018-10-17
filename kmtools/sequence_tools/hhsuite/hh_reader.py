#!/usr/bin/env python
"""
Parser for hhr result files created with hhblits|hhsearch|hhalign -o <hhr_file>

Source: https://github.com/soedinglab/hh-suite/blob/master/scripts/hh_reader.py
"""
__authors__ = (
    "Markus Meier (markus.meier@mpibpc.mpg.de), Alexey Strokach (alex.strokach@utoronto.ca)"
)
__version__ = "2.0"
__license__ = "GPL-3"

import warnings
from typing import Iterable, List, Optional

import attr


@attr.s
class HHRAlignment:
    # Identifiers
    query_id: Optional[str] = attr.ib(default=None)
    query_length: Optional[int] = attr.ib(default=None)
    query_neff: Optional[float] = attr.ib(default=None)
    template_id: Optional[str] = attr.ib(default=None)
    template_length: int = attr.ib(default=0)
    template_info: Optional[str] = attr.ib(default=None)
    template_neff: Optional[float] = attr.ib(default=None)
    # Sequence
    query_ali: str = attr.ib(default="")
    query_ss_pred: str = attr.ib(default="")
    template_ali: str = attr.ib(default="")
    template_ss_pred: str = attr.ib(default="")
    template_ss_dssp: str = attr.ib(default="")
    ali_confidence: str = attr.ib(default="")
    # Alignment definition
    query_start: Optional[int] = attr.ib(default=None)
    query_end: Optional[int] = attr.ib(default=None)
    template_start: Optional[int] = attr.ib(default=None)
    template_end: Optional[int] = attr.ib(default=None)
    # Match metrics
    probability: Optional[float] = attr.ib(default=None)
    evalue: Optional[float] = attr.ib(default=None)
    score: Optional[float] = attr.ib(default=None)
    aligned_cols: Optional[int] = attr.ib(default=None)
    identity: Optional[float] = attr.ib(default=None)
    similarity: Optional[float] = attr.ib(default=None)
    sum_probs: Optional[float] = attr.ib(default=None)

    def validate(self):
        assert self.query_id is not None
        assert self.query_length > 0
        assert self.template_id is not None
        assert self.template_length > 0
        assert len(self.query_ali) == len(self.template_ali) == len(self.ali_confidence)
        if self.template_ss_pred:
            if not self.query_ss_pred:
                warnings.warn("'query_ss_pred' is empty even though 'template_ss_pred' is not!")
            else:
                assert (
                    len(self.query_ali)
                    == len(self.query_ss_pred)
                    == len(self.template_ss_pred)
                    == len(self.template_ss_dssp)
                )


class HHRFormatError(Exception):
    def __init__(self, value):
        self.value = "ERROR: " + value

    def __str__(self):
        return repr(self.value)


def parse_hhr_data(lines: Iterable[str]) -> List[HHRAlignment]:
    """Parse an *.hhr file into a list of :any:`HHRAlignment` objects."""
    query_id = None
    query_length = None
    query_neff = None
    alignment = HHRAlignment(query_id, query_length, query_neff)
    results = []
    skipped_ali_tags = ["Consensus"]
    is_alignment_section = False
    for line in lines:
        if line.startswith("Query"):
            assert query_id is None and alignment.query_id is None
            query_id = line.split()[1]
            alignment.query_id = query_id
        elif line.startswith("Match_columns"):
            assert query_length is None and alignment.query_length is None
            query_length = int(line.split()[1])
            alignment.query_length = query_length
        elif line.startswith("Neff"):
            assert query_neff is None and alignment.query_neff is None
            query_neff = float(line.split()[1])
            alignment.query_neff = query_neff
        elif is_alignment_section and (line.startswith("No") or line.startswith("Done!")):
            if alignment.query_start is not None:
                alignment.validate()
                results.append(alignment)
            alignment = HHRAlignment(query_id, query_length, query_neff)
        elif line.startswith("Probab"):
            assert alignment.probability is None
            tokens = line.split()
            alignment.probability = float(tokens[0].split("=")[1])
            alignment.evalue = float(tokens[1].split("=")[1])
            alignment.score = float(tokens[2].split("=")[1])
            alignment.aligned_cols = int(tokens[3].split("=")[1])
            alignment.identity = float(tokens[4].split("=")[1].replace("%", "")) / 100.0
            alignment.similarity = float(tokens[5].split("=")[1])
            alignment.sum_probs = float(tokens[6].split("=")[1])
            if len(tokens) > 7:
                alignment.template_neff = float(tokens[7].split("=")[1])
            continue
        elif line.startswith(">"):
            is_alignment_section = True
            assert alignment.template_id is None and alignment.template_info is None
            alignment.template_id = line[1:].split()[0]
            alignment.template_info = line
        elif line.startswith("Q"):
            tokens = line.split()
            if tokens[1] in skipped_ali_tags:
                continue
            elif tokens[1] == "ss_pred":
                alignment.query_ss_pred += tokens[2]
                continue

            alignment.query_ali += tokens[3]

            try:
                token_2 = int(tokens[2].replace("(", "").replace(")", ""))
            except Exception:
                raise HHRFormatError(
                    ("Converting failure of start index ({}) " "of query alignment").format(
                        tokens[2]
                    )
                )
            if alignment.query_start is None:
                alignment.query_start = token_2
            else:
                alignment.query_start = min(alignment.query_start, token_2)

            try:
                token_4 = int(tokens[4].replace("(", "").replace(")", ""))
            except Exception:
                raise HHRFormatError(
                    ("Converting failure of end index ({}) " "of query alignment").format(tokens[4])
                )
            if alignment.query_end is None:
                alignment.query_end = token_4
            else:
                alignment.query_end = max(alignment.query_end, token_4)
        elif line.startswith("T"):
            tokens = line.split()
            if tokens[1] in skipped_ali_tags:
                continue
            elif tokens[1] == "ss_pred":
                alignment.template_ss_pred += tokens[2]
                continue
            elif tokens[1] == "ss_dssp":
                alignment.template_ss_dssp += tokens[2]
                continue

            alignment.template_ali += tokens[3]

            try:
                token_2 = int(tokens[2].replace("(", "").replace(")", ""))
            except Exception:
                raise HHRFormatError(
                    ("Converting failure of start index ({}) " "of template alignment").format(
                        tokens[2]
                    )
                )
            if alignment.template_start is None:
                alignment.template_start = token_2
            else:
                alignment.template_start = min(alignment.template_start, token_2)

            try:
                token_4 = int(tokens[4].replace("(", "").replace(")", ""))
            except Exception:
                raise HHRFormatError(
                    ("Converting failure of end index ({}) " "of template alignment").format(
                        tokens[4]
                    )
                )
            if alignment.template_end is None:
                alignment.template_end = token_4
            else:
                alignment.template_end = max(alignment.template_end, token_4)

            alignment.template_length += token_4 - token_2 + 1

        elif line.startswith("Confidence"):
            alignment.ali_confidence += line.split("Confidence            ")[-1].strip("\n")

    if alignment.template_id is not None and alignment.query_start is not None:
        results.append(alignment)

    return results
