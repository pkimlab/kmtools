import json
import os.path as op
import string
from typing import Dict, List, Set

import pandas as pd

DATA_DIR = op.join(op.dirname(op.abspath(__file__)), "data")

#: Mapping from 1-letter amino acid codes to 3-letter amino acid codes
A_DICT: Dict[str, str] = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "E": "GLU",
    "Q": "GLN",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
    "U": "SEC",
    "O": "PYL",
    "B": "ASX",
    "Z": "GLX",
    "J": "XLE",
    "X": "XAA",
    "*": "TER",
}

#: Mapping from 3-letter amino acid codes to 1-letter amino acid codes
AAA_DICT: Dict[str, str] = {
    **{value: key for key, value in A_DICT.items()},
    # "UNK": "X",
}

METHYLATED_LYSINES = ["MLZ", "MLY", "M3L"]
LYSINE_ATOMS = ["N", "CA", "CB", "CG", "CD", "CE", "NZ", "C", "O"]

#: Heteroatoms that are very common in PDB structures
COMMON_HETATMS = ["NH2", "HOH", "SO4"]

#: One- or two-letter identifies that should be assigned to new protein structures
CHAIN_IDS: List[str] = list(string.ascii_uppercase + string.digits + string.ascii_lowercase)
CHAIN_IDS += [(a + b) for a in CHAIN_IDS for b in CHAIN_IDS if a != b]

# Standard accessibilities for a ALA-X-ALA tripeptide (obtained from NACCESS)
with open(op.join(DATA_DIR, "standard_sasa.txt"), "rt") as fin:
    _standard_sasa = fin.read().strip()
STANDARD_SASA_ALL = [[l.strip() for l in line.split()] for line in _standard_sasa.split("\n")[1:]]
STANDARD_SASA = {x[3]: float(x[4]) for x in STANDARD_SASA_ALL}

#: Map modified residue ids to canonical residue ids
RESIDUE_MAPPING_TO_CANONICAL: Dict[str, str]
with open(op.join(DATA_DIR, "aaa_mapping_to_canonical.json"), "rt") as fin:
    RESIDUE_MAPPING_TO_CANONICAL = json.load(fin)
    RESIDUE_MAPPING_TO_CANONICAL.update({v: v for v in RESIDUE_MAPPING_TO_CANONICAL.values()})

#
RNA_MAPPING_TO_CANONICAL: Dict[str, str]
with open(op.join(DATA_DIR, "rna_mapping_to_canonical.json"), "rt") as fin:
    RNA_MAPPING_TO_CANONICAL = json.load(fin)
    RNA_MAPPING_TO_CANONICAL.update({v: v for v in RNA_MAPPING_TO_CANONICAL.values()})

DNA_MAPPING_TO_CANONICAL: Dict[str, str]
with open(op.join(DATA_DIR, "dna_mapping_to_canonical.json"), "rt") as fin:
    DNA_MAPPING_TO_CANONICAL = json.load(fin)
    DNA_MAPPING_TO_CANONICAL.update({v: v for v in DNA_MAPPING_TO_CANONICAL.values()})

RESIDUE_ATOM_NAMES: Dict[str, Set[str]] = {}
for key, group in pd.read_csv(
    op.join(DATA_DIR, "atom_nom.tbl"),
    sep="\t",
    comment="#",
    names=["AA", "BMRB", "SC", "PDB", "UCSF", "MSI", "XPLOR", "SYBYL", "MIDAS", "DIANA"],
).groupby(level=0):
    aaa = A_DICT[key]
    RESIDUE_ATOM_NAMES[aaa] = set(group["PDB"].values.tolist())
