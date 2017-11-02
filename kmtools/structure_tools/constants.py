import json
import os.path as op
import string

DATA_DIR = op.join(op.dirname(op.abspath(__file__)), 'data')
A_DICT = {
    'A': 'ALA',
    'R': 'ARG',
    'N': 'ASN',
    'D': 'ASP',
    'C': 'CYS',
    'E': 'GLU',
    'Q': 'GLN',
    'G': 'GLY',
    'H': 'HIS',
    'I': 'ILE',
    'L': 'LEU',
    'K': 'LYS',
    'M': 'MET',
    'F': 'PHE',
    'P': 'PRO',
    'S': 'SER',
    'T': 'THR',
    'W': 'TRP',
    'Y': 'TYR',
    'V': 'VAL',
    'U': 'SEC',
    'O': 'PYL',
    'B': 'ASX',
    'Z': 'GLX',
    'J': 'XLE',
    'X': 'XAA',
    '*': 'TER'
}
AAA_DICT = {**{value: key for key, value in A_DICT.items()}, 'UNK': 'X'}

METHYLATED_LYSINES = ['MLZ', 'MLY', 'M3L']
LYSINE_ATOMS = ['N', 'CA', 'CB', 'CG', 'CD', 'CE', 'NZ', 'C', 'O']

COMMON_HETATMS = [
    'NH2',
    'HOH',
    'SO4',
]

CHAIN_IDS = list(string.ascii_uppercase + string.digits + string.ascii_lowercase)
CHAIN_IDS += [(a + b) for a in CHAIN_IDS for b in CHAIN_IDS if a != b]

# Standard accessibilities for a ALA-X-ALA tripeptide (obtained from NACCESS)
with open(op.join(DATA_DIR, 'standard_sasa.txt'), 'rt') as fin:
    _standard_sasa = fin.read().strip()
STANDARD_SASA_ALL = [[l.strip() for l in line.split()] for line in _standard_sasa.split('\n')[1:]]
STANDARD_SASA = {x[3]: float(x[4]) for x in STANDARD_SASA_ALL}

# Map modified residue ids to canonical residue ids
with open(op.join(DATA_DIR, 'residue_mapping_to_canonical.json'), 'rt') as fin:
    RESIDUE_MAPPING_TO_CANONICAL = json.load(fin)
