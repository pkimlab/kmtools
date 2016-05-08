A_DICT = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'E': 'GLU',
    'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS',
    'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP',
    'Y': 'TYR', 'V': 'VAL', 'U': 'SEC', 'O': 'PYL',
    'B': 'ASX', 'Z': 'GLX', 'J': 'XLE', 'X': 'XAA', '*': 'TER'
}
AAA_DICT = dict([(value, key) for key, value in list(A_DICT.items())])
AAA_DICT['UNK'] = 'X'
AAA_DICT['MSE'] = 'M'
AAA_DICT['CSD'] = 'C'

# Phosphorylated residues
AAA_DICT['SEP'] = 'S'  # PHOSPHOSERINE
AAA_DICT['TPO'] = 'T'  # PHOSPHOTHREONINE
AAA_DICT['SEP'] = 'Y'  # O-PHOSPHOTYROSINE

# Methylated lysines
AAA_DICT['MLZ'] = 'K'
AAA_DICT['MLY'] = 'K'
AAA_DICT['M3L'] = 'K'
