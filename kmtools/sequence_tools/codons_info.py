
# Codon Usage probability for each scpecie'
USAGE_FREQ = {'E.coli':{'GGG': 0.15,'GGA': 0.11,'GGT': 0.34,'GGC': 0.4,\
                              'GAG': 0.31,'GAA': 0.69,'GAT': 0.63,'GAC': 0.37,\
                              'GTG': 0.37,'GTA': 0.15,'GTT': 0.26,'GTC': 0.22,\
                              'GCG': 0.36,'GCA': 0.21,'GCT': 0.16,'GCC': 0.27,\
                              'AGG': 0.02,'AGA': 0.04,'CGG': 0.1,'CGA': 0.06,\
                              'CGT': 0.38,'CGC': 0.4,'AAG': 0.23,'AAA': 0.77,\
                              'AAT': 0.45,'AAC': 0.55,'ATG': 1.0,'ATA': 0.07,\
                              'ATT': 0.51,'ATC': 0.42,'ACG': 0.27,'ACA': 0.13,\
                              'ACT': 0.17,'ACC': 0.44,'TGG': 1.0,'TGT': 0.45,\
                              'TGC': 0.55,'TAG': 0.07,'TAA': 0.64,'TGA': 0.29,\
                              'TAT': 0.57,'TAC': 0.43,'TTT': 0.57,'TTC': 0.43,\
                              'AGT': 0.15,'AGC': 0.28,'TCG': 0.15,'TCA': 0.12,\
                              'TCT': 0.15,'TCC': 0.15,'CAG': 0.65,'CAA': 0.35,\
                              'CAT': 0.57,'CAC': 0.43,'TTG': 0.13,'TTA': 0.13,\
                              'CTG': 0.5,'CTA': 0.04,'CTT': 0.1,'CTC': 0.1,\
                              'CCG': 0.52,'CCA': 0.19,'CCT': 0.16,'CCC': 0.12},\
                    'human':{'CTT': 0.13, 'ACC': 0.36, 'ACA': 0.28,\
                        'AAA': 0.42, 'ATC': 0.48, 'AAC': 0.54, 'ATA': 0.16,\
                        'AGG': 0.2, 'CCT': 0.28, 'ACT': 0.24, 'AGC': 0.24,\
                        'AAG': 0.58, 'AGA': 0.2, 'CAT': 0.41, 'AAT': 0.46,\
                        'ATT': 0.36, 'CTG': 0.41, 'CTA': 0.07, 'CTC': 0.2,\
                        'CAC': 0.59, 'ACG': 0.12, 'CAA': 0.25, 'AGT': 0.15,\
                        'CCA': 0.27, 'CCG': 0.11, 'CCC': 0.33, 'TAT': 0.43,\
                        'GGT': 0.16, 'TGT': 0.45, 'CGA': 0.11, 'CAG': 0.75,\
                        'TCT': 0.18, 'GAT': 0.46, 'CGG': 0.21, 'TTT': 0.45,\
                        'TGC': 0.55, 'GGG': 0.25, 'TAG': 0.2, 'GGA': 0.25,\
                        'TGG': 1.0, 'GGC': 0.34, 'TAC': 0.57, 'TTC': 0.55,\
                        'TCG': 0.06, 'TTA': 0.07, 'TTG': 0.13, 'CGT': 0.08,\
                        'GAA': 0.42, 'TAA': 0.28, 'GCA': 0.23, 'GTA': 0.11,\
                        'GCC': 0.4, 'GTC': 0.24, 'GCG': 0.11, 'GTG': 0.47,\
                        'GAG': 0.58, 'GTT': 0.18, 'GCT': 0.26, 'TGA': 0.52,\
                        'GAC': 0.54, 'TCC': 0.22, 'TCA': 0.15, 'ATG': 1.0,\
                        'CGC': 0.19}
        }


# Aminoacid to codon translation table
A2C_DICT =  {'I' : [ u'ATT',u'ATC',u'ATA' ],
                            'L' : [ u'CTT', u'CTC', u'CTA', u'CTG', u'TTA', u'TTG' ],
                            'V' : [ u'GTT', u'GTC', u'GTA', u'GTG' ],
                            'F' : [ u'TTT', u'TTC' ],
                            'M' : [ u'ATG' ],
                            'C' : [ u'TGT', u'TGC' ],
                            'A' : [ u'GCT',u'GCC', u'GCA',u'GCG' ],
                            'G' : [ u'GGT', u'GGC',u'GGA', u'GGG' ],
                            'P' : [ u'CCT', u'CCC', u'CCA', u'CCG' ],
                            'T' : [ u'ACT',u'ACC', u'ACA', u'ACG' ],
                            'S' : [ u'TCT', u'TCC', u'TCA', u'TCG', u'AGT', u'AGC' ],
                            'Y' : [ u'TAT', u'TAC' ],
                            'W' : [ u'TGG' ],
                            'Q' : [ u'CAA', u'CAG' ],
                            'N' : [ u'AAT', u'AAC' ],
                            'H' : [ u'CAT' ,u'CAC' ],
                            'E' : [ u'GAA', u'GAG' ],
                            'D' : [ u'GAT', u'GAC' ],
                            'K' : [ u'AAA', u'AAG' ],
                            'R' : [ u'CGT', u'CGC' ,u'CGA', u'CGG', u'AGA', u'AGG' ],
                            '*' : [ u'TAA', u'TAG' ,u'TGA' ]}

# Aminoacid to codon translation table
A2C_NNS_DICT =  {'I' : [u'ATC' ],
                            'L' : [ u'CTC', u'CTG', u'TTG' ],
                            'V' : [ u'GTC', u'GTG' ],
                            'F' : [ u'TTC' ],
                            'M' : [ u'ATG' ],
                            'C' : [ u'TGC' ],
                            'A' : [ u'GCC', u'GCG' ],
                            'G' : [ u'GGC', u'GGG' ],
                            'P' : [ u'CCC', u'CCG' ],
                            'T' : [ u'ACC',  u'ACG' ],
                            'S' : [ u'TCC', u'TCG', u'AGC' ],
                            'Y' : [ u'TAC' ],
                            'W' : [ u'TGG' ],
                            'Q' : [ u'CAG' ],
                            'N' : [ u'AAC' ],
                            'H' : [ u'CAC' ],
                            'E' : [ u'GAG' ],
                            'D' : [ u'GAC' ],
                            'K' : [ u'AAG' ],
                            'R' : [ u'CGC' , u'CGG', u'AGG' ],
                            '*' : [ u'TAG' ]}


# codon to Aminoacid translation table
C2A_DICT =  {u'ATT':'I', u'ATC':'I', u'ATA':'I',
             u'CTT':'L', u'CTC':'L', u'CTA':'L', u'CTG':'L', u'TTA':'L', u'TTG':'L',
             u'GTT':'V', u'GTC':'V', u'GTA':'V', u'GTG' :'V',
             u'TTT':'F', u'TTC':'F',
             u'ATG':'M',
             u'TGT':'C', u'TGC':'C',
             u'GCT':'A', u'GCC':'A', u'GCA':'A', u'GCG':'A',
             u'GGT':'G', u'GGC':'G', u'GGA':'G', u'GGG':'G',
             u'CCT':'P', u'CCC':'P', u'CCA':'P', u'CCG':'P',
             u'ACT':'T', u'ACC':'T', u'ACA':'T', u'ACG':'T',
             u'TCT':'S', u'TCC':'S', u'TCA':'S', u'TCG':'S', u'AGT':'S', u'AGC':'S',
             u'TAT':'Y', u'TAC':'Y',
             u'TGG':'W',
             u'CAA':'Q', u'CAG':'Q',
             u'AAT':'N', u'AAC':'N',
             u'CAT':'H', u'CAC':'H',
             u'GAA':'E', u'GAG':'E',
             u'GAT':'D', u'GAC':'D',
             u'AAA':'K', u'AAG':'K',
             u'CGT':'R', u'CGC':'R', u'CGA':'R', u'CGG':'R', u'AGA':'R', u'AGG':'R',
             u'TAA':'*', u'TAG':'*', u'TGA':'*'}