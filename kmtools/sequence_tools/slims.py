"""
Tools for Short Linear Motifs on Protein sequences - SLiM -

"""
import os
import re
import logging

import pandas as pd
import numpy as np
import subprocess

logger = logging.getLogger(__name__)
# To improve. Motif definitios repositories
absolute_path = os.path.dirname(__file__)


# ELM_Classes_Download_Version: 1.4
# ELM_Classes_Download_Date: 2016-05-09 10:55:27.098737
# Origin: elm.eu.org
# Type: tsv
# Num_Classes: 247
ELM = pd.read_csv(absolute_path + "/support/elms_classes.csv")


# PepX database motives
#
PEPX = pd.read_csv(absolute_path + "/support/elms_pepx.csv")


def search_motif(seq, ELMdefinitions=ELM):
    """Return all the Motif in a sequences

    Parameters
    ----------
    seq: str
        sequence to query

    ELMdefinitions: str, list, dict, pandas.Dataframe['Regex']
        Motif definitions, by default use  /support/elm_classes.tsv loaded to a pandas.DataFrame.
        Alternative you can use PepXELM_SLiM_defs - SLiM of domains in PepX.
        Custom regular expression can be inserted as str r'P.L', list [r'P.L',r'P.L.Y',r'L...P.L']
        or  dict {'PxL SHORT':r'P.L','Experimental PxLxY':r'P.L.Y'}


    Return
    ------
        Returns a Datafram with all the motif in the sequnece
        columns {Accession,Identifier,Description,Regex,Probability,start,end}

    """
    # Check input Parameter type and transform to pandas DataFrame
    if isinstance(ELMdefinitions, str):
        ELMdefinitions = list(ELMdefinitions)

    if isinstance(ELMdefinitions, list):
        ELMdefinitions = _convert_list(ELMdefinitions)

    if isinstance(ELMdefinitions, dict):
        ELMdefinitions = _convert_dict(ELMdefinitions)

    elmhits = list()
    for _, elm in ELMdefinitions.iterrows():
        # regex = Def[i][1]
        m = re.search(elm['Regex'], seq)
        if str(m) != "None":
            elm['start'] = m.start()
            elm['end'] = m.end()
            elmhits.append(elm)

    return pd.DataFrame(elmhits)


def _convert_list(slim_list):
    """Convert a list to a DataFrame.

    Parameters
    ----------

    Returns
    -------

    """
    # [r'P.L',r'P.L.Y',r'L...P.L']
    paried_list = [[i, x] for i, x in enumerate(slim_list)]
    return pd.DataFrame(paried_list, columns=['Identifier', 'Regex'])


def _convert_dict(slim_dict):
    """Convert a dictionary to a pandas DataFrame.

    Parameters
    ----------

    Returns
    -------

    """
    # {'PxL SHORT':r'P.L','Experimental PxLxY':r'P.L.Y'}
    return pd.DataFrame.from_dict(slim_dict,
                                  orient='index'
                                  ).reset_index().rename(columns={0: 'Regex',
                                                                  'index': 'Identifier'})


def _check_elm_dataframe(df):
    """Check DataFrame Format.

    Parameters
    ----------

    Returns
    -------

    """

    if 'Regex' in df.columns and df.shape[0] > 0:
        return True
    else:
        return False


# Test
# motif('A1L3X0', 'MNSVGEACTDMKREYDQCFNRWFAEKFLKGDSSGDPCTDLFKRYQQCVQKAIKEKEIPIEGLEFMGHGKEKPENSS')
# a = motifsearch('MNSVGEACTDMKREYDQCFNRWFAEKFLKGDSSGDPCTDLFKRYQQCVQKAIKEKEIPIEGLEFMGHGKEKPENSS\
    # PPPPPPPP')
# print a
# print a.shape

# A N A L Y S I S


def generate_logo(sequences, seq_len=80, filename='designs'):
    """quick logo generation.

    Parameters
    ----------
    sequences : array_like

    seq_len : int


    filename : str

    Returns
    -------

    """
    # if pass , Folder name

    ohandler = open(filename + '.fasta', 'w')
    for seq in sequences:
        print(">{}".format(seq), file=ohandler)
        print("{}".format(seq), file=ohandler)

    ohandler.close()
    # quick and dirty logo generartion
    command = subprocess.Popen('weblogo -f {} -c chemistry -o {} -F pdf -n {}'
                               '-U bits --composition equiprobable '.format(filename + '.fasta',
                                                                            filename + '.pdf',
                                                                            seq_len), shell=True)
    command.wait()

    return


#
# - PSSM
#

def get_pfm(sequences):
    """Read list of sequence  and return frequency position Matrix.

    Parameters
    ----------
    sequences : array_like
       array of list containing the sequence in str format

    Returns
    -------
    dict
       Matrix on a dict format  where the key is the position, for each position the value is other
       dict with the counts  per amino acid i.e {'A': 44, 'C': 3, }


    """
    matrix = dict()
    for seq in sequences:
        for idx, p in enumerate(seq):
            logger.debug("%i %s", idx, p)
            if idx not in matrix:
                matrix[idx] = {p: 1}
            else:
                matrix[idx][p] = 1 + matrix[idx].get(p, 0)

    return matrix


def get_ppm(sequences):
    """Generate Position Probability Matrix

    Parameters
    ----------
    sequences : array_like
        array of list containing the sequence in str format

    Returns
    -------
    pandas.DataFrame
         Position Probability Matrix , index is aminoacid, each column is the position
    """

    amino = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C",
             "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W"]

    frequencies = get_pfm(sequences)
    logger.debug(frequencies)

    N = len(sequences)
    # get first element of the dict without know keys
    seq_len = len(sequences[0])
    matrix = list()
    for a in amino:
        # row = [a]
        row = list()
        for p in range(seq_len):
            row.append(frequencies[p].get(a, 0.0) / N)
        matrix.append(row)

    # convert to pandas.Dataframe
    m = pd.DataFrame(matrix, index=amino, columns=list(range(seq_len)))
    logger.debug(m)
    return m


def get_pwm(sequences):
    """Generate Position Weights Matrix

    Parameters
    ----------
    sequences : array_like
        array of list containing the sequence in str format

    Returns
    -------
    pandas.DataFrame
         Position Probability Matrix , index is aminoacid, each column is the position
    """

    amino = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C",
             "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W"]

    frequencies = get_pfm(sequences)
    logger.debug(frequencies)

    N = len(sequences)
    # get first element of the dict without know keys
    seq_len = len(sequences[0])
    matrix = list()
    for a in amino:
        # row = [a]
        row = list()
        for p in range(seq_len):
            prob = frequencies[p].get(a, 0.0) / N
            row.append(np.log2(prob / .05))
        matrix.append(row)

    # convert to pandas.Dataframe
    m = pd.DataFrame(matrix, index=amino, columns=list(range(seq_len)))
    logger.debug(m)
    return m


def binding_score(seq, pwm):
    """Score a sequence using a PWM.

    Parameters
    ----------
    seq: str

    pwm: :class:`pandas.DataFrame`

    Returns
    -------
    float
    """
    score = list()
    for pos, symbol in enumerate(seq):
        score.append(pwm.at[symbol, pos])

    return sum(score)


def dist_PWM(pwm1, pwm2):
    """Euclidian distance between two PWM.

    $ D_{(a,b)} = \frac{1}{W} \sum_{i=1}^{w} \frac{1}{2} \sum (a_{i,L} - b_{i,L})^2 $

    Parameters
    ----------

    Returns
    -------

    """
    assert pwm1.shape == pwm2.shape

    w = len(pwm1.columns)

    colum = list()

    for c in pwm1.columns:
        rows = list()
        for r in pwm1.index:
            rows.append((pwm1.at[r, c] - pwm2.at[r, c])**2)
        colum.append(sum(rows) * .5)

    return sum(colum) / float(w)
