"""
useful Sequence Analysis tools



"""
import logging
import subprocess

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)



def generate_logo( sequences, seq_len = 80, filename='designs'):
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

    ohandler =  open(filename+'.fasta', 'w')
    for seq in sequences:
        print(">{}".format(seq),file=ohandler)
        print("{}".format(seq),file=ohandler)

    ohandler.close()
    # quick and dirty logo generartion
    command = subprocess.Popen('weblogo -f {} -c chemistry -o {} -F pdf  -n {} -U bits --composition equiprobable '.format(filename+'.fasta',
                                                                            filename+'.pdf',
                                                                            seq_len),
                    shell=True)
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
        for idx,p in enumerate(seq):
            logger.debug("%i %s",idx,p)
            if idx not in matrix:
                matrix[idx] = {p:1}
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
            row.append(frequencies[p].get(a, 0.0)/N)
        matrix.append(row)

    # convert to pandas.Dataframe
    m = pd.DataFrame(matrix, index=amino,columns=list(range(seq_len)))
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
            prob = frequencies[p].get(a, 0.0)/N
            row.append(np.log2(prob/.05))
        matrix.append(row)

    # convert to pandas.Dataframe
    m = pd.DataFrame(matrix, index=amino,columns=list(range(seq_len)))
    logger.debug(m)
    return m


def binding_score(seq, pwm):
    """Score a sequence using a PWM.

    Parameters
    ----------
    seq: str

    pwm: Pandas.DataFrame

    Returns
    -------

    float

    """

    score = list()
    for pos, symbol in enumerate(seq):
        score.append(pwm.at[symbol, pos])

    return sum(score)

def dist_PWM(pwm1, pwm2):
    """Euclidina distance in the PWM.

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
            rows.append((pwm1.at[r,c]-pwm2.at[r,c])**2)
        colum.append(sum(rows)*.5)

    return sum(colum)/float(w)
