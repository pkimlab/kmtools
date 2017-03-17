"""
Functions to deal with DNA, translate, recode to avoid restriction enzymes, ...



"""
from __future__ import print_function

import bisect
import re
import random
import logging
import operator
import functools

from codons_info import A2C_DICT, STOP_CODONS, USAGE_FREQ, A2C_NNS_DICT

logger = logging.getLogger(__name__)


def translate2aa(nseq, start=1):
    """Return AA sequence, from NA sequence (string).

    Parameters
    ----------

    nseq : str
        Nucleotide sequence
    start : int
        Start to translate from the position, by default 1

    Returns
    -------
     str with the aminoacid sequence.


    """
    bases = ['T', 'C', 'A', 'G']
    codons = [a + b + c for a in bases for b in bases for c in bases]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    codon_table = dict(zip(codons, amino_acids))
    pseq = ""

    # the start point, is where the aligment begins, not the protein, that give us the read frame
    # get the seq from this point until I found a stop codon

    i = start - 1

    while True:
        codon = nseq[i:i + 3]
        i += 3

        if len(codon) != 3:
            break

        aminoacid = codon_table[codon]
        if aminoacid == "*":
            break

        pseq += aminoacid

    return pseq


def translate2na(seq, specie='human'):
    """Return a Nucleotide seq for a aminoacid sequences. the codon will be chosed
        according the codon frequency of each specie, by default human

    Paramenters
    -----------
    seq : str
        Amino acid, in One letter code sequences

    Returns
    -------
    str
    na seq (str) codon

    """
    seq_na = []
    for a in seq:
        codons = A2C_DICT.get(a)
        seq_na.append(codon_weighted_random_choice(codons, specie))

    # print(random.choice(foo))
    return ''.join(seq_na)


def codon_weighted_random_choice(codons, specie):
    """.

    The probability for each element elem in list to be selected is weighted
    by weight(elem) weight_dictionary`` must be a callable accepting one argument,
    and returning a non-negative number. If ``weight(elem)`` is zero, ``elem``
    will not be  considered.

    Parameters
    ----------

    codons : array_like
        must be an iterable containing more than one element.

    specie : str
        Codon usage specie, human, e.coli, etc.

    Return:
    -------
    str
        Returns a  element from a list. codon base on codon usage probability

    Raises
    ------
    ValueError
        If the codon usage library is not valid

    """
    try:
        weight_dictionary = USAGE_FREQ.get(specie)
    except ValueError:
        # may be this could be a warning, and call human codon usage
        raise ValueError('{} is not a valid specie'.format(specie))

    weights = 0
    elems = []
    for elem in codons:
        w = weight_dictionary.get(elem)
        try:
            is_neg = w < 0
        except TypeError:
            raise ValueError("Weight of element '%s' is not a number (%s)" %
                             (elem, w))
        if is_neg:
            raise ValueError("Weight of element '%s' is negative (%s)" %
                             (elem, w))
        if w != 0:
            try:
                weights += w
            except TypeError:
                raise ValueError("Weight of element '%s' is not a number "
                                 "(%s)" % (elem, w))
            elems.append((weights, elem))
    if not elems:
        raise ValueError("Empty sequence")
        print('{} {}'.format(codons, weight_dictionary))
    ix = bisect.bisect(elems, (random.uniform(0, weights), None))
    # print ix
    return elems[ix][1]


def clean_restriction_sites(naseq, dict_restriction=['GAATTC', 'CCCGGG', 'GTCGAC']):
    """.

    Check if there is a restriction site for any of the enzymes in
        self.set_restriction_enzyme for a sequence.
        if it could find one, the afected codon is retranslated,
        that will generate a codon selection,  and it check again.
        Formed called _check_4_restricitions

    Parameters
    ----------

        :param naseq (str): ADN sequences to check
        :param dict_restriction (list or dict): Restriction sites defintions.
            (default) dict_restriction = ['GAATTC','CCCGGG','GTCGAC'] #or ecorI, XmaI, salI

    return
    ------
        naseq (str): if the seq contains restriction site, it will be returned recoded
                  to avoid them.


    """
    clean = False

    ilimit = 0
    while not clean:
        ilimit += 1

        matches = has_restriction_sites(naseq, dict_restriction)
        if len(matches) == 0:
            clean = True
        else:
            naseq = reshufle_seq(naseq, matches)

        if ilimit == 1000000:
            print('WARNING the cycles limit has been pass and the seq stil\
                    contain a restriction site  {}'.format(naseq))
            return naseq

    return naseq


def reshufle_seq(seq, position_pairs):
    """Resampling codon usage. Avoid restriction enzyme sites.


    Paramaters:
    -----------
        :param seq (str): Nucleotide sequence to be reshufled
        :param position_pairs (list): list of list of positions.

    """
    for restriction_match in position_pairs:
        i = 0
        afected = range(restriction_match[0], restriction_match[1] + 1)
        while i < max(afected):
            if i in afected:

                # This should return 0,1, or 2
                codon_coordinates = i % 3
                # Use the module to find the codon in this postion
                codon = seq[i - codon_coordinates:i + 3 - codon_coordinates]
                aminoacid = translate2aa(codon)
                # amino  acids with only one codon.
                if aminoacid in ['M', 'W']:
                    i += 1
                    continue
                else:

                    alternatives = list(A2C_DICT.get(aminoacid))
                    alternatives.remove(codon)
                    new_codon = random.choice(alternatives)
                    seq = seq[:i - codon_coordinates] + \
                        new_codon + seq[i + 3 - codon_coordinates:]
                    # Go to the next codon
                    i += i + 3 - codon_coordinates
            else:
                i += 1

    return seq


def has_restriction_sites(seq, dict_restriction):
    """.

    Former match restrictions. Check if there are any restriction site in the sequences,
    and if it have, return a list with the positions involved

    Paramaters
    ----------
        :param seq (str); DNA sequence
        :param dict_restriction (list or dict) with restriction enzymes

    return
    ------
        postions involved (list) [[strat,end],[strat,end]]
        the list is emptty if there is no restriction sites

    """
    matches = []

    # if the the restrictions are a dict, extract sites
    if isinstance(dict_restriction, dict):
        dict_restriction = dict_restriction.values()

    for restriction in dict_restriction:
        # print 'search',restriction,seq
        hit = re.search(restriction, seq)
        if hit:

            matches.append([hit.start(), hit.end()])

    return matches


def code_4_any(seq, filter_codons=STOP_CODONS):
    """Return True if the sequences contain any of the codons in filter_codons.

    Parameters
    ----------
    seq : str
        sequence to check, in the right frame shift

    filter_codons : array_like
        iterable with the codons to test (by default stop codons)

    Returns
    -------

    """
    reading_frame = len(seq) / 3
    jdx = 0

    for idx in range(int(reading_frame)):
        if seq[jdx:jdx + 3] in filter_codons:
            return True
        else:
            pass
        jdx += 3

    return False


def is_bias_on(seq, filter_by=A2C_NNS_DICT):
    """Return True if the sequences contain only codons in filter_codons.

    Parameters
    ----------
    seq : str
        sequence to check, in the right frame shift

    filter_by : dict
        Amino to Codon dictionary with the codons to filter (default NNS)

    Returns
    -------

    """
    filter_codons = _get_codons_from(filter_by)

    reading_frame = len(seq) / 3
    jdx = 0

    for idx in range(int(reading_frame)):
        if seq[jdx:jdx + 3] in filter_codons:
            pass
        else:
            return False
        jdx += 3

    return True


def possible_encondings(seq, codons_base=A2C_NNS_DICT):
    """Return the number of possible combinations of nucleotides to encode a AA seq.

    Parameters
    ----------
    seq : str
       Amino acid sequence

    codons_base : dict
        Dictionary with the codons list for each aminoacid, by default NNS

    Returns
    -------
    int
       The number of different Nucleotide sequence that may encode the input seq

    """

    n = list()

    for a in seq:
        n.append(len(codons_base[a]))

    return functools.reduce(operator.mul, n)


def _get_codons_from(codons_dict):
    """Extract all codons from dictionary in a single list.

    Parameters
    ----------
    codons_dict : dict

    Returns
    -------
    list


    """
    all_codons = list()
    for codons_list in codons_dict.values():
        all_codons.extend(codons_list)

    return all_codons
