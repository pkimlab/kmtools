from __future__ import print_function
from .codons_info import *
import bisect
import re
import random


def translate2aa(nseq, start=1):
    """.

    Return AA sequence, from NA sequence (string)

    Args:
    -----
        :param str nseq: Nucleotide sequence
        :param int start: Start to translate from the
         position, by default 1

        :return: str with the aminoacid sequence.
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
    """.

    Return a Nucleotide seq for a aminoacid sequences. the codon will be chosed
        according the codon frequency of each specie, by default human

        Args:
        =====
            seq (str): Amino acid, in One letter code sequences

        Return:
           na seq (str): codon
    """
    seq_na = []
    for a in seq:
        codons = codons_info.A2C_DICT.get(a)
        seq_na.append(codon_weighted_random_choice(codons, specie))

    # print(random.choice(foo))
    return ''.join(seq_na)


def codon_weighted_random_choice(codons, specie):
    """.

    Returns a random element from a list. The probability for
    each element elem in list to be selected is weighted
    by weight(elem).
    weight_dictionary`` must be a callable accepting one argument,
    and returning a  non-negative number. If ``weight(elem)`` is zero,
    ``elem`` will not be considered.

    Args:
    -----
        :param codons (list): must be an iterable containing more than one element.
        :param specie (str): Codon usage specie, human, e.coli, etc.



    Return:
    -------

        codon (str)

    """
    weight_dictionary = codons_info.USAGE_FREQ.get(specie)
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
    """.

    Resampling codon usage. Avoid restriction enzyme sites

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

                    alternatives = list(codons_info.A2C_DICT.get(aminoacid))
                    alternatives.remove(codon)
                    new_codon = random.choice(alternatives)
                    seq = seq[:i - codon_coordinates] + new_codon + seq[i + 3 - codon_coordinates:]
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
