import functools
import io
import logging
import re
import warnings

import Bio.SeqIO
import numpy as np
import pandas as pd
import requests

logger = logging.getLogger(__name__)


@functools.lru_cache(maxsize=512)
def fetch_sequence(sequence_id, database="uniprot"):
    """Fetch sequence from remote database.

    Parameters
    ----------
    sequence_id: str
        Sequence accession.

    Examples
    --------
    >>> str(fetch_sequence('P13501').seq)
    'MKVSAAALAVILIATALCAPASASPYSSDTTPCCFAYIARPLPRAHIKEYFYTSGKCSNPAVVFVTRKNRQVCANPEKKWVREYINSLEMS'
    """
    if sequence_id.startswith("UPI"):
        database = "uniparc"
        url_template = "http://www.uniprot.org/uniparc/{}.fasta"
    elif sequence_id.startswith("UniRef"):
        database = "uniref"
        url_template = "http://www.uniprot.org/uniref/{}.fasta"
    else:
        database = "uniprot"
        url_template = "http://www.uniprot.org/uniprot/{}.fasta"

    url = url_template.format(sequence_id)
    logger.debug("Downloading sequence {} from {}...".format(sequence_id, url))

    r = requests.get(url)
    if r.status_code != 200:
        raise Exception(
            "Failed to fetch sequence with return code: {}".format(r.status_code)
        )

    seq = Bio.SeqIO.read(io.StringIO(r.text), "fasta")
    if database == "uniprot":
        seq.annotations["db"], seq.id, seq.name = re.split("[\| ]", seq.id)
    return seq


def mutation_matches_sequence(mutation, sequence):
    """Make sure that mutation matches protein sequence.

    Examples
    --------
    >>> mutation_matches_sequence('A1B', 'AAAAA')
    True
    >>> mutation_matches_sequence('B1A', 'AAAAA')
    False
    >>> mutation_matches_sequence('A10B', 'AAAAA')
    False
    >>> mutation_matches_sequence('A1B', None)
    nan
    >>> mutation_matches_sequence(None, 'AAAAA')
    nan
    >>> mutation_matches_sequence('A10GG', 'AAAAA')
    nan
    """
    if pd.isnull(mutation) or pd.isnull(sequence):
        return np.nan
    try:
        mutation_pos = int(mutation[1:-1])
    except ValueError:
        logger.warning("Could not extract position from mutation '{}'".format(mutation))
        return np.nan
    return mutation_pos < len(sequence) and sequence[mutation_pos - 1] == mutation[0]


def mutation_in_domain(mutation, domain):
    """Make sure that mutation falls inside the specified domain.

    Examples
    --------
    >>> mutation_in_domain('A10B', '10:20')
    True
    >>> mutation_in_domain('A100B', '10:100')
    True
    >>> mutation_in_domain('A10B', '1:20')
    True
    >>> mutation_in_domain('A10B', '11:20')
    False
    """
    if pd.isnull(mutation) or pd.isnull(domain):
        return np.nan
    try:
        mutation_pos = int(mutation[1:-1])
    except ValueError:
        logger.error("Could not extract position from mutation '{}'!".format(mutation))
        return np.nan
    domain_range = range(int(domain.split(":")[0]), int(domain.split(":")[1]) + 1)
    return mutation_pos in domain_range


def format_hgvs_mutation(mutation_refseq_aa):
    """Convert HGVS mutation into (Protein ID, Mutation, Position) tuple.

    Examples
    --------
    >>> format_hgvs_mutation('NP_000005:p.C972Y')
    ('NP_000005', 'C972Y', 972)
    >>> format_hgvs_mutation('NP_000005.2:p.V1000I')
    ('NP_000005', 'V1000I', 1000)
    >>> format_hgvs_mutation(None)
    (nan, nan, nan)
    >>> format_hgvs_mutation(np.nan)
    (nan, nan, nan)
    """
    if pd.isnull(mutation_refseq_aa):
        return (np.nan, np.nan, np.nan)
    refseq_base_id = mutation_refseq_aa.split(":")[0].split(".")[0]
    refseq_mutation = mutation_refseq_aa.split(":")[-1].lstrip("p.")
    refseq_mutation_pos = int(refseq_mutation[1:-1])
    return refseq_base_id, refseq_mutation, refseq_mutation_pos
