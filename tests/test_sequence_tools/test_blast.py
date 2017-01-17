import concurrent.futures
import functools
import logging
import os
import os.path as op
import tarfile
import time
import urllib.request

import pandas as pd

from kmtools import sequence_tools

logger = logging.getLogger(__name__)

BLAST_DB = op.abspath(op.join(op.splitext(__file__)[0], 'pdbaa'))
BLAST_DB_ARCHIVE = BLAST_DB + '.tar.gz'
BLAST_DB_FILES = [
    'pdbaa.00.phr', 'pdbaa.00.pin', 'pdbaa.00.pnd', 'pdbaa.00.pni', 'pdbaa.00.pog', 'pdbaa.00.ppd',
    'pdbaa.00.ppi', 'pdbaa.00.psd', 'pdbaa.00.psi', 'pdbaa.00.psq', 'pdbaa.pal',
]


def setup_module():
    """Download and extract PDB blast database for testing."""
    # Download blast PDB database for testing
    if not op.isfile(BLAST_DB_ARCHIVE):
        logger.debug('Downloading blast database...')
        urllib.request.urlretrieve(
            url='ftp://ftp.ncbi.nlm.nih.gov/blast/db/pdbaa.tar.gz',
            filename=BLAST_DB_ARCHIVE)
    with tarfile.open(BLAST_DB_ARCHIVE, 'r:gz') as ifh:
        logger.debug('Extracting blast database...')
        ifh.extractall(path=op.abspath(op.splitext(__file__)[0]))
    for filename in BLAST_DB_FILES:
        assert op.isfile(op.join(op.splitext(__file__)[0], filename))


def teardown_module(module):
    """Remove extracted files of the PDB blast database."""
    for filename in BLAST_DB_FILES:
        os.remove(op.join(op.splitext(__file__)[0], filename))
    for filename in BLAST_DB_FILES:
        assert not op.isfile(op.join(op.splitext(__file__)[0], filename))


def test_1():
    NUM_BLASTS = 100
    blastp = functools.partial(sequence_tools.blastp, db=BLAST_DB)
    #
    df = pd.read_csv(
        op.join(op.splitext(__file__)[0], 'elaspic_training_set_core.tsv'),
        sep='\t')
    logger.debug(df.head())
    #
    with concurrent.futures.ThreadPoolExecutor() as p:
        results_1 = p.map(blastp, df['uniprot_sequence'][:NUM_BLASTS].values)
        results_2 = p.map(blastp, df['domain_sequence'][:NUM_BLASTS].values)
    return results_1, results_2


def _blast_sequences(df, column):
    N_BLASTS = 100
    counter = 0
    time_0 = time.time()
    for i, row in df.iterrows():
        sequence = row[column]
        result_df = sequence_tools.blastp(sequence, BLAST_DB)
        assert not result_df.empty
        counter += 1
        if counter > 100:
            break
    logger.info(
        "Performed {} blast '{}' searches in {:.1f} seconds..."
        .format(N_BLASTS, column, time.time() - time_0))
