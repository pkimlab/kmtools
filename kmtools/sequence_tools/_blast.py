import io
import logging
import shlex
import subprocess

import pandas as pd

import kmtools.sequence_tools

logger = logging.getLogger(__name__)

BLAST_OUTFMT6 = """\
'6 qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq'\
"""
BLAST_OUTFMT6_COLUMN_NAMES = [
    'query_id', 'subject_id', 'pc_identity', 'alignment_length', 'mismatches', 'gap_opens',
    'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bitscore', 'qseq', 'sseq',
]


# @lru_cache(maxsize=1024, typed=False)
def blastp(sequence, db, evalue=0.001, max_target_seqs=100000):
    """Run `blastp`.

    .. note::

        It is often useful to wrap this function inside `functools.lru_cache`
        in order to speed up repeated queries.

    Parameters
    ----------
    sequence : str
        Query sequence.
    blast_db : str
        Full path to the blast database.
    """
    system_command = (
        'blastp -db {db} -outfmt {outfmt} -evalue {evalue} -max_target_seqs {max_target_seqs}'
        .format(db=db, outfmt=BLAST_OUTFMT6, evalue=evalue, max_target_seqs=max_target_seqs)
    )
    logger.debug(system_command)
    cp = subprocess.Popen(
        shlex.split(system_command),
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        universal_newlines=True)
    result, error_message = cp.communicate(sequence)
    if error_message:
        logger.error(error_message)
    result_df = pd.read_csv(io.StringIO(result), sep='\t', names=BLAST_OUTFMT6_COLUMN_NAMES)
    return result_df


def expand_blast_results(result_df: pd.DataFrame) -> pd.DataFrame:
    result_df = result_df.copy()
    result_df['a2b'], result_df['b2a'] = list(zip(*(
        kmtools.sequence_tools.get_crossmapping(*x, skip_mismatch=False)
        for x in result_df[['qseq', 'sseq']].values)))
    assert (result_df['qseq'].str.len() == result_df['sseq'].str.len()).all()
    return result_df
