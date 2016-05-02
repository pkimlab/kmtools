# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 22:54:20 2014

@author: alexey
"""

import subprocess
import shlex
import datetime
import six

from fastcache import clru_cache

import pandas as pd
pd.options.mode.chained_assignment = None


#%%
blastp_outfmt6_column_names = [
    'query_id', 'subject_id', 'pc_identity', 'alignment_length', 'mismatches', 'gap_opens',
    'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bitscore', 'qseq', 'sseq']
blast_outfmt = "'6 qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq'"


@clru_cache(maxsize=1024, typed=False)
def call_blast(domain_sequence, db_path):
    """ blast a given sequence against a specified library
    """
    system_command = 'blastp -db {db} -evalue 0.001 -outfmt {outfmt} -max_target_seqs 100000'.format(outfmt=blast_outfmt, db=db_path)
    args = shlex.split(system_command)
    cp = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result, error_message = cp.communicate(six.b(domain_sequence))
    if error_message:
        print(error_message)
        return pd.DataFrame(columns=blastp_outfmt6_column_names), system_command
    result_buf = six.StringIO(result.decode())
    result_df = pd.read_csv(result_buf, sep='\t', names=blastp_outfmt6_column_names)
    return result_df, system_command



def annotate_blast_results(result_df, domain_start, domain_sequence_length):
    """ used by `run_blast` to add additional information to the blast results df
    """
    if result_df is None or len(result_df) == 0:
        return

    if result_df['subject_id'].values[0].count('|') == 2:
        HEADER_TYPE = 1
    elif result_df['subject_id'].values[0].count('|') == 5:
        HEADER_TYPE = 2
    else:
        raise Exception

    # Parse blast results
    result_df['pdb_id'] = result_df['subject_id'].apply(lambda x: x.split('|')[0].split('_')[0])
    result_df['pdb_chain'] = result_df['subject_id'].apply(lambda x: x.split('|')[0].split('_')[1])
    result_df['pdb_pdbfam_name'] = result_df['subject_id'].apply(lambda x: x.split('|')[1])
    result_df['pdb_pdbfam_idx'] = result_df['subject_id'].apply(lambda x: int(x.split('|')[2]))
    if HEADER_TYPE == 2:
        result_df['pdb_pfam_clan'] = result_df['subject_id'].apply(lambda x: x.split('|')[3])
        result_df['pdb_domain_def'] = result_df['subject_id'].apply(lambda x: x.split('|')[4])
        result_df['pdb_cath_id'] = result_df['subject_id'].apply(lambda x: x.split('|')[5])

    # Score sequence identity
    alpha = 0.95
    result_df['alignment_identity'] = result_df['pc_identity'] / 100.0
    result_df['alignment_coverage'] = (result_df['q_end'] - result_df['q_start'] + 1) / float(domain_sequence_length)
    result_df['alignment_score'] = (
        alpha * result_df['alignment_identity'].values * result_df['alignment_coverage'].values +
        (1 - alpha) * result_df['alignment_coverage'].values
    )

    # Domain definitions
    domain_start_new = domain_start + result_df['q_start'].astype(int) - 1
    domain_end_new = domain_start + result_df['q_end'].astype(int) - 1
    result_df['domain_start_new'] = domain_start_new
    result_df['domain_end_new'] = domain_end_new
    result_df['domain_def_new'] = domain_start_new.astype(str) + ':' + domain_end_new.astype(str)
    result_df['t_date_modified'] = datetime.datetime.now()



@clru_cache(maxsize=512, typed=False)
def run_blast(uniprot_sequence, pfam_clan, pdbfam_name, domain_def, blastp_libraries_path):
    """
    .. note:: Deprecated
        It is recommended to use the same database for all templates. Thus, querying the databaes
        by pdbfam name of pfam clan name is no longer necessary.
    """
    domain_start, domain_end = [int(x) for x in domain_def.split(':')]
    domain_sequence = uniprot_sequence[domain_start-1:domain_end]

    # Run blast
    result_df, system_command = call_blast(domain_sequence, blastp_libraries_path + pdbfam_name + '/' + pdbfam_name)
    annotate_blast_results(result_df, domain_start, len(domain_sequence))

    if pfam_clan != pdbfam_name and (result_df is None or len(result_df) == 0):
        print(system_command)
        print("No blast resutls for pdbfam name {}, domain def {}!".format(pdbfam_name, domain_def))
        print("Retryign using the blast library of the pfam clan {}...".format(pfam_clan))
        result_df, system_command = call_blast(domain_sequence, blastp_libraries_path + pfam_clan + '/' + pfam_clan)
        annotate_blast_results(result_df, domain_start, len(domain_sequence))
    elif (pfam_clan != pdbfam_name and len(result_df[
            (result_df['alignment_identity'] > 0.35) &
            (result_df['alignment_coverage'] > 0.35)]) == 0):
        print(system_command)
        print("Bad blast results for pdbfam name {}, domain def {}!".format(pdbfam_name, domain_def))
        print("Retryign using the blast library of the pfam clan {}...".format(pfam_clan))
        result_df, system_command = call_blast(domain_sequence, blastp_libraries_path + pfam_clan + '/' + pfam_clan)
        annotate_blast_results(result_df, domain_start, len(domain_sequence))

    if result_df is None or len(result_df) == 0:
        print(system_command)
        print("No blast resutls for pfam clan {}, pdbfam name {}, domain def {}!".format(pfam_clan, pdbfam_name, domain_def))

    return result_df


@clru_cache(maxsize=512, typed=False)
def run_blast_basic(sequence, blast_db_path, domain_def=None):
    """
    .. note::
        Previously, input parameters were:
        (domain_sequence, pfam_clan, pdbfam_name, domain_def, blast_db_path)
    """
    if domain_def is not None:
        domain_start, domain_end = [int(x) for x in domain_def.split(':')]
        domain_sequence = sequence[domain_start-1 : domain_end]
    else:
        domain_start = 1
        domain_sequence = sequence
    result_df, system_command = call_blast(domain_sequence, blast_db_path)
    annotate_blast_results(result_df, domain_start, len(domain_sequence))
    return result_df, system_command


#%%
call_blast.cache_clear()
run_blast.cache_clear()
run_blast_basic.cache_clear()
