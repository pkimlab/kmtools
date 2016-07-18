from __future__ import print_function
import subprocess
import random
import os


def align_pairwise(sequence_ref, sequence_alt):
    """.

    Examples
    --------
    >>> align_pairwise('AAAGGVVV', 'AAAGGVVV')
    ('AAAGGVVV', 'AAAGGVVV')
    >>> align_pairwise('AAAGGVVVAAA', 'AAAGGAAA')
    ('AAAGGVVVAAA', 'AAAGG---AAA')
    """
    seqs = ">1\n{}\n>2\n{}\n".format(sequence_ref, sequence_alt)
    p = subprocess.run(
        ['muscle'],
        input=seqs,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        check=True)
    __, alignment_ref, alignment_alt = p.stdout.split('>')
    alignment_ref = alignment_ref.strip('\n 12').replace('\n', '')
    alignment_alt = alignment_alt.strip('\n 12').replace('\n', '')
    assert len(alignment_ref) == len(alignment_alt)
    return alignment_ref, alignment_alt


def get_crossmapping(alignment_ref, alignment_alt):
    """.

    Examples
    --------
    >>> get_crossmapping('ABCDF', 'A-CEF')
    ('1,3,,5', '1,,2,,4')
    """
    assert len(alignment_ref) == len(alignment_alt)
    mapping_ref2alt, mapping_alt2ref = [], []
    a_gaps, b_gaps = 0, 0
    for i, (a, b) in enumerate(zip(alignment_ref, alignment_alt)):
        if a == b:
            mapping_ref2alt.append(str(i - a_gaps + 1))
            mapping_alt2ref.append(str(i - b_gaps + 1))
        elif a == '-':
            a_gaps += 1
            mapping_ref2alt.append('')
            # mapping_alt2ref.append('')
        elif b == '-':
            b_gaps += 1
            # mapping_ref2alt.append('')
            mapping_alt2ref.append('')
        else:
            mapping_ref2alt.append('')
            mapping_alt2ref.append('')
    return ','.join(mapping_ref2alt), ','.join(mapping_alt2ref)

##
# Alignment with needle


def align_sw(seqa, seqb, **kargs):
    """.

    Smith-Waterman aligment
    central function, call  to get score, or identiy, or aligment etc
    """
    outfile = 'aln.needle2.{}.out'.format(random.randint(1, 30000))

    call_needle(seqa, seqb, outfile=outfile, **kargs)

    return load_scores(outfile, **kargs)


def call_needle(seqa, seqb, gapopen=12.0, gapextend=1.0, outfile='aln.needle2.out', **kargs):

    # overwrite_files
    filename_a = 'file_{}.seq'.format(random.randint(1, 30000))
    filename_b = 'file_{}.seq'.format(random.randint(30001, 90000))
    filea = open(filename_a, 'w')
    fileb = open(filename_b, 'w')
    print(seqa.strip(), file=filea)
    print(seqb.strip(), file=fileb)
    filea.close()
    fileb.close()
    # exec_needle
    status = subprocess.call([
        'needle', '-asequence', filename_a, '-bsequence', filename_b,
        '-gapopen', str(gapopen), '-gapextend', str(gapextend), '-outfile', outfile
    ])
    os.remove(filename_a)
    os.remove(filename_b)

    return status


def load_scores(outfile, lenght=16.0):
    ident = 0.0
    sim = 0

    with open(outfile, 'r') as input_file:

        for line in input_file:
            # print line
            tabdata = line.split()
            if len(tabdata) > 3:
                if tabdata[1] == 'Identity:':
                    matches, target = tabdata[2].split('/')
                    ident = float(matches) / lenght
                    # match = tabdata[3][1:-2]
                if tabdata[1] == 'Similarity:':
                    matches, target = tabdata[2].split('/')
                    sim = float(matches) / lenght
                    # sim = tabdata[3][1:-2]
                    input_file.close()
                    return ident * 100, sim * 100

    return
