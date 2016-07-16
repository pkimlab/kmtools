import subprocess


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
