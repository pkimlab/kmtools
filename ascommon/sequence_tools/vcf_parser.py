import sys
import re
import gzip

header = [
    'chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info',
    'gene', 'strand', 'cds', 'aa', 'cnt', 'aaref', 'aapos', 'aaalt', 'mutation',
]

RE_MISSENSE_MUTATION = re.compile('p\.([A-Za-z][0-9]+[A-Za-z])$')


def format_mutation_missense(mutation):
    r""".

    Examples
    --------
    >>> format_mutation_missense('p.R985Q')
    ['R', '985', 'Q', 'R985Q']
    >>> format_mutation_missense('p.R985*')
    ['\\N', '\\N', '\\N', '\\N']
    """
    match = RE_MISSENSE_MUTATION.match(mutation)
    if match is None:
        return ['\\N', '\\N', '\\N', '\\N']
    mutation = match.groups()[0]
    return [mutation[0], mutation[1:-1], mutation[-1], mutation]


def parse_line(line):
    row = [v if v != '.' else '\\N' for v in line.strip().split('\t')]
    extra_values = [v.split('=')[-1] for v in row[-1].split(';') if '=' in v]
    mutation_values = format_mutation_missense(extra_values[-2])
    return row + extra_values + mutation_values


def main(infile, outfile):
    if infile.endswith('.gz'):
        ifh = gzip.open(infile, 'rt')
    else:
        ifh = open(infile, 'rt')
    ofh = open(outfile, 'wt')
    ofh.write('\t'.join(header) + '\n')
    for line in ifh:
        if line[0] == '#':
            continue
        row = parse_line(line)
        ofh.write('\t'.join(row) + '\n')
    ifh.close()
    ofh.close()
    print('DONE!')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', type=str)
    parser.add_argument('-o', '--outfile', type=str)
    args = parser.parse_args()
    sys.exit(main(args.infile, args.outfile))
