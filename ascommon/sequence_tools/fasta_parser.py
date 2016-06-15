"""Convert fasta files into tsv.

This file should be kept compatible with Python 2.7 / PyPy!
"""
import sys
import os.path as op
import logging
import gzip
import re

logger = logging.getLogger(__name__)

header = [
    'db',
    'uniprot_acc',
    'uniprot_id',
    'protein_name',
    'organism_name',
    'gene_name',
    'protein_existence',
    'sequence_version',
    'uniprot_sequence',
]

RE_HEADER_LINE = (
    re.compile('^>(sp|tr)\|([A-Z0-9-]+)\|(\w+)(.*?)(OS=.*?)?(GN=.*?)?(PE=.*?)?(SV=.*?)?$')
)


def main(infile, outfile):
    if infile.endswith('.gz'):
        ifh = gzip.open(infile)
    else:
        ifh = open(infile)
    ofh = open(outfile, 'wt')
    row = []
    for line in ifh:
        if line[0] == '>':
            # Header
            try:
                row[-1] = ''.join(row[-1])
            except IndexError:
                logger.info('Writing header...')
                row = header[:]
            ofh.write('\t'.join(row) + '\n')
            try:
                row = parse_header_line(line)
            except AttributeError:
                print(line)
                raise
        else:
            # Sequence
            row[-1].append(line.strip())
    # Last line
    row[-1] = ''.join(row[-1])
    ofh.write('\t'.join(row) + '\n')
    # Cleanup
    ifh.close()
    ofh.close()
    logger.info('DONE!')


def parse_header_line(line):
    """.

    Examples
    --------
    >>> header = (
        '>tr|A0A024R161|A0A024R161_HUMAN Guanine nucleotide-binding protein subunit gamma '
        'OS=Homo sapiens GN=DNAJC25-GNG10 PE=3 SV=1'
    )
    >>> parse_header_line(header) == {
        'db': 'tr',
        'gene_name': 'DNAJC25-GNG10',
        'organism_name': 'Homo sapiens',
        'protein_existence': '3',
        'protein_name': 'Guanine nucleotide-binding protein subunit gamma',
        'sequence_version': '1',
        'uniprot_id': 'A0A024R161',
        'uniprot_name': 'A0A024R161_HUMAN'
    }
    """
    db, uniprot_acc, uniprot_id, protein_name, organism_name, gene_name, protein_existence, \
        sequence_version = RE_HEADER_LINE.match(line).groups()
    row = [
        db.strip(),
        uniprot_acc.strip(),
        uniprot_id.strip(),
        protein_name.strip(),
        organism_name.strip().replace('OS=', '') if organism_name is not None else '\\N',
        gene_name.strip().replace('GN=', '') if gene_name is not None else '\\N',
        protein_existence.strip().replace('PE=', '') if protein_existence is not None else '\\N',
        sequence_version.strip().replace('SV=', '') if sequence_version is not None else '\\N',
        [],
    ]
    return row


if __name__ == '__main__':
    import argparse
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    # logger.setLevel(logging.INFO)
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', type=str, required=True)
    parser.add_argument('-o', '--outfile', type=str, default=None)
    args = parser.parse_args()
    if args.outfile is None:
        name, ext = op.splitext(args.infile)
        if ext in ['.fasta']:
            args.outfile = name + '.tsv'
        else:
            args.outfile = args.infile + '.tsv'
    sys.exit(main(args.infile, args.outfile))
