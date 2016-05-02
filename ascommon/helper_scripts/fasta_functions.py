#!/usr/bin/python

def split_fasta_file(infile):
    """
    Split fasta files into more manageable chunks

    split_fasta_file('uniprot_sprot.fasta')
    """
    outfile_prefix = infile[:infile.rfind('.')]
    with open(infile, 'r') as ifh:
        output_file_counter = 1
        ofh = open(outfile_prefix + '_' + str(output_file_counter) + '.fasta', 'w')
        fasta_sequence_counter = 0
        for line in ifh:
            if line[0] == '>':
                fasta_sequence_counter += 1
                if fasta_sequence_counter % 50000 == 0:
                    ofh.close()
                    output_file_counter += 1
                    ofh = open(outfile_prefix + '_' + str(output_file_counter) + '.fasta', 'w')
            ofh.writelines(line)
        ofh.close()



def fasta_to_flatfile(infile):
    """
    Convert a fasta file into a tsv file that can be imported into a database

    >>> fasta_to_flatfile('uniprot_sprot.fasta')
    >>> fasta_to_flatfile('uniprot_sprot_varsplic.fasta')
    >>> fasta_to_flatfile('uniprot_trembl.fasta')
    """

    outfile = infile + '.csv'
    with open(outfile, 'w') as ofh:
        column_names = ['db','uniprot_id', 'uniprot_name', 'protein_name', 'organism_name', 'organism_subname', 'gene_name', 'protein_existence', 'sequence_version', 'uniprot_sequence']
        ofh.writelines('\t'.join(column_names) + '\n')
        data = {'uniprot_sequence': '\N'}
        first_sequence = True
        with open(infile, 'r') as ifh:
            for line in ifh:
                if line[0] == '>':
                    if not first_sequence:
                        ofh.writelines('\t'.join([data[x] for x in column_names]) + '\n')
                    data['db'], data['uniprot_id'], data['uniprot_name'] = line[1:line.find(' ')].split('|')
                    desc = line[line.find(' ')+1:].strip()
                    data['protein_name'] = desc[:desc.find(' OS=')].strip()
                    if desc.find(' GN=') != -1:
                        if desc.find(' PE=') != -1:
                            if desc.find(' SV=') != -1:
                                data['organism_name'] = desc[desc.find(' OS=')+4:desc.find(' GN=')].strip()
                                data['gene_name'] = desc[desc.find(' GN=')+4:desc.find(' PE=')].strip()
                                data['protein_existence'] = desc[desc.find(' PE=')+4:desc.find(' SV=')].strip()
                                data['sequence_version'] = desc[desc.find(' SV=')+4:].strip()
                            else:
                                data['organism_name'] = desc[desc.find(' OS=')+4:desc.find(' GN=')].strip()
                                data['gene_name'] = desc[desc.find(' GN=')+4:desc.find(' PE=')].strip()
                                data['protein_existence'] = desc[desc.find(' PE=')+4:].strip()
                                data['sequence_version'] = '\N'
                        else:
                            if desc.find(' SV=') != -1:
                                data['organism_name'] = desc[desc.find(' OS=')+4:desc.find(' GN=')].strip()
                                data['gene_name'] = desc[desc.find(' GN=')+4:desc.find(' SV=')].strip()
                                data['protein_existence'] = '\N'
                                data['sequence_version'] = desc[desc.find(' SV=')+4:].strip()
                            else:
                                data['organism_name'] = desc[desc.find(' OS=')+4:desc.find(' GN=')].strip()
                                data['gene_name'] = desc[desc.find(' GN=')+4:].strip()
                                data['protein_existence'] = '\N'
                                data['sequence_version'] = '\N'
                    else:
                        if desc.find(' PE=') != -1:
                            if desc.find(' SV=') != -1:
                                data['organism_name'] = desc[desc.find(' OS=')+4:desc.find(' PE=')].strip()
                                data['gene_name'] = '\N'
                                data['protein_existence'] = desc[desc.find(' PE=')+4:desc.find(' SV=')].strip()
                                data['sequence_version'] = desc[desc.find(' SV=')+4:].strip()
                            else:
                                data['organism_name'] = desc[desc.find(' OS=')+4:desc.find(' GN=')].strip()
                                data['gene_name'] = '\N'
                                data['protein_existence'] = desc[desc.find(' PE=')+4:].strip()
                                data['sequence_version'] = '\N'
                        else:
                            if desc.find(' SV=') != -1:
                                data['organism_name'] = desc[desc.find(' OS=')+4:desc.find(' SV=')].strip()
                                data['gene_name'] = '\N'
                                data['protein_existence'] = '\N'
                                data['sequence_version'] = desc[desc.find(' SV=')+4:].strip()
                            else:
                                data['organism_name'] = desc[desc.find(' OS=')+4:].strip()
                                data['gene_name'] = '\N'
                                data['protein_existence'] = '\N'
                                data['sequence_version'] = '\N'
                    first_sequence = False
                    data['uniprot_sequence'] = '\N'
                else:
                    data['uniprot_sequence'] = data['uniprot_sequence'] + line.strip()
            ofh.writelines('\t'.join([data[x] for x in column_names]) + '\n')



def select_human(infile, outfile):
    """
    Parse a fasta file and select only sequences that belong to a particular organism

    >>> select_human('uniprot_sprot.fasta', 'uniprot_sprot_human.fasta')
    """
    with open(outfile, 'w') as ofh:
        is_human = False
        with open(infile, 'r') as ifh:
            for line in ifh:
                if line[0] == '>':
                    seq_id = line.split()[0]
                    species = seq_id.split('|')[-1]
                    species = species.split('_')[-1]
                    if species == 'HUMAN':
                        is_human = True
                    else:
                        is_human = False
                if is_human:
                    ofh.writelines(line)
