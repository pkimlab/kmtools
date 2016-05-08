"""Process PDB SIFTS info."""
import os.path as op
import gzip
import re
import logging
from lxml import etree
import numpy as np
import pandas as pd
import lxml.etree
from .. import settings, system_tools, sequence_tools, pdb_tools


logger = logging.getLogger(__name__)

RE_TAG = re.compile('({.+})?(.+)')
MUTATION_SPLIT_RE = re.compile(' +')


class SIFTSError(Exception):
    pass


# #################################################################################################
# Obtain SIFTS data


def _iter_residues_xml(xml_data):
    """Interate over all residue entires in the XML-formatted SIFTS file.

    Parameters
    ----------
    xml_data : bytes
        SIFTS data in XML format.

    Yields
    ------
    residue : lxml.etree._Element
    """
    for entity in etree.fromstring(xml_data):
        # Entries
        if entity.tag.split('}')[-1] == 'entity':
                # Chain segments
                for segment in entity:
                    if segment.tag.split('}')[-1] == 'segment':
                        # Lists of residues
                        for listResidue in segment:
                            if listResidue.tag.split('}')[-1] == 'listResidue':
                                # Residues
                                for residue in listResidue:
                                    if residue.tag.split('}')[-1] == 'residue':
                                        yield residue


def _get_residue_info_xml(residue):
    """Get cross-reference data associated with `residue` element.

    Parameters
    ----------
    residue : lxml.etree._Element

    Returns
    -------
    residue_data : dict
    """
    residue_data = {'is_observed': True, 'comments': []}

    # Go over all database crossreferences keeping only those
    # that come from uniprot and match the given uniprot_id.
    for residue_child in residue:
        residue_child_tag = RE_TAG.findall(residue_child.tag)
        assert len(residue_child_tag) == 1 and len(residue_child_tag[0]) == 2
        residue_child_tag = residue_child_tag[0][1]
        # Some more details about the residue
        if residue_child_tag == 'residueDetail':
            residue_data['comments'].append(residue_child.text)
            if residue_child.text == 'Not_Observed':
                residue_data['is_observed'] = False
        # Mappings to other databases
        if residue_child_tag == 'crossRefDb':
            if residue_child.attrib.get('dbSource') == 'PDB':
                residue_data['pdb_id'] = residue_child.attrib.get('dbAccessionId')
                residue_data['pdb_chain'] = residue_child.attrib.get('dbChainId')
                residue_data['resnum'] = residue_child.attrib.get('dbResNum')
                resname = residue_child.attrib.get('dbResName')
                if resname in pdb_tools.AAA_DICT:
                    residue_data['pdb_aa'] = pdb_tools.AAA_DICT[resname]
                else:
                    logger.warning(
                        'Could not convert amino acid {} to a one letter code!'
                        .format(resname))
                    residue_data['pdb_aa'] = resname
            elif residue_child.attrib.get('dbSource') == 'UniProt':
                residue_data['uniprot_id'] = residue_child.attrib.get('dbAccessionId')
                residue_data['uniprot_position'] = residue_child.attrib.get('dbResNum')
                residue_data['uniprot_aa'] = residue_child.attrib.get('dbResName')
            elif residue_child.attrib.get('dbSource') == 'Pfam':
                residue_data['pfam_id'] = residue_child.attrib.get('dbAccessionId')

    residue_data['comments'] = ','.join(residue_data['comments'])

    return residue_data


def _get_sifts_data(pdb_id, cache_dir=None, cache_dict={}):
    """Return SIFTS data for a particular PDB.

    Parameters
    ----------
    cache_dict : dict
        Defaults to memoisation cache.

    Download the xml file for the pdb file with the pdb id pdb_id, parse that
    xml file, and return a dictionry which maps pdb resnumbing to uniprot
    numbering for the chain specified by pdb_chain and uniprot specified by
    uniprot_id.
    """
    cache_dir = settings.get_cache_dir(cache_dir)
    if pdb_id in cache_dict:
        return cache_dict[pdb_id]

    # Download the sifts file if it is not in cache
    sifts_filename = pdb_id.lower() + '.xml.gz'
    if not op.isfile(op.join(cache_dir, sifts_filename)):
        url = 'ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/{}'.format(sifts_filename)
        system_tools.download(url, op.join(cache_dir, sifts_filename))

    # Go over the xml file and find all cross-references to uniprot
    pdb_sifts_data = []
    with gzip.open(op.join(cache_dir, sifts_filename)) as ifh:
        xml_data = ifh.read()
    for residue in _iter_residues_xml(xml_data):
        residue_data = _get_residue_info_xml(residue)
        if residue_data is None:
            continue
        pdb_sifts_data.append(residue_data)

    # Convert data to a DataFrame and make sure we have no duplicates
    pdb_sifts_data_df = pd.DataFrame(pdb_sifts_data)
    assert len(pdb_sifts_data_df) == len(
        pdb_sifts_data_df.drop_duplicates(subset=['pdb_chain', 'resnum']))

    # TODO: should optimise the code above instead of simply removing duplicates
    # pdb_sifts_data_df = pdb_sifts_data_df.drop_duplicates()

    cache_dict[pdb_id] = pdb_sifts_data_df
    return pdb_sifts_data_df


def get_sifts_data(pdb_id, pdb_mutations, cache_dir=None, cache_dict={}):
    """Wrapper around _get_sifts_data."""
    cache_dir = settings.get_cache_dir(cache_dir)
    try:
        sifts_df = cache_dict[pdb_id]
    except KeyError:
        sifts_df = _get_sifts_data(pdb_id, cache_dir)
        cache_dict[pdb_id] = sifts_df

    # Sometimes SIFTS does not contain any uniprot information for the protein in question
    subset_columns = ['uniprot_position', 'uniprot_aa', 'pdb_chain']
    missing_columns = sorted(set(subset_columns) - set(sifts_df.columns))
    if missing_columns == ['uniprot_aa', 'uniprot_position']:
        raise SIFTSError('SIFTS has no UniProt annotation for this protein')
    elif missing_columns:
        raise SIFTSError("SIFTS information missing: {}".format(missing_columns))
    sifts_df = sifts_df.dropna(subset=subset_columns)

    # residx counts the index of the residue in the pdb, starting from 1
    residx = []
    for pdb_chain in set(sifts_df['pdb_chain']):
        residx.extend(range(1, len(sifts_df[sifts_df['pdb_chain'] == pdb_chain]) + 1))
    sifts_df['residx'] = residx

    assert sifts_df['resnum'].dtype.char == 'O'
    assert sifts_df['residx'].dtype.char in ['l', 'd']
    return sifts_df


# #################################################################################################
# Convert PDB to UniProt using SIFTS


def split_mutation(mutation):
    """Split mutation where `wt`, `pos`, and `mut` are separated by spaces.

    Examples
    --------
    >>> split_mutation('A 10 C')
    ('A', '10', 'C')

    """
    mutation_split = MUTATION_SPLIT_RE.split(mutation)
    if len(mutation_split) == 3:
        wt, pos, mut = mutation_split
    elif len(mutation_split) == 2:
        if any(c.isdigit() for c in mutation_split[0]):
            wt, pos, mut = '', *mutation_split
        elif any(c.isdigit() for c in mutation_split[1]):
            wt, pos, mut = *mutation_split, ''
        else:
            raise Exception(mutation)
    else:
        raise Exception(mutation)
    return wt, pos, mut


def validate_mutations(mutations, seqrecord):
    """Make sure that mutations match seqrecord.

    Returns
    -------
    validated_mutations : str
        Comma-separated list of mutations where each unmatched mutation has been replaced by a '?'.
    """
    validated_mutations = []
    for mutation in mutations.split(','):
        wt, pos, mut = split_mutation(mutation)
        wt_valid = ''
        for i, aa in enumerate(wt):
            try:
                mutations_match = str(seqrecord.seq)[int(pos) - 1 + i] == aa
            except (IndexError, ValueError):
                mutations_match = False
            if mutations_match:
                wt_valid += aa
            else:
                wt_valid += '?'
        validated_mutations.append(wt_valid)
    return ','.join(validated_mutations)


def convert_amino_acid(
        pdb_id, pdb_chain, pdb_aa, pdb_resnum, uniprot_id, sifts_df, pdb_resnum_offset=0):
    """Convert a single amino acid from PDB to UniProt coordinates.

    Parameters
    ----------
    pdb_id : str
        PDB ID.
    pdb_chain : str | None
        PDB chain of the residue to be converted.
    pdb_aa : str | ''
        PDB amino acid of the residue to be converted.
    pdb_resnum : str
        PDB RESNUM of the residue to be converted.
    uniprot_id : str
        UniProt ID that is expected.
    sifts_df : DataFrame
        SIFTS data to use for conversion.
    pdb_resnum_offset : int
        Move `pdb_resnum` forward or backward by a certain number of residues.

    Returns
    -------
    dict
    """
    if pdb_resnum_offset:
        pdb_resnum_idx = sifts_df[sifts_df['resnum'] == pdb_resnum].index[0]
        pdb_resnum = sifts_df.loc[pdb_resnum_idx + pdb_resnum_offset, 'resnum']

    pdb_residx = int(''.join(c for c in pdb_resnum if c.isdigit()))

    if pdb_aa:
        if 'uniprot_aa' in sifts_df:
            sifts_df_pdb_aa_match = (sifts_df['uniprot_aa'] == pdb_aa)
        else:
            sifts_df_pdb_aa_match = (sifts_df['pdb_aa'] == pdb_aa)
    else:
        sifts_df_pdb_aa_match = True

    if pdb_chain is None:
        sifts_df_pdb_chain_match = True
    else:
        sifts_df_pdb_chain_match = (sifts_df['pdb_chain'] == pdb_chain)

    try:
        # Get the subset of rows that we are interested in
        sifts_df_subset_0 = sifts_df[
            (sifts_df['pdb_id'] == pdb_id) &
            (sifts_df_pdb_chain_match) &
            (sifts_df['resnum'] == pdb_resnum) &
            (sifts_df_pdb_aa_match) &
            (sifts_df['uniprot_id'] == uniprot_id)]
        # Try using residx instead of resnum
        sifts_df_subset_2 = sifts_df[
            (sifts_df['pdb_id'] == pdb_id) &
            (sifts_df_pdb_chain_match) &
            (sifts_df['residx'] == pdb_residx) &
            (sifts_df_pdb_aa_match) &
            (sifts_df['uniprot_id'] == uniprot_id)]
    except KeyError as e:
        print(e)
        sifts_df_subset_0 = []
        sifts_df_subset_2 = []

    # Try mapping to a wildcard uniprot
    sifts_df_subset_1 = sifts_df[
        (sifts_df['pdb_id'] == pdb_id) &
        (sifts_df_pdb_chain_match) &
        (sifts_df['resnum'] == pdb_resnum) &
        (sifts_df_pdb_aa_match)]
    # Or residx and wildcard uniprot
    sifts_df_subset_3 = sifts_df[
        (sifts_df['pdb_id'] == pdb_id) &
        (sifts_df_pdb_chain_match) &
        (sifts_df['residx'] == pdb_residx) &
        (sifts_df_pdb_aa_match)]

    # Choose the best availible subset
    if len(sifts_df_subset_0):
        sifts_df_subset = sifts_df_subset_0
    if len(sifts_df_subset_1):
        sifts_df_subset = sifts_df_subset_1
    elif len(sifts_df_subset_2):
        sifts_df_subset = sifts_df_subset_2
    elif len(sifts_df_subset_3):
        sifts_df_subset = sifts_df_subset_3
    else:
        error_message = """\
SIFTS failed to match residue ({}, {}, {}, {})\
""".format(pdb_id, pdb_aa, pdb_resnum, uniprot_id)
        raise SIFTSError(error_message)

    # Result
    uniprot_id = sifts_df_subset.iloc[0].get('uniprot_id', uniprot_id)
    uniprot_aa = sifts_df_subset.iloc[0]['uniprot_aa']
    uniprot_pos = int(sifts_df_subset.iloc[0]['uniprot_position'])
    pdb_chain = sifts_df_subset.iloc[0]['pdb_chain']

    uniprot_seqrecord = sequence_tools.get_uniprot_sequence(uniprot_id)
    if str(uniprot_seqrecord.seq)[uniprot_pos - 1] != uniprot_aa:
        uniprot_aa = '?'
        uniprot_pos = np.nan

    return dict(
        uniprot_id=uniprot_id, uniprot_aa=uniprot_aa, uniprot_pos=uniprot_pos, pdb_chain=pdb_chain
    )


def convert_pdb_mutations_to_uniprot(pdb_id, pdb_chains, pdb_mutations, uniprot_id, sifts_df):
    """Convert a PDB mutations to UniProt.

    Works for a list of mutations joined with ',' but only one row at a time.

    Parameters
    ----------
    pdb_id : str
    pdb_chains : str | None
        Comma-separated list of PDB chain IDs, with one ID for every mutation

    Returns
    -------
    uniprot_id_out : str
        `uniprot_id` extracted from SIFTS.
    uniprot_mutations : str
        Comma-separated list of mutations in UniProt coordinates.
    """
    uniprot_mutations = []

    def get_pdb_chain_mutation():
        if ',' not in pdb_mutations:
            yield pdb_chains, pdb_mutations
        elif pd.isnull(pdb_chains) or ',' not in pdb_chains:
            for pdb_mutation in pdb_mutations.split(','):
                yield pdb_chains, pdb_mutation
        else:
            yield from zip(pdb_chains.split(','), pdb_mutations.split(','))

    for pdb_chain, pdb_mutation in get_pdb_chain_mutation():
        pdb_wt, pdb_resnum, pdb_mut = split_mutation(pdb_mutation)
        # Convert each of the mutant residues
        uniprot_aa_data = []
        for i, pdb_wt_aa in enumerate(pdb_wt if pdb_wt else ['']):
            uniprot_aa_data.append(
                convert_amino_acid(
                    pdb_id, pdb_chain, pdb_wt_aa, pdb_resnum, uniprot_id, sifts_df,
                    pdb_resnum_offset=i))
        uniprot_id = uniprot_aa_data[0]['uniprot_id']
        pdb_chains = ','.join(x['pdb_chain'] for x in uniprot_aa_data)
        uniprot_wt = ''.join(x['uniprot_aa'] for x in uniprot_aa_data)
        uniprot_pos = uniprot_aa_data[0]['uniprot_pos']
        uniprot_mutation = '{}{}{}'.format(uniprot_wt, uniprot_pos, pdb_mut)
        uniprot_mutations.append(uniprot_mutation)

    return uniprot_id, ','.join(uniprot_mutations), pdb_chains


def get_uniprot_id_mutation(pdb_id, pdb_chains, pdb_mutations, uniprot_id):
    """Convert a whole column of PDB mutations to UniProt.

    Parameters
    ----------
    pdb_id : str
    pdb_mutations : str
    uniprot_id : str

    Returns
    -------
    uniprot_id : str
    uniprot_mutations : str
        Comma-separated list of mutations.
    pdb_mutations : str
        Comma-separated list of cleaned-up PDB mutations.
    """
    if pd.isnull(pdb_mutations) or pdb_mutations == '-' or pdb_mutations.startswith('WILD'):
        return uniprot_id, '-', np.nan, np.nan

    def _uniprot_fallback(uniprot_id, pdb_mutations):
        """Try treating PDB mutations as though they are from UniProt."""
        if pd.isnull(uniprot_id):
            return np.nan, np.nan, np.nan, np.nan
        else:
            uniprot_seqrecord = sequence_tools.get_uniprot_sequence(uniprot_id)
            validated_mutations = validate_mutations(pdb_mutations, uniprot_seqrecord)
            return uniprot_id, validated_mutations, np.nan, np.nan

    # If pdb_id is null, make sure pdb_mutations matches UniProt
    if pd.isnull(pdb_id):
        return _uniprot_fallback(uniprot_id, pdb_mutations)

    # Get SIFTS data
    try:
        sifts_df = get_sifts_data(pdb_id, pdb_mutations)
    except lxml.etree.XMLSyntaxError:
        return _uniprot_fallback(uniprot_id, pdb_mutations)

    uniprot_id, uniprot_mutations, pdb_chains = convert_pdb_mutations_to_uniprot(
        pdb_id, pdb_chains, pdb_mutations, uniprot_id, sifts_df)
    return uniprot_id, uniprot_mutations, pdb_chains, pdb_mutations
