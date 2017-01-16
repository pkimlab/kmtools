"""Process PDB SIFTS info.

"""
import gzip
import logging
import os
import os.path as op
import re
import tempfile
import urllib.request

import lxml.etree
import numpy as np
import pandas as pd
from lxml import etree

from kmtools import sequence_tools, structure_tools

logger = logging.getLogger(__name__)

CACHE_DIR = None
RE_TAG = re.compile('({.+})?(.+)')
MUTATION_SPLIT_RE = re.compile(' +')


class SIFTSError(Exception):
    pass


def get_cache_dir():
    global CACHE_DIR
    if CACHE_DIR is None:
        CACHE_DIR = op.join(tempfile.gettempdir(), 'sifts_cache')
    return CACHE_DIR


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
                if resname in structure_tools.AAA_DICT:
                    residue_data['pdb_aa'] = structure_tools.AAA_DICT[resname]
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


def get_sifts_data(pdb_id, cache_dict={}):
    """Return SIFTS data for a particular PDB.

    Parameters
    ----------
    pdb_id : str

    cache_dict : dict
        Defaults to memoisation cache.

    Download the xml file for the pdb file with the pdb id pdb_id, parse that
    xml file, and return a dictionry which maps pdb resnumbing to uniprot
    numbering for the chain specified by pdb_chain and uniprot specified by
    uniprot_id.
    """
    if pdb_id in cache_dict:
        return cache_dict[pdb_id]

    cache_dir = get_cache_dir()

    # Download the sifts file if it is not in cache
    sifts_filename = pdb_id.lower() + '.xml.gz'
    if not op.isfile(op.join(cache_dir, sifts_filename)):
        url = 'ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/{}'.format(sifts_filename)
        os.makedirs(cache_dir, exist_ok=True)
        urllib.request.urlretrieve(url, op.join(cache_dir, sifts_filename))

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
    sifts_df = pd.DataFrame(pdb_sifts_data)
    assert len(sifts_df) == len(
        sifts_df.drop_duplicates(subset=['pdb_chain', 'resnum']))

    # residx counts the index of the residue in the pdb, starting from 1
    residx = []
    for pdb_chain in set(sifts_df['pdb_chain']):
        residx.extend(range(1, len(sifts_df[sifts_df['pdb_chain'] == pdb_chain]) + 1))
    sifts_df['residx'] = residx

    # TODO: should optimise the code above instead of simply removing duplicates
    # sifts_df = sifts_df.drop_duplicates()

    assert sifts_df['resnum'].dtype.char == 'O'
    assert sifts_df['residx'].dtype.char in ['l', 'd']

    cache_dict[pdb_id] = sifts_df
    return sifts_df


# #################################################################################################
# Convert PDB to UniProt using SIFTS


def split_mutation(mutation):
    """Split mutation where `wt`, `pos`, and `mut` are separated by spaces.

    Examples
    --------
    >>> split_mutation('A 10 C')
    ('A', '10', 'C')
    >>> split_mutation('A10C')
    ('A', '10', 'C')
    >>> split_mutation(None)
    (None, None, None)
    """
    if pd.isnull(mutation):
        return None, None, None
    if ' ' in mutation:
        mutation_split = MUTATION_SPLIT_RE.split(mutation)
    else:
        mutation_split = (mutation[0], mutation[1:-1], mutation[-1])
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
                mutations_match = str(seqrecord)[int(pos) - 1 + i] == aa
            except (IndexError, ValueError):
                mutations_match = False
            if mutations_match:
                wt_valid += aa
            else:
                wt_valid += '?'
        validated_mutations.append(wt_valid)
    return ','.join(validated_mutations)


def convert_amino_acid(
        pdb_id, pdb_chain, pdb_aa, pdb_resnum, sifts_df, uniprot_id=None, pdb_resnum_offset=0):
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
    sifts_df : DataFrame
        SIFTS data to use for conversion.
    uniprot_id : str
        UniProt ID that is expected.
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
        if uniprot_id is not None and 'uniprot_aa' in sifts_df:
            sifts_df_pdb_aa_match = (sifts_df['uniprot_aa'] == pdb_aa)
        else:
            sifts_df_pdb_aa_match = (sifts_df['pdb_aa'] == pdb_aa)
    else:
        sifts_df_pdb_aa_match = True

    if pd.isnull(pdb_chain):
        sifts_df_pdb_chain_match = True
    else:
        sifts_df_pdb_chain_match = (sifts_df['pdb_chain'] == pdb_chain)

    if uniprot_id is not None and 'uniprot_id' in sifts_df:
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
        if 'uniprot_position' in sifts_df.columns:
            sifts_df_subset_4 = sifts_df[
                (sifts_df['pdb_id'] == pdb_id) &
                (sifts_df_pdb_chain_match) &
                (sifts_df['uniprot_position'] == pdb_residx) &
                (sifts_df_pdb_aa_match) &
                (sifts_df['uniprot_id'] == uniprot_id)]
        else:
            sifts_df_subset_4 = []
    else:
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
    if 'uniprot_position' in sifts_df.columns:
        sifts_df_subset_5 = sifts_df[
            (sifts_df['pdb_id'] == pdb_id) &
            (sifts_df_pdb_chain_match) &
            (sifts_df['uniprot_position'] == pdb_residx) &
            (sifts_df_pdb_aa_match)]
    else:
        sifts_df_subset_5 = []

    # Choose the best availible subset
    if len(sifts_df_subset_0):
        sifts_df_subset = sifts_df_subset_0
    if len(sifts_df_subset_1):
        sifts_df_subset = sifts_df_subset_1
    elif len(sifts_df_subset_2):
        sifts_df_subset = sifts_df_subset_2
    elif len(sifts_df_subset_3):
        sifts_df_subset = sifts_df_subset_3
    elif len(sifts_df_subset_4):
        sifts_df_subset = sifts_df_subset_4
    elif len(sifts_df_subset_5):
        sifts_df_subset = sifts_df_subset_5
    else:
        error_message = """\
SIFTS failed to match residue ({}, {}, {}, {}, {}, {})\
""".format(pdb_id, pdb_chain, pdb_aa, pdb_resnum, uniprot_id, pdb_resnum_offset)
        raise SIFTSError(error_message)

    # Result
    sifts_row = sifts_df_subset.iloc[0]
    uniprot_id_sifts = sifts_row.get('uniprot_id', uniprot_id)
    uniprot_aa = sifts_row.get('uniprot_aa')
    uniprot_pos = sifts_row.get('uniprot_position')
    uniprot_pos = int(uniprot_pos) if pd.notnull(uniprot_pos) else uniprot_pos
    pdb_chain = sifts_row['pdb_chain']
    pfam_id = sifts_row.get('pfam_id', np.nan)

    if pd.notnull(uniprot_aa) and (pdb_aa != uniprot_aa):
        error_message = """\
PDB AA and UniProt AA are not the same! ({}, {}, {}, {}, {}, {}) ({} != {} ({}))\
""".format(pdb_id, pdb_chain, pdb_aa, pdb_resnum, uniprot_id, pdb_resnum_offset,
           pdb_aa, uniprot_aa, uniprot_id_sifts)
        raise SIFTSError(error_message)

    uniprot_seqrecord = sequence_tools.fetch_sequence(uniprot_id_sifts)
    if pd.notnull(uniprot_aa) and pd.notnull(uniprot_pos):
        if not uniprot_seqrecord:
            error_message = """\
    Could not fetch uniprot sequence to make sure that AA '{}{}' is correct!\
    """.format(uniprot_aa, uniprot_pos)
            raise SIFTSError(error_message)
        elif uniprot_seqrecord and (str(uniprot_seqrecord.seq)[uniprot_pos - 1] != uniprot_aa):
            error_message = """\
    Sequence mismatch: {}{} in {}!\
    """.format(uniprot_aa, uniprot_pos, uniprot_seqrecord.seq)
            raise SIFTSError(error_message)

    return dict(
        uniprot_id=uniprot_id_sifts,
        uniprot_aa=uniprot_aa,
        uniprot_pos=uniprot_pos,
        pdb_chain=pdb_chain,
        pfam_id=pfam_id
    )


PDB_MUTATION_RE = re.compile(
    '^([a-zA-Z0-9]{0,1})_?([GVALICMFWPDESTYQNKRH]{1}[1-9][0-9]*[GVALICMFWPDESTYQNKRH]*)$'
)


def format_pdb_mutation(pdb_mutation):
    """
    Examples
    --------
    >>> format_pdb_mutation('A_M1C')
    ('A', 'M1C')
    >>> format_pdb_mutation('AM1C')
    ('A', 'M1C')
    >>> format_pdb_mutation('M1C')
    ('', 'M1C')
    >>> format_pdb_mutation('XXXM1C')
    (None, None)
    """
    match = PDB_MUTATION_RE.match(pdb_mutation)
    if match is None:
        return None, None
    else:
        return match.groups()


def _get_pdb_chain_mutation(pdb_chains, pdb_mutations):
    """
    Examples
    --------
    >>> list(_get_pdb_chain_mutation('A', 'M1C'))
    [('A', 'M1C')]
    >>> list(_get_pdb_chain_mutation('A', 'A_M1C'))
    [('A', 'M1C')]
    >>> list(_get_pdb_chain_mutation(None, 'A_M1C'))
    [('A', 'M1C')]
    >>> list(_get_pdb_chain_mutation('I', 'I_T19L.I_E21K'))
    [('I', 'T19L'), ('I', 'E21K')]
    >>> list(_get_pdb_chain_mutation('B', 'A_M1C'))
    Traceback (most recent call last):
    kmtools.structure_tools.sifts.SIFTSError: PDB chains do not match!
    """
    # sep
    sep = None
    for s in [',', '.']:
        if s in pdb_mutations:
            sep = s
            break
    # pdb_mutations
    if sep is not None:
        pdb_mutations = pdb_mutations.split(sep)
    else:
        pdb_mutations = [pdb_mutations]
    # pdb_chains
    if sep is not None and pd.notnull(pdb_chains) and sep in pdb_chains:
        pdb_chains = pdb_chains.split(sep)
    else:
        pdb_chains = [pdb_chains] * len(pdb_mutations)
    assert len(pdb_chains) == len(pdb_mutations)
    # yield
    for pdb_chain, pdb_mutation in zip(pdb_chains, pdb_mutations):
        c, m = format_pdb_mutation(pdb_mutation)
        if not c:
            yield pdb_chain, m
        else:
            if pdb_chain is not None and pdb_chain != c:
                raise SIFTSError(
                    "pdb_chain ('{}') and c ('{}') do not match!\n({}, {})"
                    .format(pdb_chain, c, pdb_chains, pdb_mutations))
            yield c, m


def convert_pdb_mutations_to_uniprot(
        pdb_id, pdb_chains=None, pdb_mutations=None, *, uniprot_id=None, sifts_df=None):
    """Convert a PDB mutations to UniProt.

    Works for a list of mutations joined with ',' but only one row at a time.

    Fails if can't map at mutation level.

    Parameters
    ----------
    pdb_id : str
    pdb_chains : str | None
        Comma-separated list of PDB chain IDs, with one ID for every mutation

    Returns
    -------
    dict

    Examples
    --------
    >>> from pprint import pprint
    >>> sifts_df = get_sifts_data('1jrh')
    >>> pprint(convert_pdb_mutations_to_uniprot('1jrh', 'I', 'I_T19L.I_E21K', sifts_df=sifts_df))
    {'pdb_mutations_sifts': 'I_T19L.I_E21K',
     'pfam_id_sifts': 'PF01108',
     'uniprot_id_sifts': 'P15260',
     'uniprot_mutations_sifts': 'T36L.E38K'}
    """
    pdb_mutations_sifts = []
    uniprot_mutations_sifts = []
    for pdb_chain, pdb_mutation in _get_pdb_chain_mutation(pdb_chains, pdb_mutations):
        if pd.isnull(pdb_mutation):
            error_message = """\
Failed to get pdb_mutation ({}) for input: ({}, {}, {}, {})\
""".format(pdb_chain, pdb_mutation, pdb_id, pdb_chains, pdb_mutations, uniprot_id)
            raise SIFTSError(error_message)
        pdb_wt, pdb_resnum, pdb_mut = split_mutation(pdb_mutation)
        # Convert each of the mutant residues
        uniprot_aa_data = []
        for i, pdb_wt_aa in enumerate(pdb_wt if pdb_wt else ['']):
            sifts_aa_data = convert_amino_acid(
                pdb_id, pdb_chain, pdb_wt_aa, pdb_resnum, sifts_df,
                uniprot_id=uniprot_id, pdb_resnum_offset=i)
            uniprot_aa_data.append(sifts_aa_data)
        #
        if (any(pd.isnull(x['uniprot_aa']) for x in uniprot_aa_data) or
                any(pd.isnull(x['uniprot_pos']) for x in uniprot_aa_data)):
            error_message = """\
No mutation mapping available! ({}, {}, {}, {}):\n    {}\
""".format(pdb_id, pdb_mutations, pdb_chain, pdb_mutation, uniprot_aa_data)
            raise SIFTSError(error_message)
        uniprot_id_sifts = uniprot_aa_data[0]['uniprot_id']
        pfam_id_sifts = uniprot_aa_data[0]['pfam_id']
        uniprot_wt = ''.join(x['uniprot_aa'] for x in uniprot_aa_data)
        uniprot_pos = uniprot_aa_data[0]['uniprot_pos']
        pdb_chain = uniprot_aa_data[0]['pdb_chain']
        uniprot_mutation = '{}{}{}'.format(uniprot_wt, uniprot_pos, pdb_mut)

        pdb_mutations_sifts.append('{}_{}{}{}'.format(pdb_chain, pdb_wt, pdb_resnum, pdb_mut))
        uniprot_mutations_sifts.append(uniprot_mutation)
        # pdb_chains_sifts = ','.join(x['pdb_chain'] for x in uniprot_aa_data)

    return dict(
        uniprot_id_sifts=uniprot_id_sifts,
        pfam_id_sifts=pfam_id_sifts,
        pdb_mutations_sifts=(
            '.'.join(pdb_mutations_sifts) if pdb_mutations_sifts else np.nan
        ),
        uniprot_mutations_sifts=(
            '.'.join(uniprot_mutations_sifts) if uniprot_mutations_sifts else np.nan
        )
    )


def get_uniprot_id_mutation(pdb_id, pdb_chains, pdb_mutations, uniprot_id):
    """Convert a whole column of PDB mutations to UniProt.

    This function gets the SIFTS data for you.

    Parameters
    ----------
    pdb_id : str
    pdb_mutations : str
    uniprot_id : str

    Returns
    -------
    dict with keys:
        uniprot_id_sifts : str
        uniprot_mutations_sifts : str
            Comma-separated list of mutations.
        pfam_id_sifts : str
            Comma-separated list of cleaned-up PDB mutations.

    Examples
    --------
    >>>
    """
    if pd.isnull(pdb_mutations) or pdb_mutations == '-' or pdb_mutations.startswith('WILD'):
        return {}

    def _uniprot_fallback(uniprot_id, pdb_mutations):
        """Try treating PDB mutations as though they are from UniProt."""
        if pd.isnull(uniprot_id):
            return {}
        else:
            uniprot_seqrecord = sequence_tools.fetch_sequence(uniprot_id)
            if pd.isnull(uniprot_seqrecord):
                print(uniprot_id)
                logger.error(uniprot_id)
            validated_mutations = validate_mutations(pdb_mutations, uniprot_seqrecord)
            return {
                'uniprot_id_sifts': uniprot_id,
                'uniprot_mutations_sifts': validated_mutations,
            }

    # If pdb_id is null, make sure pdb_mutations matches UniProt
    if pd.isnull(pdb_id):
        return _uniprot_fallback(uniprot_id, pdb_mutations)

    # Get SIFTS data
    try:
        sifts_df = get_sifts_data(pdb_id)
    except lxml.etree.XMLSyntaxError:
        return _uniprot_fallback(uniprot_id, pdb_mutations)

    return convert_pdb_mutations_to_uniprot(
        pdb_id, pdb_chains, pdb_mutations, uniprot_id, sifts_df)
