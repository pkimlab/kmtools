"""Tools for dealing with PDB files."""
import os
import os.path as op
import gzip
import re
import logging
from functools import lru_cache
from lxml import etree
import pandas as pd
from .system_tools import download


logger = logging.getLogger(__name__)

RE_TAG = re.compile('({.+})?(.+)')

A_DICT = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'E': 'GLU',
    'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS',
    'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP',
    'Y': 'TYR', 'V': 'VAL', 'U': 'SEC', 'O': 'PYL',
    'B': 'ASX', 'Z': 'GLX', 'J': 'XLE', 'X': 'XAA', '*': 'TER'
}
AAA_DICT = dict([(value, key) for key, value in list(A_DICT.items())])
AAA_DICT['UNK'] = 'X'
AAA_DICT['MSE'] = 'M'
AAA_DICT['CSD'] = 'C'

# Phosphorylated residues
AAA_DICT['SEP'] = 'S'  # PHOSPHOSERINE
AAA_DICT['TPO'] = 'T'  # PHOSPHOTHREONINE
AAA_DICT['SEP'] = 'Y'  # O-PHOSPHOTYROSINE

# Methylated lysines
AAA_DICT['MLZ'] = 'K'
AAA_DICT['MLY'] = 'K'
AAA_DICT['M3L'] = 'K'


class SIFTSError(Exception):
    pass


def iter_residues(xml_data):
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


def get_residue_data(residue):
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
                if resname in AAA_DICT:
                    residue_data['pdb_aa'] = AAA_DICT[resname]
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


@lru_cache(maxsize=512)
def get_sifts_data(pdb_id, cache_dir, cache_dict={}):
    """Return SIFTS data from a particular PDB.

    Download the xml file for the pdb file with the pdb id pdb_id, parse that
    xml file, and return a dictionry which maps pdb resnumbing to uniprot
    numbering for the chain specified by pdb_chain and uniprot specified by
    uniprot_id.
    """
    if pdb_id in cache_dict:
        return cache_dict[pdb_id]

    # Download the sifts file if it is not in cache
    sifts_filename = pdb_id.lower() + '.xml.gz'
    if not os.path.isfile(op.join(cache_dir, sifts_filename)):
        url = 'ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/{}'.format(sifts_filename)
        download(url, op.join(cache_dir, sifts_filename))

    # Go over the xml file and find all cross-references to uniprot
    pdb_sifts_data = []
    with gzip.open(op.join(cache_dir, sifts_filename)) as ifh:
        xml_data = ifh.read()
    for residue in iter_residues(xml_data):
        residue_data = get_residue_data(residue)
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
