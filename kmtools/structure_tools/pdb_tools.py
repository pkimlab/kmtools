import os.path as op
import io
import logging
import tempfile
import subprocess

import Bio.PDB

from kmtools import system_tools

# __all__ = ['fetch_structure', 'load_structure']

logger = logging.getLogger(__name__)


def _guess_pdb_id(pdb_file):
    """Extract the PDB id from a PDB file.

    Examples
    --------
    >>> _guess_pdb_id('4dkl.pdb')
    '4dkl'
    >>> _guess_pdb_id('/data/structures/divided/pdb/26/pdb126d.ent.gz')
    '126d'
    >>> _guess_pdb_id('/tmp/100d.cif.gz')
    '100d'
    """
    pdb_id = op.basename(pdb_file)
    pdb_id = system_tools.remove_extensions(pdb_id, ['.gz', '.pdb', '.ent', '.cif'])
    if len(pdb_id) == 7 and (pdb_id.startswith('ent') or pdb_id.startswith('pdb')):
        pdb_id = pdb_id[3:]
        assert len(pdb_id) == 4
    pdb_id = pdb_id.lower()
    pdb_id = pdb_id.replace('.', '')
    return pdb_id


def _guess_pdb_type(pdb_file):
    """Guess PDB file type from file name.

    Examples
    >>> _guess_pdb_type('4dkl.pdb')
    'pdb'
    >>> _guess_pdb_type('/tmp/4dkl.cif.gz')
    'cif'
    """
    file_name, file_ext = op.splitext(pdb_file)
    while file_ext:
        file_ext = file_ext.lower()
        if file_ext in ['.pdb']:
            return 'pdb'
        if file_ext in ['.cif', '.mmcif']:
            return 'cif'
        file_name, file_ext = op.splitext(file_name)
    return None


def get_wwpdb_url(
        pdb_id, pdb_type='cif', biounit=False, wwpdb_base_url=None):
    """Get PDB filename (inlcuding path) from a local mirror of the PDB repository.

    This assumes that you have mirrored the PDB ftp repository:
    ftp://ftp.wwpdb.org/pub/pdb/data/

    Examples
    --------
    >>> get_wwpdb_url('4DKL', 'pdb', biounit=False, wwpdb_base_url='')
    'structures/divided/pdb/dk/pdb4dkl.ent.gz'
    >>> get_wwpdb_url('4dkl', 'cif', biounit=False, wwpdb_base_url='')
    'structures/divided/mmCIF/dk/4dkl.cif.gz'
    >>> get_wwpdb_url('4DKL', 'pdb', biounit=True, wwpdb_base_url='/tmp')
    '/tmp/biounit/PDB/divided/dk/4dkl.pdb1.gz'
    >>> get_wwpdb_url('4NWR', 'cif', biounit=True)
    'ftp://ftp.wwpdb.org/pub/pdb/data/biounit/mmCIF/divided/nw/4nwr-assembly1.cif.gz'
    """
    if wwpdb_base_url is None:
        wwpdb_base_url = 'ftp://ftp.wwpdb.org/pub/pdb/data/'

    if pdb_type.startswith('pdb') and not biounit:
        subfolder = op.join('structures', 'divided', 'pdb')
        filename = 'pdb{}.ent.gz'.format(pdb_id.lower())
    elif pdb_type.startswith('pdb') and biounit:
        subfolder = op.join('biounit', 'PDB', 'divided')
        filename = '{}.pdb1.gz'.format(pdb_id.lower())
    elif pdb_type.startswith('cif') and not biounit:
        subfolder = op.join('structures', 'divided', 'mmCIF')
        filename = '{}.cif.gz'.format(pdb_id.lower())
    elif pdb_type.startswith('cif') and biounit:
        subfolder = op.join('biounit', 'mmCIF', 'divided')
        filename = '{}-assembly1.cif.gz'.format(pdb_id.lower())
    else:
        raise Exception

    pdb_filename = op.join(wwpdb_base_url, subfolder, pdb_id[1:3].lower(), filename)
    return pdb_filename


def _get_pdb_url(pdb_id, pdb_type='cif', mirror='rcsb'):
    """Download PDB from RCSB or EBI website.

    .. todo:: Add an option for biounit urls.

    .. note:: Not used at the moment.
    """
    if mirror == 'rcsb':
        return 'http://www.rcsb.org/pdb/files/{}.pdb'.format(pdb_id.lower())
    elif mirror == 'ebi':
        return 'http://www.ebi.ac.uk/pdbe/entry-files/download/pdb{}.ent'.format(pdb_id.lower())
    else:
        raise Exception("Unsupported mirror: {}.".format(mirror))


def _get_pdb_parser(pdb_type):
    """Get BioPython PDB parser appropriate for `pdb_type`."""
    if pdb_type == 'pdb':
        return Bio.PDB.PDBParser()
    elif pdb_type == 'cif':
        return Bio.PDB.MMCIFParser()
    else:
        raise Exception


def _gen_assembly(data):
    with tempfile.TemporaryDirectory() as tmpdirname:
        cif_file = op.join(tmpdirname, 'xxxx.cif')
        with open(cif_file, 'wb') as ofh:
            ofh.write(data)
        subprocess.run(['Assemblies', cif_file], cwd=op.dirname(cif_file))
        with open(cif_file.replace('.cif', '-1.cif'), 'rb') as ifh:
            data = ifh.read()
    return data


def fetch_structure(pdb_id, pdb_type='cif', biounit=False, pdb_mirror=None):
    """Fetch remote PDB file.

    .. warning::

       For `cif` biounits, we have to run the ``Assembly`` binary
       to generate the biounit structure.

       Even for `pdb` biounits, some structures don't have a biounit file
       (probably because they are NMR structures, etc.)

    Returns
    -------
    :class:`Bio.PDB.Structure.Structure`
        Protein structure.

    Examples
    --------
    >>> fetch_structure('4dkl')
    Structure(id='4dkl')
    >>> fetch_structure('4dkl', pdb_type='cif')
    Structure(id='4dkl')
    >>> fetch_structure('4NWR', pdb_type='cif', biounit=True)
    Structure(id='4NWR')
    """
    gen_assembly = False
    if pdb_type == 'cif' and biounit:
        biounit = False
        gen_assembly = True

    url = get_wwpdb_url(pdb_id, pdb_type=pdb_type, biounit=biounit, wwpdb_base_url=pdb_mirror)
    logger.debug("url: %s", url)

    pdb_data = system_tools.read_url(url)

    if gen_assembly:
        pdb_data = _gen_assembly(pdb_data)

    parser = _get_pdb_parser(pdb_type)
    structure = parser.get_structure(pdb_id, io.StringIO(pdb_data.decode('utf-8')))

    return structure


def load_structure(pdb_file, pdb_id=None, pdb_type=None):
    """Load local PDB file.

    Parameters
    ----------
    pdb_file : :obj:`str`
        File to load
    pdb_id : :obj:`str`, optional
        PDB_ID to assign to the loaded structure.
    pdb_type: :obj:`str`, one of: `{'pdb', 'cif'}`, optional
        Type of PDB file. If ``None``, try all available types until one works.

    Returns
    -------
    :class:`Bio.PDB.Structure.Structure`
        Protein structure.

    Examples
    --------
    >>> import urllib.request

    >>> pdb_file = op.join(tempfile.gettempdir(), '4dkl.pdb')
    >>> r = urllib.request.urlretrieve('https://files.rcsb.org/download/4dkl.pdb', pdb_file)
    >>> load_structure(pdb_file)
    Structure(id='4dkl')
    >>> os.remove(pdb_file)

    >>> pdb_file = op.join(tempfile.gettempdir(), '3K1Q.cif')
    >>> r = urllib.request.urlretrieve('https://files.rcsb.org/download/3K1Q.cif', pdb_file)
    >>> load_structure(pdb_file)
    Structure(id='3k1q')
    >>> os.remove(pdb_file)
    """
    if pdb_id is None:
        pdb_id = _guess_pdb_id(pdb_file)

    if pdb_type is not None:
        pdb_types = [pdb_type]
    else:
        pdb_types = ['cif', 'pdb']

    with system_tools.open_compressed(pdb_file, mode='rt') as ifh:
        for parser in (_get_pdb_parser(pdb_type) for pdb_type in pdb_types):
            try:
                structure = parser.get_structure(pdb_id, ifh)
            except KeyError:
                logger.info("Count not load structure using the %s parser.", parser)

    return structure
