""".

When using `xml.etree`, you have to clear the each element as you walk along the tree:
http://www.ibm.com/developerworks/xml/library/x-hiperfparse/
"""
import os.path as op
import gzip
import xml.etree.ElementTree as ET


def multi_open(filepath, *args, **kwargs):
    extension = op.splitext(filepath)[-1]
    if extension == 'gz':
        return gzip.open(*args, **kwargs)
    elif extension == 'bz2':
        raise NotImplementedError
    else:
        return open(*args, **kwargs)


def parse_uniparc_xml(filepath):
    ifh = multi_open(filepath)
    for context in ET.iterparse(ifh, events=('start', 'end')):
        pass
    ifh.close()
