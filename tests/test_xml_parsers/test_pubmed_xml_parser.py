import logging
import os.path as op

import pandas as pd
import pytest
import yaml

from kmtools.xml_parsers import pubmed_xml_parser

logger = logging.getLogger(__name__)


def _iter_test_data(basename):
    basename = op.join(op.splitext(__file__)[0], basename)
    xml_children = pubmed_xml_parser._iter_root_children(basename + '.xml')
    with open(op.join(op.dirname(op.splitext(__file__)[0]), basename + '.yaml')) as ifh:
        xml_data = yaml.load(ifh)
    yield from zip(xml_children, xml_data)


@pytest.mark.parametrize("attribute, xml_child, xml_data_",
                         [(attribute, xml_child, xml_data[attribute])
                          for (xml_child, xml_data) in _iter_test_data('test_data')
                          for attribute in xml_data])
def test_finders(attribute, xml_child, xml_data_):
    fn = getattr(pubmed_xml_parser, "_find_{}".format(attribute))
    xml_data = fn(xml_child)
    assert xml_data == xml_data_


@pytest.mark.parametrize("xml_file", [op.join(op.splitext(__file__)[0], 'medline_sample.xml.gz')])
def test_parse_pubmed_xml_file(xml_file):
    data = pubmed_xml_parser.parse_pubmed_xml_file(xml_file)
    df = pd.DataFrame(data)
    assert df.doi.notnull().sum() >= 72
    assert df.abstract.notnull().sum() >= 967
    assert df.shape == (1861, 8)
    logger.info(df.tail())
    logger.info(df.year_published.isnull().sum())
