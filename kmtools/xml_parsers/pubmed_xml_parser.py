import xml.etree.ElementTree as ET
from collections import namedtuple
from typing import List

from kmtools import system_tools

PubmedRow = namedtuple("PubmedRow", [
    'pmid',
    'title',
    'authors',
    'journal',
    'year_published',
    'abstract',
    'mesh_terms',
    'doi',
    'pmc',
])


def parse_pubmed_xml_file(xml_file: str) -> List[PubmedRow]:
    """Parse PubMed XML File."""
    data = []
    for child in _iter_root_children(xml_file):
        row = _process_child(child)
        data.append(row)
    return data


def _iter_root_children(xml_file):
    with system_tools.open_compressed(xml_file) as ifh:
        tree = ET.parse(ifh)
    root = tree.getroot()
    yield from root


def _process_child(child):
    row = PubmedRow(
        _find_pmid(child),
        _find_title(child),
        _find_authors(child),
        _find_journal(child),
        _find_year_published(child),
        _find_abstract(child),
        _find_mesh_terms(child),
        _find_doi(child),
        _find_pmc(child),)
    return row


def _find_pmid(child) -> int:
    pmid = child.find('MedlineCitation/PMID')
    return None if pmid is None else int(pmid.text.strip())


def _find_title(child) -> str:
    title = child.find('MedlineCitation/Article/ArticleTitle')
    return None if title is None else title.text.strip()


def _find_authors(child):
    authors = child.findall('MedlineCitation/Article/AuthorList/Author/LastName')
    return [a.text.strip() for a in authors]


def _find_journal(child):
    journal = child.find('MedlineCitation/Article/Journal/ISOAbbreviation')
    return None if journal is None else journal.text.strip()


def _find_year_published(child):
    year_published = child.find('MedlineCitation/Article/Journal/JournalIssue/PubDate/Year')
    return None if year_published is None else int(year_published.text.strip())


def _find_abstract(child) -> str:
    abstract = child.find('MedlineCitation/Article/Abstract/AbstractText')
    return None if abstract is None else abstract.text.strip()


def _find_mesh_terms(child):
    mesh_terms = child.findall('MedlineCitation/MeshHeadingList/MeshHeading/DescriptorName')
    return [v.text.strip() for v in mesh_terms]


def _find_doi(child) -> str:
    ids = child.findall('PubmedData/ArticleIdList/ArticleId')
    ids = [id_ for id_ in ids if id_.attrib['IdType'] == 'doi']
    return None if not ids else ids[0].text.strip()


def _find_pmc(child) -> str:
    ids = child.findall('PubmedData/ArticleIdList/ArticleId')
    ids = [id_ for id_ in ids if id_.attrib['IdType'] == 'pmc']
    return None if not ids else ids[0].text.strip()
