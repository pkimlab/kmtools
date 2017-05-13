import json
import os.path as op
from tempfile import TemporaryDirectory

import pandas as pd
import pytest
from lxml import etree

from kmtools import structure_tools


class TestPresent:
    """Test the case where the SIFTS xml file is already present in the cache directory."""

    @pytest.mark.parametrize('pdb_id', ['1arr', '3mbp', '1jrh'])
    def test_1(self, pdb_id):
        structure_tools.sifts.CACHE_DIR = op.abspath(op.splitext(__file__)[0])
        sifts_data = structure_tools.sifts.get_sifts_data(pdb_id)
        assert isinstance(sifts_data, pd.DataFrame)
        assert not sifts_data.empty

    def test_2(self):
        pdb_id = '1dvf'
        structure_tools.sifts.CACHE_DIR = op.abspath(op.splitext(__file__)[0])
        sifts_data = structure_tools.sifts.get_sifts_data(pdb_id)
        assert isinstance(sifts_data, pd.DataFrame)
        assert not sifts_data.empty
        sifts_data_subset = (
            sifts_data[
                (sifts_data['pdb_id'] == pdb_id) &
                (sifts_data['pdb_chain'] == 'A') &
                (sifts_data['resnum'] == '98')
            ]
        )
        assert not sifts_data_subset.empty
        return sifts_data


class TestAbscent:
    """Test the case where the SIFTS xml file is missing from the cache directory."""

    @classmethod
    def setup_class(cls):
        cls.temp_dir = TemporaryDirectory()
        structure_tools.sifts.CACHE_DIR = cls.temp_dir.name

    @classmethod
    def teardown_class(cls):
        cls.temp_dir.cleanup()

    def test_1(self):
        sifts_data = structure_tools.sifts.get_sifts_data('1arr')
        assert isinstance(sifts_data, pd.DataFrame)
        assert not sifts_data.empty

    def test_2(self):
        sifts_data = structure_tools.sifts.get_sifts_data('3mbp')
        assert isinstance(sifts_data, pd.DataFrame)
        assert not sifts_data.empty
        del sifts_data


@pytest.mark.parametrize("pdb_id,pdb_chains,pdb_mutations", [
    ('2qja', 'C', 'T74A,L78M,G79K,L80Y'),
])
def test_sifts_exception(pdb_id, pdb_chains, pdb_mutations):
    """Test the case where PDB AA and UniProt AA are different."""
    sifts_df = structure_tools.sifts.get_sifts_data(pdb_id)
    with pytest.raises(structure_tools.sifts.SIFTSError):
        structure_tools.sifts.convert_pdb_mutations_to_uniprot(
            pdb_id, pdb_chains, pdb_mutations, sifts_df=sifts_df)


@pytest.fixture(scope='session', params=range(2))
def residue_in_out(request):
    path_template = op.join(op.splitext(__file__)[0], 'residue_data_{}'.format(request.param))
    with open(path_template + '.xml', 'rb') as ifh:
        # Residue entries are wrapped inside <entry></entry> blocks, so take [0]
        xml_data = etree.fromstring(ifh.read())[0]
    with open(path_template + '.json', 'r') as ifh:
        json_data = json.load(ifh)
    return xml_data, json_data


def test__get_residue_info_xml(residue_in_out):
    xml_data, json_data = residue_in_out
    print(xml_data)
    assert structure_tools.sifts._get_residue_info_xml(xml_data) == json_data
