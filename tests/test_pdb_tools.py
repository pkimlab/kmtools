import os.path as op
import json
from lxml import etree
from tempfile import TemporaryDirectory
import pandas as pd
import pytest
from ascommon import pdb_tools


class TestPresent:
    """Test the case where the SIFTS xml file is already present in the cache directory."""

    def test_1(self):
        sifts_data = pdb_tools.get_sifts_data('1arr', op.splitext(__file__)[0])
        assert isinstance(sifts_data, pd.DataFrame)
        assert not sifts_data.empty

    def test_2(self):
        sifts_data = pdb_tools.get_sifts_data('3mbp', op.splitext(__file__)[0])
        assert isinstance(sifts_data, pd.DataFrame)
        assert not sifts_data.empty


class TestAbscent:
    """Test the case where the SIFTS xml file is missing from the cache directory."""

    @classmethod
    def setup_class(cls):
        cls.temp_dir = TemporaryDirectory()

    @classmethod
    def teardown_class(cls):
        cls.temp_dir.cleanup()

    def test_1(self):
        sifts_data = pdb_tools.get_sifts_data('1arr', self.temp_dir)
        assert isinstance(sifts_data, pd.DataFrame)
        assert not sifts_data.empty

    def test_2(self):
        sifts_data = pdb_tools.get_sifts_data('3mbp', self.temp_dir)
        assert isinstance(sifts_data, pd.DataFrame)
        assert not sifts_data.empty
        del sifts_data


@pytest.fixture(scope='session', params=range(2))
def residue_in_out(request):
    path_template = op.join(op.splitext(__file__)[0], 'residue_data_{}'.format(request.param))
    with open(path_template + '.xml', 'rb') as ifh:
        # Residue entries are wrapped inside <entry></entry> blocks, so take [0]
        xml_data = etree.fromstring(ifh.read())[0]
    with open(path_template + '.json', 'r') as ifh:
        json_data = json.load(ifh)
    return xml_data, json_data


def test_get_residue_data(residue_in_out):
    xml_data, json_data = residue_in_out
    print(xml_data)
    assert pdb_tools.get_residue_data(xml_data) == json_data


if __name__ == '__main__':
    import pytest
    pytest.main([__file__, '-svx'])
