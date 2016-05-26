import os.path as op
from tempfile import TemporaryDirectory
import pytest
from ascommon import system_tools


# %%
urls = [
    'http://files.rcsb.org/download/4DKL.pdb.gz',
    'http://files.rcsb.org/download/4DKL.cif.gz',

    'http://files.rcsb.org/download/1QNV.pdb.gz',
    'http://files.rcsb.org/download/1QNV.cif.gz',

    'http://files.rcsb.org/download/1DSY.pdb.gz',
    'http://files.rcsb.org/download/1DSY.cif.gz',

    'http://files.rcsb.org/download/3TWY.pdb.gz',
    'http://files.rcsb.org/download/3TWY.cif.gz',
]


@pytest.fixture(scope='session', params=urls)
def url(request):
    return request.param


# %%
class Test:

    @classmethod
    def setup_class(cls):
        cls.temp_dir = TemporaryDirectory()

    @classmethod
    def teardown_class(cls):
        cls.temp_dir.cleanup()

    def test(self, url):
        output_file = op.join(self.temp_dir.name, op.basename(url))
        assert not op.isfile(output_file)

        system_tools.download(url, output_file)
        assert op.isfile(output_file)


# %%
if __name__ == '__main__':
    import pytest
    pytest.main([__file__, '-svx'])
