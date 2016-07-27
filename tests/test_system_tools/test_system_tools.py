import os.path as op
from tempfile import TemporaryDirectory
import pytest
import tempfile
import shutil
import logging
from kmtools import system_tools

logger = logging.getLogger(__name__)


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


class TestOpenExclusively:

    @classmethod
    def setup_class(cls):
        cls.tempdir = tempfile.mkdtemp()
        cls.tempfile = op.join(cls.tempdir, 'deleteme.txt')

    def teardown_class(cls):
        shutil.rmtree(cls.tempdir)

    def _worker(self, x):
        with system_tools.open_exclusively(self.tempfile) as ofh:
            ofh.write('\t'.join(str(x) for i in range(10)) + '\n')

    def test_open_exclusively(self):
        logger.info("Running test_open_exclusively({})".format(self))
        # Set up
        NUMBER_OF_WRITES = 5000
        # Run
        from multiprocessing import Pool
        with Pool(processes=36) as pool:
            pool.map(self._worker, range(NUMBER_OF_WRITES))
        # Results
        with open(self.tempfile, 'r') as ofh:
            data = ofh.readlines()
        assert len(data) == NUMBER_OF_WRITES
