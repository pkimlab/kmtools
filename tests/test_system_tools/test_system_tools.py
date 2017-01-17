import logging
import os.path as op
import shutil
import tempfile

from kmtools import system_tools

logger = logging.getLogger(__name__)


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
