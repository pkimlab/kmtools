import logging
import os
import os.path as op
import shutil
import tempfile

import pytest

from kmtools import system_tools

logger = logging.getLogger(__name__)


@pytest.mark.parametrize("host", [None, "127.0.0.1"])
def test_run_command(host):
    """
    Used by `odbo`.
    """
    root_files_ = sorted(os.listdir(op.expanduser("~")))
    stdout, stderr, returncode = system_tools.run_command("ls -A1 ~")
    assert stderr == ""
    assert returncode == 0
    root_files = sorted(stdout.split())
    assert root_files == root_files_


class TestOpenExclusively:
    @classmethod
    def setup_class(cls):
        cls.tempdir = tempfile.mkdtemp()
        cls.tempfile = op.join(cls.tempdir, "deleteme.txt")

    def teardown_class(cls):
        shutil.rmtree(cls.tempdir)

    def _worker(self, x):
        with system_tools.open_exclusively(self.tempfile) as ofh:
            ofh.write("\t".join(str(x) for i in range(10)) + "\n")

    def test_open_exclusively(self):
        logger.info("Running test_open_exclusively({})".format(self))
        # Set up
        NUMBER_OF_WRITES = 5000
        # Run
        from multiprocessing import Pool

        with Pool(processes=36) as pool:
            pool.map(self._worker, range(NUMBER_OF_WRITES))
        # Results
        with open(self.tempfile, "r") as ofh:
            data = ofh.readlines()
        assert len(data) == NUMBER_OF_WRITES
