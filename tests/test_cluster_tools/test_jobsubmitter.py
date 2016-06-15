""".

.. note::
    Tests in this file require access to compute clusters on hard-coded IPs.
"""
import sys
import os
import os.path as op
import logging
import tempfile
import time
import json
import paramiko
import pytest
import tarfile
import shutil
from collections import Counter
import ascommon
import datapkg

ascommon.py_tools.configure_logging('debug')
logger = logging.getLogger(__name__)
logger.info('Hello world')


def _parse_connection_string(connection_string):
    _db_info = datapkg.parse_connection_string(connection_string)
    return _db_info['db_type'], _db_info['host_ip']


def _test_ssh_connection(connection_string):
    _, ip = _parse_connection_string(connection_string)
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    try:
        ssh.connect(ip)
        return True
    except Exception as e:
        logger.error("{}: '{}'".format(type(e), str(e)))
        return False


test_input = [
    (connection_string, concurrent_job_limit)
    for connection_string in [('sge://:@192.168.6.201'), ('pbs://:@192.168.233.150')]
    for concurrent_job_limit in [None, 50]
    if _test_ssh_connection(connection_string)
]
logger.info("Collected {} test inputs".format(len(test_input)))


def get_iterable(script_filename):
    return [
        (str(job_id), '{python} "{script}" -i {job_id}'.format(
            python=sys.executable,
            script=op.join(op.abspath(op.dirname(__file__)), script_filename),
            job_id=job_id))
        for job_id in range(99)
    ]


@pytest.mark.skipif(pytest.config.getvalue("quick"), reason="Tests take several minutes.")
@pytest.mark.parametrize("connection_string, concurrent_job_limit", test_input)
def test_1(connection_string, concurrent_job_limit):
    """Test on tasks that finish successfully."""
    job_name = 'test_1'
    iterable = get_iterable('_test_1.py')
    # tempdir = tempfile.TemporaryDirectory(dir=op.expanduser('~/tmp'))
    # lrp = tempdir.name
    lrp = tempfile.mkdtemp(dir=op.expanduser('~/tmp'))
    # Submit jobs
    js = ascommon.cluster_tools.JobSubmitter(
        job_name, connection_string, lrp,
        walltime='01:00:00',
        concurrent_job_limit=concurrent_job_limit,
        queue='short',
        priority="0")
    logger.info('Submitting...')
    with js.connect():
        js.submit(iterable)
    logger.info('Finished submitting...')
    # Make sure that jobs finish successfully
    time_0 = time.time()
    for job_id, system_command in iterable:
        while True:
            stderr_file = op.join(js.log_path, '{}.err'.format(job_id))
            try:
                with open(stderr_file, 'rt') as ifh:
                    stderr = ifh.read().strip()
            except FileNotFoundError:
                logger.debug("File '{}' not found!".format(stderr_file))
                time.sleep(5)
                continue
            if stderr.endswith('DONE!'):
                stdout_file = op.join(js.log_path, '{}.out'.format(job_id))
                with open(stdout_file, 'rt') as ifh:
                    stdout = ifh.read().strip()
                try:
                    stdout_dict = json.loads(stdout)
                except json.JSONDecodeError:
                    logger.error(stdout)
                    assert False, "Bad output"
                assert stdout_dict['job_id'] == job_id
                break
            elif stderr.endswith('ERROR!'):
                logger.error(stderr)
                assert False, (job_id, system_command)
            else:
                logger.debug('Waiting for job to finish...')
                time.sleep(5)
                if (time.time() - time_0) <= 10 * 60:
                    continue
                else:
                    logger.error(stderr)
                    assert False, "Timeout!"


@pytest.mark.skipif(pytest.config.getvalue("quick"), reason="Tests take several minutes.")
@pytest.mark.parametrize("connection_string, concurrent_job_limit", test_input)
def test_2(connection_string, concurrent_job_limit):
    """Test on tasks that finish crash."""
    job_name = 'test_2'
    iterable = get_iterable('_test_2.py')
    # tempdir = tempfile.TemporaryDirectory(dir=op.expanduser('~/tmp'))
    # lrp = tempdir.name
    lrp = tempfile.mkdtemp(dir=op.expanduser('~/tmp'))
    # Submit jobs
    js = ascommon.cluster_tools.JobSubmitter(
        job_name, connection_string, lrp,
        walltime='01:00:00',
        concurrent_job_limit=concurrent_job_limit,
        queue='short',
        priority="0")
    logger.info('Submitting...')
    with js.connect():
        js.submit(iterable)
    logger.info('Finished submitting...')
    # Make sure that jobs finish successfully
    time_0 = time.time()
    for job_id, system_command in iterable:
        while True:
            stderr_file = op.join(js.log_path, '{}.err'.format(job_id))
            try:
                with open(stderr_file, 'rt') as ifh:
                    stderr = ifh.read().strip()
            except FileNotFoundError:
                logger.debug("File '{}' not found!".format(stderr_file))
                time.sleep(5)
                continue
            if stderr.endswith('DONE!'):
                logger.error(stderr)
                assert False, (job_id, system_command)
            elif stderr.endswith('ERROR!'):
                break
            else:
                logger.debug('Waiting for job to finish...')
                time.sleep(5)
                if (time.time() - time_0) <= 10 * 60:
                    continue
                else:
                    logger.error(stderr)
                    assert False, "Timeout!"


class TestJobStatus:

    @classmethod
    def setup_class(cls):
        cls.job_name = 'test_logs_1'
        cls.connection_string = 'sge://:@192.168.0.1'
        cls.log_base_dir = op.abspath(op.splitext(__file__)[0])
        cls.log_dir = op.join(cls.log_base_dir, cls.job_name)
        os.makedirs(cls.log_dir, exist_ok=True)
        with tarfile.open(cls.log_dir + '.tar.gz') as t:
            t.extractall(cls.log_dir)

    @classmethod
    def teardown_class(cls):
        shutil.rmtree(cls.log_dir)

    def test_job_status(self):
        js = ascommon.cluster_tools.JobSubmitter(
            self.job_name,
            self.connection_string,
            self.log_base_dir,
            force_new_folder=False)
        results_df = js.job_status([(i, i) for i in range(3360)])
        assert (Counter(results_df['status']) ==
                Counter({'done': 2650, 'running': 387, 'missing': 323}))
