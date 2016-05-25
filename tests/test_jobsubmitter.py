""".

.. note::
    Tests in this file require access to compute clusters on hard-coded IPs.
"""
import sys
import os.path as op
import logging
import tempfile
import time
import json
import paramiko
import pytest
import ascommon

ascommon.configure_logging('debug')
logger = logging.getLogger(__name__)
logger.info('Hello world')


def _test_ssh_connection(ip):
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    try:
        ssh.connect(ip)
        return True
    except Exception as e:
        logger.error("{}: '{}'".format(type(e), str(e)))
        return False


test_input = [
    (head_node_info, concurrent_job_limit)
    for head_node_info in [('sge', '192.168.6.201'), ('pbs', '192.168.233.150')]
    for concurrent_job_limit in [None, 50]
    if _test_ssh_connection(head_node_info[1])
]
logger.info("Collected {} test inputs".format(len(test_input)))


def get_iterable(script_filename):
    return [
        (str(job_id), '{python} "{script}" -i {job_id}'.format(
            python=sys.executable,
            script=op.join(op.splitext(__file__)[0], script_filename),
            job_id=job_id))
        for job_id in range(99)
    ]


@pytest.mark.skipif(pytest.config.getvalue("quick"), reason="Tests take several minutes.")
@pytest.mark.parametrize("head_node_info, concurrent_job_limit", test_input)
def test_1(head_node_info, concurrent_job_limit):
    """Test on tasks that finish successfully."""
    job_name = 'test_1'
    iterable = get_iterable('test_1.py')
    head_node_type, head_node_ip = head_node_info
    # tempdir = tempfile.TemporaryDirectory(dir=op.expanduser('~/tmp'))
    # lrp = tempdir.name
    lrp = tempfile.mkdtemp(dir=op.expanduser('~/tmp'))
    # Submit jobs
    with ascommon.cluster_tools.JobSubmitter(
            job_name, head_node_type, head_node_ip, lrp,
            walltime='01:00:00',
            queue='short',
            concurrent_job_limit=concurrent_job_limit) as js:
        logger.info('Submitting...')
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
@pytest.mark.parametrize("head_node_info, concurrent_job_limit", test_input)
def test_2(head_node_info, concurrent_job_limit):
    """Test on tasks that finish crash."""
    job_name = 'test_2'
    iterable = get_iterable('test_2.py')
    head_node_type, head_node_ip = head_node_info
    # tempdir = tempfile.TemporaryDirectory(dir=op.expanduser('~/tmp'))
    # lrp = tempdir.name
    lrp = tempfile.mkdtemp(dir=op.expanduser('~/tmp'))
    # Submit jobs
    with ascommon.cluster_tools.JobSubmitter(
            job_name, head_node_type, head_node_ip, lrp,
            walltime='01:00:00',
            queue='short',
            concurrent_job_limit=concurrent_job_limit) as js:
        logger.info('Submitting...')
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
