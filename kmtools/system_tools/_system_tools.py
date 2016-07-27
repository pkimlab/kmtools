import os
import os.path as op
import logging
import contextlib
import fcntl
import pycurl
import paramiko
from retrying import retry
from contextlib import contextmanager
import subprocess

logger = logging.getLogger(__name__)


@contextmanager
def open_exclusively(filename, mode='a'):
    fd = os.open(filename, os.O_CREAT | os.O_RDWR)
    fcntl.lockf(fd, fcntl.LOCK_EX)
    try:
        f = os.fdopen(fd, mode)
        yield f
    except:
        raise
    finally:
        f.close()


@contextmanager
def decompress(file):
    """Temporarly decompress a file."""
    try:
        print("Gunzipping file '{}'...".format(file))
        subprocess.check_call("gunzip '{}'".format(file), shell=True)
    except Exception as e:
        print('Unzipping the file failed with an error: {}'.format(e))
        raise e
    else:
        yield
    finally:
        print("Gzipping the file back again...")
        subprocess.check_call("gzip '{}'".format(file.rstrip('.gz')), shell=True)


@contextmanager
def switch_paths(working_path):
    current_path = os.getcwd()
    try:
        os.chdir(working_path)
        yield
    except:
        raise
    finally:
        os.chdir(current_path)


@retry(
    retry_on_exception=lambda exc: isinstance(exc, pycurl.error),
    wait_exponential_multiplier=1000,  # milliseconds
    wait_exponential_max=60000,  # milliseconds
    stop_max_attempt_number=7)
def download(url, output_file):
    """Download file from 'url' into 'output_file'.

    .. deprecated::
       Use `import urllib.request` instead.
    """
    with open(output_file, 'wb') as ofh:
        c = pycurl.Curl()
        c.setopt(c.URL, url)
        c.setopt(c.WRITEDATA, ofh)
        c.perform()
        c.close()


@contextlib.contextmanager
def ssh_client(ip):
    try:
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(ip)
        yield ssh
    except Exception as e:
        logger.error(type(e))
        logger.error(e)


def make_tarfile(source_dir, output_filename):
    """Compress folder into a `*.tar.gz` file."""
    import tarfile
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=op.basename(source_dir))
