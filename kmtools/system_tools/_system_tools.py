import bz2
import contextlib
import fcntl
import functools
import gzip
import logging
import lzma
import os
import os.path as op
import shlex
import shutil
import string
import subprocess
import time
import urllib.request
from contextlib import contextmanager

import paramiko
import sqlalchemy as sa
from retrying import retry

from kmtools import system_tools

from . import exc

logger = logging.getLogger(__name__)


def slugify(filename_string):
    valid_chars = "-_.()" + string.ascii_letters + string.digits
    return ''.join(c if c in valid_chars else '_' for c in filename_string)


def remove_extensions(filename, extensions):
    """Remove extensions from file.

    Examples
    --------
    >>> remove_extensions('/tmp/a/b/c.d.e.f.g', ['.e', '.f', '.g'])
    '/tmp/a/b/c.d'
    >>> remove_extensions('do.re.mi', ['do', 'mi'])
    'do.re'
    """
    # Add missing '.'
    extensions = [(ext if ext.startswith('.') else '.' + ext) for ext in extensions]
    # Strip extensions
    while True:
        file_name, file_ext = op.splitext(filename)
        if file_ext in extensions:
            filename = file_name
        else:
            break
    return filename


def format_unprintable(string):
    r"""Escape tabs (\t), newlines (\n), etc. for system commands and printing.

    Examples
    --------
    >>> format_unprintable('\t')
    '\\t'
    """
    return repr(string).strip("'")


# Subprocess
def _set_process_group(parent_process_group_id):
    """Set group_id of the child process to the group_id of the parent process.

    This way when you delete the parent process you also delete all the children.
    """
    child_process_id = os.getpid()
    os.setpgid(child_process_id, parent_process_group_id)


@functools.wraps(subprocess.run)
def run(system_command, **vargs):
    if not isinstance(system_command, (list, tuple)):
        system_command = shlex.split(system_command)
    p = subprocess.run(
        system_command, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        preexec_fn=lambda: _set_process_group(os.getpgrp()),
        **vargs)
    p.stdout = p.stdout.strip()
    p.stderr = p.stderr.strip()
    return p


def get_hostname():
    return run('hostname | cut -d. -f1').stdout


def which(bin_name):
    return run('which ' + bin_name).stdout


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
def switch_paths(working_path):
    current_path = os.getcwd()
    try:
        os.chdir(working_path)
        yield
    except:
        raise
    finally:
        os.chdir(current_path)


# Retry
def check_exception(exc, valid_exc):
    logger.error('The following exception occured:\n{}'.format(exc))
    to_retry = isinstance(exc, valid_exc)
    if to_retry:
        logger.error('Retrying...')
    return to_retry


def retry_database(fn):
    """Decorator to keep probing the database untill you succeed."""
    _check_exception = functools.partial(check_exception, valid_exc=sa.exc.OperationalError)
    r = retry(
        retry_on_exception=_check_exception,
        wait_exponential_multiplier=1000,
        wait_exponential_max=60000,
        stop_max_attempt_number=7)
    return r(fn)


def retry_subprocess(fn):
    _check_exception = functools.partial(check_exception, valid_exc=subprocess.SubprocessError)
    r = retry(
        retry_on_exception=_check_exception,
        wait_exponential_multiplier=1000,
        wait_exponential_max=60000,
        stop_max_attempt_number=7)
    return r(fn)


def retry_archive(fn):
    """Decorator to keep probing the database untill you succeed."""
    _check_exception = functools.partial(check_exception, valid_exc=system_tools.exc.ArchiveError)
    r = retry(
        retry_on_exception=_check_exception,
        wait_fixed=2000,
        stop_max_attempt_number=2)
    return r(fn)


@contextmanager
def decompress(filename):
    """Temporarly decompress a file."""
    try:
        print("Gunzipping file '{}'...".format(filename))
        subprocess.check_call("gunzip '{}'".format(filename), shell=True)
    except Exception as e:
        print('Unzipping the file failed with an error: {}'.format(e))
        raise e
    else:
        yield
    finally:
        print("Gzipping the file back again...")
        subprocess.check_call("gzip '{}'".format(filename.rstrip('.gz')), shell=True)


@contextmanager
def open_compressed(filename, mode='rb'):
    """Open a potentially compressed file."""
    if filename.endswith('.gz'):
        fh = gzip.open(filename, mode)
    elif filename.endswith('.bz2'):
        fh = bz2.open(filename, mode)
    elif filename.endswith('.xz'):
        fh = lzma.open(filename, mode)
    else:
        fh = open(filename, mode)
    try:
        yield fh
    finally:
        # Close the filehandle even if an error occurs
        fh.close()


def read_url(url):
    # Read PDB data
    if url.startswith(('ftp:', 'http:', 'https:', )):
        with urllib.request.urlopen(url) as ifh:
            data = ifh.read()
    else:
        with open(url, 'rb') as ifh:
            data = ifh.read()

    # Uncompress
    if url.endswith('.gz'):
        data = gzip.decompress(data)

    return data


# chmod
def copyfile(infile, outfile, mode=None):
    shutil.copyfile(infile, outfile)
    if mode is not None:
        os.chmod(outfile, mode)


def makedirs(path, mode=None, exist_ok=True):
    if mode is None:
        os.makedirs(path, exist_ok=exist_ok)
    else:
        original_umask = os.umask(0)
        try:
            os.makedirs(path, mode=mode, exist_ok=exist_ok)
        finally:
            os.umask(original_umask)


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


def start_subprocess(system_command):
    p = subprocess.Popen(
        shlex.split(system_command),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        bufsize=1)
    return p


def iter_stdout(p):
    for line in p.stdout:
        line = line.strip()
        if ' [Note] ' in line:
            line = line.partition(' [Note] ')[-1]
        if not line:
            if p.poll() is None:
                time.sleep(0.1)
                continue
            else:
                # logger.debug("DONE! (reached an empty line)")
                return
        yield line


# === Run command ====
class MySSHClient:

    def __init__(self, ssh_host):
        self.ssh_host = ssh_host
        self.ssh = paramiko.SSHClient()
        self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    def __enter__(self):
        logger.debug("Initializing SSH client: '{}'".format(self.ssh_host))
        self.ssh.connect(self.ssh_host)
        return self.ssh

    def __exit__(self, exc_type, exc_value, exc_tb):
        self.ssh.close()
        if exc_type or exc_value or exc_tb:
            import traceback
            logger.error(exc_type)
            logger.error(exc_value)
            logger.error(traceback.print_tb(exc_tb))
            return True
        else:
            return False


@retry_subprocess
def run_command(system_command, host=None, *, shell=False, allowed_returncodes=[0]):
    """Run system command either locally or over ssh."""
    # === Run command ===
    logger.debug(system_command)
    if host is not None:
        logger.debug("Running on host: '{}'".format(host))
        with MySSHClient(host) as ssh:
            _stdin, _stdout, _stderr = ssh.exec_command(system_command)
            stdout = _stdout.read().decode()
            stderr = _stderr.read().decode()
            returncode = _stdout.channel.recv_exit_status()
    else:
        logger.debug("Running locally")
        sp = subprocess.run(
            system_command if shell else shlex.split(system_command),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=shell, universal_newlines=True,
        )
        stdout = sp.stdout
        stderr = sp.stderr
        returncode = sp.returncode
    # === Process results ===
    stdout_lower = stdout.lower()
    if returncode not in allowed_returncodes:
        error_message = (
            "Encountered an error: '{}'\n".format(stderr) +
            "System command: '{}'\n".format(system_command) +
            "Output: '{}'\n".format(stdout) +
            "Return code: {}".format(returncode)
        )
        logger.error(error_message)
        raise exc.SubprocessError(
            command=system_command,
            host=host,
            stdout=stdout,
            stderr=stderr,
            returncode=returncode,
        )
    elif 'warning' in stdout_lower or 'error' in stdout_lower:
        logger.warning("Command ran with warnings / errors:\n{}".format(stdout.strip()))
    else:
        logger.debug("Command ran successfully:\n{}".format(stdout.strip()))
    return stdout, stderr, returncode
