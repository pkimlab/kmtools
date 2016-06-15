import logging
import contextlib
import pycurl
import paramiko
from retrying import retry

logger = logging.getLogger(__name__)


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
