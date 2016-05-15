import os
import tempfile
import pycurl
from retrying import retry


def set_temp_dir(tempdir=None):
    """Set the temporary directory to be used by `tempfile`."""
    if tempdir is None:
        if 'TMPDIR' in os.environ:
            tempdir = os.environ['TMPDIR']
        else:
            # use tempfile defaults
            return
    os.makedirs(tempdir, exist_ok=True)
    tempfile.tempdir = tempdir


@retry(
    retry_on_exception=lambda exc: isinstance(exc, pycurl.error),
    wait_exponential_multiplier=1000,  # milliseconds
    wait_exponential_max=60000,  # milliseconds
    stop_max_attempt_number=7)
def download(url, output_file):
    """Download file from 'url' into 'output_file'."""
    with open(output_file, 'wb') as ofh:
        c = pycurl.Curl()
        c.setopt(c.URL, url)
        c.setopt(c.WRITEDATA, ofh)
        c.perform()
        c.close()
