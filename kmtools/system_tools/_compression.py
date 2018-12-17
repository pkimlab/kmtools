import bz2
import gzip
import logging
import lzma
import os.path as op
import shlex
import subprocess
from contextlib import contextmanager

logger = logging.getLogger(__name__)


@contextmanager
def decompress(filename):
    """Temporarly decompress a file."""
    try:
        print("Gunzipping file '{}'...".format(filename))
        subprocess.check_call(shlex.split("gunzip '{}'".format(filename)))
    except Exception as e:
        print("Unzipping the file failed with an error: {}".format(e))
        raise e
    else:
        yield op.splitext(filename)[0]
    finally:
        print("Gzipping the file back again...")
        subprocess.check_call(shlex.split("gzip '{}'".format(filename.rstrip(".gz"))))


@contextmanager
def open_compressed(filename, mode="rb"):
    """Open a potentially compressed file.

    Note
    ----
    It's probably best to use `smart_open` for .gz files.
    """
    if filename.endswith(".gz"):
        fh = gzip.open(filename, mode)
    elif filename.endswith(".bz2"):
        fh = bz2.open(filename, mode)
    elif filename.endswith(".xz"):
        fh = lzma.open(filename, mode)
    else:
        fh = open(filename, mode)
    try:
        yield fh
    finally:
        # Close the filehandle even if an error occurs
        fh.close()
