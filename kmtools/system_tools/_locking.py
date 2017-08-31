import fcntl
import logging
import os
from contextlib import contextmanager

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
