import fcntl
import logging
import os
import time
from contextlib import contextmanager

logger = logging.getLogger(__name__)


@contextmanager
def open_exclusively(filename, mode="a"):
    fd = os.open(filename, os.O_CREAT | os.O_RDWR)
    fcntl.lockf(fd, fcntl.LOCK_EX)
    try:
        f = os.fdopen(fd, mode)
        yield f
    except Exception:
        raise
    finally:
        f.close()


@contextmanager
def get_lock(name):
    lock = None

    def close_lock(lock):
        if lock is not None:
            lock.close()
            os.remove(lock.name)

    while True:
        try:
            lock = open(name + ".lock", "x")
            yield lock
            close_lock(lock)
            break
        except FileExistsError:
            time.sleep(60)
        except Exception:
            close_lock(lock)
            raise
