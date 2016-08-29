import os
import time
import logging
import inspect
import signal
from contextlib import contextmanager

logger = logging.getLogger(__name__)


def unique(l):
    """Return a list of unique elements of `l`, preserving order.

    Examples
    --------
    >>> list(unique([1, 2, 3, 2, 1]))
    [1, 2, 3]
    """
    seen = set()
    for x in l:
        if x not in seen:
            seen.add(x)
            yield x


@contextmanager
def print_heartbeats():
    """Spawn a fork that prints a message every minute.

    This is required for travis-ci.
    """
    from elaspic import conf
    # Don't print random stuff if not testing
    if not conf.CONFIGS['testing']:
        yield
        return
    # Print a heartbeat to keep travis happy.
    pid = os.fork()
    if pid == 0:
        while True:
            time.sleep(60)
            logger.info("Subprocess is still running...")
        os._exit()
    try:
        yield
    finally:
        os.kill(pid, 15)
        os.waitpid(pid, 0)


@contextmanager
def get_lock(name):
    lock = None

    def close_lock(lock):
        if lock is not None:
            lock.close()
            os.remove(lock.name)

    while True:
        try:
            lock = open(name + '.lock', 'x')
            yield lock
            close_lock(lock)
            break
        except FileExistsError:
            time.sleep(60)
        except:
            close_lock(lock)
            raise


def decorate_all_methods(decorator):
    """Decorate all methods of a class with `decorator`."""
    def apply_decorator(cls):
        for k, f in cls.__dict__.items():
            if inspect.isfunction(f):
                setattr(cls, k, decorator(f))
        return cls
    return apply_decorator


def kill_child_process(child_process):
    if child_process.poll() is not None:
        print('Child process with pid {} already terminated with return code {}'
              .format(child_process.pid, child_process.returncode))
        return
    try:
        print('Trying to terminate gracefully child process with pid: {}'.child_process.pid)
        os.killpg(child_process.pid, signal.SIGTERM)
#        child_process.terminate()
    except Exception as e:
        print("Didn't work because of error: {}".format(e.__str__()))
        try:
            print('Trying to kill child process...')
            os.killpg(child_process.pid, signal.SIGKILL)
#            child_process.kill()
        except:
            print("Didn't work because of error: {}".format(e.__str__()))
            print("Letting it go...")
            pass
    print('OK')
