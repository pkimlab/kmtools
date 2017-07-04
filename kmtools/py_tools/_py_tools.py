import functools
import importlib
import inspect
import logging
import os
import signal
import time
import pkgutil
from contextlib import contextmanager

logger = logging.getLogger(__name__)


def iter_submodules(package):
    """Import all submodules of a module, recursively, including subpackages.

    Adapted from https://stackoverflow.com/a/25562415/2063031
    """
    yield package.__name__, package
    for loader, name, ispkg in pkgutil.walk_packages(package.__path__):
        module = importlib.import_module(package.__name__ + '.' + name)
        if ispkg:
            yield from iter_submodules(module)
        else:
            yield module.__name__, module


def log_calls(fn):
    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        logger.warning(
            fn.__name__ + '(' + ', '.join(args) +
            ', '.join('{}={}'.format(k, v) for k, v in kwargs.items()) + ')')
        return fn(*args, **kwargs)
    return wrapper


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


def strip_ps(name, prefix=None, suffix=None):
    """Remove `prefix` and / or `suffix` from `name`.

    Examples
    --------
    >>> strip_ps('good_god_gomer', 'good', 'gomer')
    '_god_'
    """
    if prefix and name.startswith(prefix):
        name = name[len(prefix):]
    if suffix and name.endswith(suffix):
        name = name[:-len(suffix)]
    return name


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
