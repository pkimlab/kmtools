import codecs
import functools
import importlib
import inspect
import logging
import os
import pkgutil
import signal
import time
from contextlib import contextmanager

import dill
import psutil

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
    """Print a logging message whenever the decorated function or module is called."""
    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        logger.warning(
            fn.__name__ + '(' + ', '.join(args) +
            ', '.join('{}={}'.format(k, v) for k, v in kwargs.items()) + ')')
        return fn(*args, **kwargs)
    return wrapper


def serialize_object_to_string(obj):
    """Serialize a Python object into a string.

    Useful for when you want to pass an object as a command-line
    argument to a file.

    Parameters
    ----------
    obj : object
        The object that you want to serialize.

    Returns
    -------
    data_hex_string : str
        String representation of the object.
    """
    data_bytes = dill.dumps(obj)
    data_hex = codecs.encode(data_bytes, "hex")
    data_hex_string = data_hex.decode('utf-8')
    return data_hex_string


def deserialize_object_from_string(data_hex_string):
    """Deserialize a Python object from a string.

    Parameters
    ----------
    data_hex_string : str
        String representation of the Python object.

    Returns
    -------
    obj : object
        The object that was encoded inside the string argument.
    """
    data_hex = data_hex_string.encode('utf-8')
    data_bytes = codecs.decode(data_hex, "hex")
    obj = dill.loads(data_bytes)
    return obj


def uniquify(l):
    """Return a list of unique elements of `l`, preserving order.

    Examples
    --------
    >>> list(uniquify([1, 2, 3, 2, 1]))
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


def kill_child_processes(parent_pid, sig=signal.SIGTERM):
    try:
        parent = psutil.Process(parent_pid)
    except psutil.NoSuchProcess:
        logger.debug("Counld not find parent process with pid '%s'", parent_pid)
        return
    children = parent.children(recursive=True)
    for process in children:
        logger.debug("Killing process '%s'...", process)
        process.send_signal(sig)
