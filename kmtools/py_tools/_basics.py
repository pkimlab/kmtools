import importlib
import inspect
import logging
import pkgutil

logger = logging.getLogger(__name__)


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


def decorate_all_methods(decorator):
    """Decorate all methods of a class with `decorator`."""
    def apply_decorator(cls):
        for k, f in cls.__dict__.items():
            if inspect.isfunction(f):
                setattr(cls, k, decorator(f))
        return cls
    return apply_decorator


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
