import importlib
import logging
import pkgutil
from typing import Any, Callable, Generator

logger = logging.getLogger(__name__)


def iter_forever(iterable: Callable[[], Generator[Any, Any, Any]]) -> Generator[Any, Any, Any]:
    """Iterate over an iterable forever.

    Like `itertools.cycle`, but without storing the seen elements in memory.

    Examples:
        >>> import itertools
        >>> def foo():
        ...     yield from range(3)
        >>> list(itertools.islice(iter_forever(foo), 7))
        [0, 1, 2, 0, 1, 2, 0]
    """
    while True:
        yield from iterable()


def iter_submodules(package):
    """Import all submodules of a module, recursively, including subpackages.

    Adapted from https://stackoverflow.com/a/25562415/2063031
    """
    yield package.__name__, package
    for loader, name, ispkg in pkgutil.walk_packages(package.__path__):
        try:
            module = importlib.import_module(package.__name__ + "." + name)
        except ModuleNotFoundError as e:
            logger.error("Could not import module '%s' (%s)", module, e)
            continue
        if ispkg:
            yield from iter_submodules(module)
        else:
            yield module.__name__, module
