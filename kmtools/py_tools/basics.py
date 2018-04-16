import inspect
import logging

logger = logging.getLogger(__name__)


class irange:
    """Infinite range containing `[start, stop)`.

    Examples:
        >>> 1.5 in irange(1, 2)
        True
        >>> 11.1 in irange(11.1, 12)
        True
        >>> -10.0 in irange(-20, -10.0)
        False
    """
    start: float
    stop: float

    def __init__(self, start: float, stop: float) -> None:
        self.start = start
        self.stop = stop

    def __contains__(self, value: float) -> bool:
        if self.start is not None and value < self.start:
            return False
        if self.stop is not None and value >= self.stop:
            return False
        return True


def uniquify(l):
    """Return a list of unique elements of `l`, preserving order.

    Examples:
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
