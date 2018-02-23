from typing import Any, Callable, Generator


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
