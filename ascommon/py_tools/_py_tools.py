def unique(l):
    """Return a list of unique elements of `l`, preserving order.

    Examples
    --------
    >>> list(unique([1, 2, 3, 2, 1]))
    [1, 2, 3]
    """
    seen = set()
    yield from (x for x in l if x not in seen and not seen.add(x))


def format_unprintable(string):
    r"""Escape tabs (\t), newlines (\n), etc. for system commands and printing.

    Examples
    --------
    >>> format_unprintable('\t')
    '\\t'
    """
    return repr(string).strip("'")
