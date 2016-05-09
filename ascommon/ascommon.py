def unique(l):
    """Return a list of unique elements of `l`, preserving order."""
    seen = set()
    return [x for x in l if x not in seen or seen.add(x)]
