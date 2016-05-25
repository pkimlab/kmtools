def iterate_params(param_grid, _params=None):
    """

    Parameters
    ----------
    parameter_grid : dict
        Keys are parameters. Values are lists of values that the parameters can take.

    Examples
    --------
    >>> from collections import OrderedDict
    >>> list(iter_params(OrderedDict([('a', [1, 2]), ('b': [3, 4])])
    [{'a': 1, 'b': 3}, {'a': 2, 'b': 3}, {'a': 1, 'b': 4}, {'a': 2, 'b': 4}]
    """
    param_grid = param_grid.copy()
    # Don't modify dictionaries in-place
    if _params is None:
        _params = dict()
    # Terminal case
    if not param_grid:
        yield _params
        return
    # Recurse
    key, values = param_grid.popitem()
    for value in values:
        _params[key] = value
        try:
            yield from iter_params(param_grid.copy(), _params.copy())
        except StopIteration:
            continue

assert list(iter_params({'a': [1,2], 'b': [3,4]})) == list(iter_params({'a': [1,2], 'b': [3,4]}))
assert len(list(iter_params(param_grid))) == 2 * 3 * 4 * 4 * 4
