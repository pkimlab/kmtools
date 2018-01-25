import pytest

from kmtools import py_tools


@pytest.mark.parametrize("value, bounds, result", [(0, (0, 1), True), (1, (0, 1), False)])
def test_irange(value, bounds, result):
    assert result == (value in py_tools.irange(*bounds))


@pytest.mark.parametrize("duplicated, unique", [
    ([1, 2, 3, 2, 1], [1, 2, 3]),
    ([(1), (1, 2), (1)], [(1), (1, 2)]),
])
def test_uniquify(duplicated, unique):
    assert tuple(py_tools.uniquify(duplicated)) == tuple(unique)
