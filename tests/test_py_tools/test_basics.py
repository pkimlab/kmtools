import pytest

from kmtools import py_tools


@pytest.mark.parametrize("duplicated, unique", [
    ([1, 2, 3, 2, 1], [1, 2, 3]),
    ([(1), (1, 2), (1)], [(1), (1, 2)]),
])
def test_uniquify(duplicated, unique):
    assert tuple(py_tools.uniquify(duplicated)) == tuple(unique)
