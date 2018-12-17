import logging

import pytest

import kmtools

logger = logging.getLogger(__name__)


@pytest.mark.parametrize("s1, s2", [("ABCDF", "A-CEF"), ("AAAGGVVVAAA", "AAAGG---AAA")])
def test_1(s1, s2):
    a2b, b2a = kmtools.sequence_tools.get_crossmapping(s1, s2, skip_mismatch=True)
    logger.debug("(a2b, b2a): ('{}', '{}')".format(a2b, b2a))
    s1_test = _map_a2b(s1, b2a)
    assert s1_test == "".join(a for (a, b) in zip(s1, s2) if a == b)
    s2_test = _map_a2b(s2, a2b)
    assert s2_test == "".join(b for (a, b) in zip(s1, s2) if a == b)


@pytest.mark.parametrize("s1, s2", [("ABCDF", "A-CEF"), ("AAAGGVVVAAA", "AAAGG---AAA")])
def test_2(s1, s2):
    a2b, b2a = kmtools.sequence_tools.get_crossmapping(s1, s2, skip_mismatch=False)
    s1_test = _map_a2b(s1, b2a)
    assert s1_test == "".join(a for a, b in zip(s1, s2) if a != "-" and b != "-")
    s2_test = _map_a2b(s2, a2b)
    assert s2_test == "".join(b for a, b in zip(s1, s2) if a != "-" and b != "-")


def _map_a2b(s, a2b):
    return "".join(
        s.replace("-", "")[kmtools.sequence_tools.find_in_set(i + 1, a2b) - 1]
        for i in range(len(s))
        if kmtools.sequence_tools.find_in_set(i + 1, a2b)
    )
