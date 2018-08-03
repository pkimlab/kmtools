from hypothesis import given
from hypothesis.strategies import text

from kmtools.sequence_tools import crc64, crc64_slow


@given(text())
def test_crc64(s):
    assert crc64(s) == crc64_slow(s)
