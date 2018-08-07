import pytest

import kmtools


@pytest.mark.parametrize('attribute', ['__version__'])
def test_attribute(attribute):
    assert getattr(kmtools, attribute)


def test_main():
    import kmtools

    assert kmtools
