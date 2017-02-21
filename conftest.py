import os.path as op
import sys
import pytest

# Use the installed version of the package instead of the current directory
sys.path.remove(op.dirname(op.abspath(__file__)))


def pytest_addoption(parser):
    parser.addoption("--quick", action="store_true", help="Run only quick tests.")


@pytest.fixture(autouse=True)
def set_doctest_namespace(doctest_namespace):
    import os
    import os.path as op
    import tempfile
    import numpy as np
    doctest_namespace['os'] = os
    doctest_namespace['op'] = op
    doctest_namespace['tempfile'] = tempfile
    doctest_namespace['np'] = np
