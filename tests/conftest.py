import os.path as op
import sys

# Use the installed version of the package instead of the current directory
sys.path.remove(op.dirname(op.abspath(__file__)))


def pytest_addoption(parser):
    parser.addoption("--quick", action="store_true", help="Run only quick tests.")
