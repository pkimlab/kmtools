import doctest
import logging
import os
import os.path as op
import tempfile

import numpy as np
import pandas as pd
import pytest

import kmtools.structure_tools
from kmtools import py_tools

logger = logging.getLogger(__name__)

DOCTEST_OPTIONFLAGS = (
    doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS | doctest.IGNORE_EXCEPTION_DETAIL)

DOCTEST_EXTRAGLOBS = {'os': os, 'op': op, 'tempfile': tempfile, 'np': np, 'pd': pd}


@pytest.mark.parametrize("module_name, module", py_tools.iter_submodules(kmtools.structure_tools))
def test_doctest(module_name, module):
    pd.options.display.width = 80
    pd.options.display.max_rows = 100
    pd.options.display.max_columns = 100
    failure_count, test_count = doctest.testmod(
        module, optionflags=DOCTEST_OPTIONFLAGS, extraglobs=DOCTEST_EXTRAGLOBS)
    assert failure_count == 0
