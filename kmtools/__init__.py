import os
import os.path as op

_blacklist = [
    'graph_tools',
    'omnia_tools',
    'pymol_tools',
    'jupyter_tools',
    'plotting_tools',
]

BASE_DIR = op.dirname(__file__)

__all__ = [
    op.splitext(f)[0]  # remove .py extension
    for f in os.listdir(BASE_DIR)  # list contents of current dir
    if not f.startswith('_') and
    ((op.isfile(op.join(BASE_DIR, f)) and f.endswith('.py')) or
     (op.isdir(op.join(BASE_DIR, f)) and op.isfile(op.join(BASE_DIR, f, '__init__.py')))) and
    not any(f.startswith(pre) for pre in _blacklist)
]

from . import *  # noqa
