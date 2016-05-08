import os
import os.path as op

blacklist = [
    'graph_tool',
    'g2d',
    'ipython_tools',
    'openmm_tools',
    'optimized',
    'pymol_tools',
]

__all__ = [
    op.splitext(f)[0]  # remove .py extension
    for f in os.listdir(op.dirname(__file__))  # list contents of current dir
    if op.isfile(f) and f.endswith('.py') and not f.startswith('_') and
    not any(f.startswith(pre) for pre in blacklist)
]

from . import *  # noqa
