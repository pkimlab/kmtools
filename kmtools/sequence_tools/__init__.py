"""Sequence tools

.. autosummary::
    :toctree: _modules

    codons_info
    dna_util
    elms
    run_psipred
    read_psipred
"""
try:
    from ._hashes import crc64
except ImportError:
    import warnings

    warnings.warn("Could not import Cythonized functions!")
    from .hashes import crc64_slow as crc64

    del warnings

from .hashes import crc64_slow
from .alignment_tools import *
from .blast import *
from .sequence_tools import *
from .codons_info import *
from .dna_util import *
from .slims import *
from .psipred import run_psipred, read_psipred
