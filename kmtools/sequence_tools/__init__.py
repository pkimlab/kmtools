"""

Sequences
---------

.. autosummary::
    :toctree: _modules

    fetch_sequence
    mutation_matches_sequence
    mutation_in_domain
    format_hgvs_mutation

Alignments
----------

.. autosummary::
    :toctree: _modules

    align_pairwise
    get_crossmapping
    find_in_set
    align_sw

HH-Suite
--------

.. autosummary::
    :toctree: _modules

    hhblits
    hhfilter
    addss
    hhmake
    hhsearch
    hhmakemodel

    parse_hhr_data
"""
try:
    from .hashes_fast import crc64
except ImportError:
    import warnings

    warnings.warn("Could not import Cythonized functions!")
    from .hashes_slow import crc64

from .alignment_tools import *
from .blast import *
from .sequence_tools import *
from .codons_info import *
from .dna_util import *
from .slims import *
from .fasta_tools import *
from .fastq_tools import *
from .hhsuite import *
from . import *
