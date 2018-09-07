from pathlib import Path
from typing import Optional

import attr


@attr.s
class HHInput:
    sequence_file: Path = attr.ib()
    temp_dir: Path = attr.ib()


@attr.s
class HHBlits(HHInput):
    hhblits_a3m_file: Path = attr.ib()
    hhblits_hhm_file: Path = attr.ib()
    hhblits_hhr_file: Path = attr.ib()
    #
    hhblits_database: Path = attr.ib()
    hhblits_preset: str = attr.ib()


@attr.s
class HHFilter(HHBlits):
    pass


@attr.s
class HHAddSS(HHBlits):
    pass


@attr.s
class HHMake(HHBlits):
    hhmake_hhm_file: Path = attr.ib()


@attr.s
class HHSearch(HHMake):
    hhsearch_tab_file: Path = attr.ib()
    hhsearch_hhr_file: Path = attr.ib()
    #
    hhsearch_database: Path = attr.ib()


@attr.s
class HHMakeModel(HHSearch):
    hhmakemodel_pir_file: Path = attr.ib()
