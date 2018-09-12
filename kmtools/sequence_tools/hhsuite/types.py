from pathlib import Path

import attr


@attr.s
class HHInput:
    #: Path to the fasta file with the input sequence.
    sequence_file: Path = attr.ib()
    #: Directory for storing temporary files.
    temp_dir: Path = attr.ib()


@attr.s(kw_only=True)
class HHBlits(HHInput):
    a3m_file: Path = attr.ib()
    hhm_file: Path = attr.ib()
    hhblits_hhr_file: Path = attr.ib()
    hhblits_database_dir: Path = attr.ib()
    hhblits_extra_args: str = attr.ib(default="")


@attr.s(kw_only=True)
class HHFilter(HHBlits):
    hhfilter_extra_args: str = attr.ib(default="")


@attr.s(kw_only=True)
class AddSS(HHFilter):
    addss_extra_args: str = attr.ib(default="")


@attr.s(kw_only=True)
class HHMake(AddSS):
    hhmake_extra_args: str = attr.ib(default="")


@attr.s(kw_only=True)
class HHSearch(HHMake):
    hhsearch_hhr_file: Path = attr.ib()
    hhsearch_tab_file: Path = attr.ib()
    hhsearch_database_dir: Path = attr.ib()
    hhsearch_extra_args: str = attr.ib(default="")


@attr.s(kw_only=True)
class HHMakeModel(HHSearch):
    hhmakemodel_pir_file: Path = attr.ib()
    hhmakemodel_extra_args: str = attr.ib(default="")
