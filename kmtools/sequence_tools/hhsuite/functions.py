import logging
import shlex
import subprocess
from contextlib import closing
from pathlib import Path
from typing import Optional
from kmtools import py_tools

from .types import HHBlits, HHInput, HHFilter, HHAddSS, HHMake, HHSearch, HHMakeModel

logger = logging.getLogger(__name__)


def run(cmd: str, cwd: Path) -> None:
    logger.debug("system_commmand: '%s'", cmd)
    with closing(py_tools.LogPipe(logger.debug)) as pipe:
        subprocess.run(  # type: ignore
            shlex.split(cmd), stdout=pipe, stderr=pipe, universal_newlines=True, cwd=cwd, check=True
        )


def hhblits(input_: HHInput, database: Path, preset: Optional[str] = None) -> HHBlits:
    data = HHBlits(
        hhblits_a3m_file=input_.sequence_file.with_suffix(".a3m"),
        hhblits_hhm_file=input_.sequence_file.with_suffix(".hhm"),
        hhblits_hhr_file=input_.sequence_file.with_suffix(".hhr"),
        hhblits_database=database,
        hhblits_preset=preset,
        **vars(input_),
    )

    cmd = (
        f"hhblits"
        f" -i {data.sequence_file}"
        f" -o {data.hhblits_hhr_file}"
        f" -ohhm {data.hhblits_hhm_file}"
        f" -oa3m {data.hhblits_a3m_file}"
        f" -d {data.hhblits_database}{data.hhblits_database.name}"
    )
    if data.hhblits_preset is None:
        cmd += (
            #
            f" -oalis hhblits_alis"
            f" -n 3"
            f" -mact 0.5"
            f" -cpu 1"
        )
    elif data.hhblits_preset == "hhpred":
        cmd += (
            #
            f" -all"
            f" -maxfilt 100000"
            f" -realign_max 100000"
            f" -B 100000"
            f" -Z 100000"
        )
    run(cmd, data.temp_dir)

    return data


def hhfilter(input_: HHBlits) -> HHFilter:
    data = HHFilter(**vars(input_))
    data.hhblits_a3m_file = input_.hhblits_a3m_file.with_suffix(".filt.a3m")
    cmd = (
        f"hhfilter"
        f" -i '{input_.hhblits_a3m_file}'"
        f" -o '{data.hhblits_a3m_file}'"
        f" -id 90"
        f" -neff 15"
        f" -qsc -30"
    )
    run(cmd, data.temp_dir)
    return data


def addss(input_: HHBlits) -> HHAddSS:
    data = HHAddSS(**vars(input_))
    data.hhblits_a3m_file = input_.hhblits_a3m_file.with_suffix(".withss.a3m")
    cmd = (
        #
        f"addss.pl"
        f" '{input_.hhblits_a3m_file}'"
        f" '{data.hhblits_a3m_file}'"
        f" -a3m"
    )
    run(cmd, data.temp_dir)
    return data


def hhmake(input_: HHBlits) -> HHMake:
    data = HHMake(
        hhmake_hhm_file=input_.hhblits_hhm_file.with_suffix(".hhmake.hhm"), **vars(input_)
    )
    cmd = (
        #
        f"hhmake"
        f" -i '{data.hhblits_hhm_file}'"
        f" -o '{data.hhmake_hhm_file}'"
    )
    run(cmd, data.temp_dir)
    return data


def hhsearch(input_: HHMake, database: Path) -> HHSearch:
    data = HHSearch(
        hhsearch_tab_file=..., hhsearch_hhr_file=..., hhsearch_database=database, **vars(input_)
    )
    cmd = (
        f"hhsearch"
        f" -i {data.hhmake_hhm_file}"
        f" -o {data.hhsearch_hhr_file}"
        f" -d {data.hhsearch_database}/{data.hhsearch_database.name}"
        f" -mact 0.05"
        f" -cpu 1"
        f" -atab {data.hhsearch_tab_file}"
        # AS additions
        f" -e 0.1"
        f" -glob"
    )
    run(cmd, data.temp_dir)
    return data


def hhmakemodel(input_: HHSearch) -> HHMakeModel:
    data = HHMakeModel(hhmakemodel_pir_file=..., **vars(input_))
    cmd = (
        f"hhmakemodel.pl"
        f" {data.hhsearch_hhr_file}"
        f" -m 1"
        f" -q {data.hhblits_a3m_file}"
        f" -v 1"
        f" -d {data.hhsearch_database}/{data.hhsearch_database.name}"
        f" -pir {data.hhmakemodel_pir_file}"
    )
    run(cmd, data.temp_dir)
    return data
