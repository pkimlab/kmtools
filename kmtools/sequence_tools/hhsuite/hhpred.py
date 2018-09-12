import logging
import os
import shlex
import subprocess
from contextlib import closing
from pathlib import Path

import attr
from elapsam.extern import SequenceTool, ToolCache

from kmtools import py_tools

logger = logging.getLogger(__name__)


@attr.s
class HHPredCache(ToolCache):
    # hhblits
    hhblits_a3m_file: Path = attr.ib(default=Path("hhblits.a3m"))
    hhblits_hhm_file: Path = attr.ib(default=Path("hhblits.hhm"))
    hhblits_hhr_file: Path = attr.ib(default=Path("hhblits.hhr"))
    # addss
    addss_a3m_file: Path = attr.ib(default=Path("addss.a3m"))
    # hhmake
    hhmake_hhm_file: Path = attr.ib(default=Path("hhmake.hhm"))
    # hhsearch
    hhsearch_tab_file = attr.ib(default=Path("hhsearch.tab"))
    hhsearch_hhr_file = attr.ib(default=Path("hhsearch.hhr"))
    # hhmakemodel
    # hhmakemodel_pir_file: Path = attr.ib(default=Path('hhmakemodel.pir'))


class HHPred(SequenceTool):

    cache: HHPredCache

    def _build(self):
        """

        Notes:
            - ``hhsearch`` and ``hhmake`` produce ``{name}.hhr`` files,
              where ``{name}`` corresponds to the *input* filename!
        """
        # hhblits
        proc = self._run_system_command(self.hhblits_cmd)
        # addss
        proc = self._run_system_command(self.addss_cmd)
        # hhmake
        proc = self._run_system_command(self.hhmake_cmd)
        # hhsearch
        proc = self._run_system_command(self.hhsearch_cmd)
        # ---
        # shutil.copyfile(self.cache.hhsearch_hhr_file,
        #                 self.cache.hhsearch_hhr_file.with_name(
        #                     self.cache.hhsearch_hhr_file.name + '-start.hhr'))
        # proc = self._run_system_command(self.hhmakemodel_cmd)
        return proc

    def _run_system_command(self, system_command: str) -> subprocess.CompletedProcess:
        with closing(py_tools.LogPipe(logger.info)) as log_pipe:
            logger.info(f"system_command: '{system_command}'")
            proc = subprocess.run(
                shlex.split(system_command),
                stdout=subprocess.PIPE,
                stderr=log_pipe,
                universal_newlines=True,
                cwd=self.tempdir,
                env={"PATH": os.environ["PATH"], "HHLIB": os.environ["HHLIB"]},
            )
        logger.info(proc.stdout.strip())
        proc.check_returncode()
        return proc

    @property
    def hhblits_cmd(self, add_ss=True) -> str:
        cmd = "hhblits"
        cmd += f" -i {self.sequence_file}"
        cmd += f" -o {self.cache.hhblits_hhr_file}"
        cmd += f" -oalis hhblits_alis"
        cmd += f" -ohhm {self.cache.hhblits_hhm_file}"
        cmd += f" -n 3"
        cmd += f" -mact 0.5"
        cmd += f" -oa3m {self.cache.hhblits_a3m_file}"
        cmd += f" -d {os.environ['UNIPROT20_DB']}"
        cmd += f" -cpu 1"
        return cmd

    @property
    def addss_cmd(self) -> str:
        cmd = "addss.pl"
        cmd += f" {self.cache.hhblits_a3m_file}"
        cmd += f" {self.cache.addss_a3m_file}"
        cmd += f" -a3m"
        return cmd

    @property
    def hhmake_cmd(self) -> str:
        cmd = "hhmake"
        cmd += f" -i {self.cache.addss_a3m_file}"
        cmd += f" -o {self.cache.hhmake_hhm_file}"
        return cmd

    @property
    def hhsearch_cmd(self) -> str:
        cmd = "hhsearch"
        cmd += f" -i {self.cache.hhmake_hhm_file}"
        cmd += f" -o {self.cache.hhsearch_hhr_file}"
        cmd += f" -d {os.environ['PDB70_DB']}"
        cmd += f" -mact 0.05"
        cmd += f" -cpu 1"
        cmd += f" -atab {self.cache.hhsearch_tab_file}"
        # AS additions
        cmd += f" -e 0.1"
        cmd += f" -glob"
        return cmd

    @property
    def hhmakemodel_cmd(self) -> str:
        cmd = "hhmakemodel.pl"
        cmd += f" {self.cache.hhsearch_hhr_file}"
        cmd += f" -m 1"
        cmd += f" -q {self.cache.hhblits_a3m_file}"
        cmd += f" -v 1"
        cmd += f" -d {os.environ['PDB70_DB']}"
        cmd += f" -pir {self.cache.hhmakemodel_pir_file}"
        return cmd
