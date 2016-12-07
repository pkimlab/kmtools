import os
import os.path as op
import shutil
import time
import tempfile
import subprocess
import pytest
import pandas as pd

from kmtools import system_tools, sequence_tools

TEMPDIR = op.abspath(op.splitext(__file__)[0])
os.makedirs(TEMPDIR, exist_ok=True)


class TestUniParcXMLParser:

    @classmethod
    def setup_method(self, method):
        self.file_path = op.join(op.splitext(__file__)[0], 'test_uniparc_all.xml.gz')
        self.output_dir = tempfile.mkdtemp(dir=TEMPDIR)

    @classmethod
    def teardown_method(self, method):
        shutil.rmtree(self.output_dir)

    @pytest.mark.parametrize("match,output", [
        ('''type="EMBL" id="AAF63732" version_i="1" active="Y" version="1" '''
         '''created="2003-03-12" last="2015-11-15"''',
         {'type': 'EMBL', 'id': 'AAF63732', 'version_i': '1', 'active': 'Y', 'version': '1',
          'created': '2003-03-12', 'last': '2015-11-15'}),
        ('''type="protein_name" value="aminotransferase family protein\\"''',
         {'type': 'protein_name', 'value': 'aminotransferase family protein'})
    ])
    def test__parse_match(self, match, output):
        self.parser = sequence_tools.UniParcXMLParser(self.file_path, self.output_dir)
        assert self.parser._parse_match(match) == output

    def test__append_to_file(self):
        self.parser = sequence_tools.UniParcXMLParser(
            self.file_path, self.output_dir, 'csv')
        data = [
            {'uniparc_id': 1},
            {'uniparc_id': 2, 'dataset': 'uniparc'},
        ]
        self.parser._append_to_file('uniparc.tsv', data, self.parser._uniparc_columns)

    @pytest.mark.skipif(
        pytest.config.getvalue("quick"), reason="Tests take several minutes.")
    @pytest.mark.parametrize("writer", ['pandas', 'csv'])
    def test_run(self, writer):
        with system_tools.open_compressed(self.file_path, 'rt') as ifh:
            self.parser = sequence_tools.UniParcXMLParser(
                ifh, self.output_dir, 'pandas')
            t0 = time.time()
            self.parser.parse()
            t1 = time.time()
        print("Finished in {:.2f} seconds".format(t1 - t0))
        self._assert_dataframes_match()

    @pytest.mark.skipif(
        shutil.which('pypy') is None, reason="`pypy` must be installed.")
    @pytest.mark.parametrize("optimize", ['', '-O'])
    def test_run_pypy(self, optimize):
        proc = subprocess.run(['which', 'pypy'], stdout=subprocess.PIPE, universal_newlines=True)
        proc.check_returncode()
        if not proc.stdout.strip():
            print("No `pypy` installed, but running with `pypy` is over 15 times faster!")
            return

        system_command = (
            "zcat {} | pypy {} '{}' --output_dir '{}'"
            .format(
                self.file_path,
                optimize,
                sequence_tools.uniparc_xml_parser.__file__,
                self.output_dir)
        )
        print(system_command)
        t0 = time.time()
        proc = subprocess.run(
            system_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True, shell=True)
        t1 = time.time()
        print("Finished in {:.2f} seconds".format(t1 - t0))
        proc.check_returncode()
        self._assert_dataframes_match()

    def _assert_dataframes_match(self):
        filenames = [
            'uniparc.tsv', 'uniparc_sequence.tsv', 'uniparc_xref.tsv', 'uniparc_xref_prop.tsv'
        ]
        for filename in filenames:
            df_ref = pd.read_csv(
                op.join(op.splitext(__file__)[0], 'test_uniparc_all_output', filename + '.gz'),
                sep='\t')
            df_test = pd.read_csv(
                op.join(self.output_dir, filename),
                sep='\t')
            assert (df_ref.fillna(0) == df_test.fillna(0)).all().all()
