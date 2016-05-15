import os
import os.path as op
import shutil
import time
import tempfile
import subprocess
import shlex
import pytest
import pandas as pd
import ascommon

TEMPDIR = op.abspath(op.splitext(__file__)[0])
os.makedirs(TEMPDIR, exist_ok=True)


class TestHackyXMLParser:

    @classmethod
    def setup_method(self, method):
        self.file_path = op.join(op.splitext(__file__)[0], 'test_uniparc_all.xml.gz')
        self.output_dir = tempfile.mkdtemp(dir=TEMPDIR)

    @classmethod
    def teardown_method(self, method):
        shutil.rmtree(self.output_dir)

    def test__parse_match(self):
        match = """\
type="EMBL" id="AAF63732" version_i="1" active="Y" version="1" created="2003-03-12" \
last="2015-11-15"\
"""
        output = {
            'active': 'Y',
            'created': '2003-03-12',
            'id': 'AAF63732',
            'last': '2015-11-15',
            'type': 'EMBL',
            'version': '1',
            'version_i': '1'
        }
        self.parser = ascommon.sequence_tools.HackyXMLParser(self.file_path, self.output_dir)
        assert self.parser._parse_match(match) == output

    def test__append_to_file(self):
        self.parser = ascommon.sequence_tools.HackyXMLParser(
            self.file_path, self.output_dir, 'csv')
        data = [
            {'uniparc_id': 1},
            {'uniparc_id': 2, 'dataset': 'uniparc'},
        ]
        self.parser._append_to_file('uniparc.tsv', data, self.parser._uniparc_columns)

    @pytest.mark.parametrize("writer", ['pandas', 'csv'])
    def test_run(self, writer):
        self.parser = ascommon.sequence_tools.HackyXMLParser(
            self.file_path, self.output_dir, 'pandas')
        t0 = time.time()
        self.parser.parse()
        t1 = time.time()
        print("Finished in {:.2f} seconds".format(t1 - t0))
        self._assert_dataframes_match()

    @pytest.mark.parametrize("optimize", ['', '-O'])
    def test_run_pypy(self, optimize):
        proc = subprocess.run(['which', 'pypy'], stdout=subprocess.PIPE, universal_newlines=True)
        proc.check_returncode()
        if not proc.stdout.strip():
            print("No `pypy` installed, but running with `pypy` is over 15 times faster!")
            return
        import ascommon.sequence_tools.hacky_xml_parser
        system_command = (
            "pypy {} '{}' --file_path '{}' --output_dir '{}'"
            .format(
                optimize,
                ascommon.sequence_tools.hacky_xml_parser.__file__,
                self.file_path,
                self.output_dir)
        )
        print(system_command)
        t0 = time.time()
        proc = subprocess.run(
            shlex.split(system_command), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True)
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

if __name__ == '__main__':
    import pytest
    pytest.main([__file__, '-sv'])
