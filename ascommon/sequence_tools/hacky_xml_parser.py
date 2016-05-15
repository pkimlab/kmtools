import sys
import os
import os.path as op
import gzip
import re
import shlex


class HackyXMLParser:
    """
    """
    def __init__(self, file_path, output_dir, writer='csv'):
        self.file_path = file_path
        self.file_handle = self._get_file_handle(file_path)
        #
        if not op.isdir(output_dir):
            raise Exception("`output_dir` must exist!")
        elif os.listdir(output_dir):
            raise Exception("`output_dir` must be empty!")
        self.output_dir = output_dir
        #
        if writer not in ['pandas', 'csv']:
            raise Exception("Only 'csv' and 'pandas' writers are supported!")
        self.writer = writer
        #
        self._uniparc = {'uniparc_id': None}
        self._uniparc_sequence = {'uniparc_id': None}
        self._uniparc_xref = {'uniparc_xref_id': None}
        self._uniparc_xref_prop = {'uniparc_xref_prop_id': None}
        #
        self._uniparc_xref_id = 0
        self._uniparc_xref_prop_id = 0
        #
        self._uniparc_cache = []
        self._uniparc_sequence_cache = []
        self._uniparc_xref_cache = []
        self._uniparc_xref_prop_cache = []
        #
        self._uniparc_columns = [
            'uniparc_id', 'dataset', 'UniProtKB_exclusion'
        ]
        self._uniparc_sequence_columns = [
            'uniparc_id', 'length', 'checksum', 'sequence'
        ]
        self._uniparc_xref_columns = [
            'uniparc_xref_id', 'uniparc_id', 'type', 'id', 'version_i', 'active', 'version',
            'created', 'last'
        ]
        self._uniparc_xref_prop_columns = [
            'uniparc_xref_prop_id', 'uniparc_xref_id', 'type', 'value'
        ]

    def _get_file_handle(self, file_path):
        """Return an iterator over compressed or uncompressed files."""
        extension = op.splitext(file_path)[-1]
        if extension == '.gz':
            return gzip.open(self.file_path, mode='rt')
        elif extension == '.bz2':
            raise NotImplementedError
        else:
            return open(file_path)

    @property
    def _file_iterator(self):
        for line in self.file_handle:
            line = line.strip(' \n')
            yield line
        self.file_handle.close()

    def parse(self):
        """Parse UniParc XML file in a hacky non-validating manner."""
        for line in self._file_iterator:
            match = []
            for parser, fn in self._parsers:
                match = parser.findall(line)
                if match:
                    assert len(match) == 1
                    fn(match[0])
                    break
            if not match:
                print("Did not match the following line: '{}'\n".format(line))
        # Flush the last chunk
        self._flush_cache()

    @property
    def _parsers(self):
        return [
            # === UniParc ===
            (re.compile('<entry (.*)>'), self._parse_uniparc_start),
            (re.compile('<accession>(\w+)</accession>'), self._parse_uniparc_accession),
            (re.compile('<sequence (.*)>'), self._parse_uniparc_sequence),
            (re.compile('</entry>'), self._parse_uniparc_end),

            # === UniParc XRef ===
            (re.compile('<dbReference (.*)>'), self._parse_uniparc_xref),
            (re.compile('</dbReference>'), self._parse_uniparc_xref_end),

            # === UniParc XRef Prop ===
            (re.compile('<property (.*)/>'), self._parse_uniparc_xref_prop),

            # === Junk ===
            (re.compile('<?xml.*|<uniparc.*|</uniparc>'), lambda x: None)
        ]

    def _parse_match(self, match):
        def split_kv(kv):
            k, _, v = kv.partition('=')
            v = v.strip('"')
            return k, v
        data = dict(split_kv(kv) for kv in shlex.split(match))
        return data

    # === UniParc ===
    def _parse_uniparc_start(self, match):
        assert self._uniparc['uniparc_id'] is None
        #
        data = self._parse_match(match)
        assert all(c in self._uniparc_columns for c in data)
        self._uniparc.update(data)

    def _parse_uniparc_accession(self, match):
        assert match.startswith('UPI')
        self._uniparc['uniparc_id'] = match

    def _parse_uniparc_sequence(self, match):
        assert (self._uniparc['uniparc_id'] is not None and
                self._uniparc_sequence['uniparc_id'] is None)
        self._uniparc_sequence['uniparc_id'] = self._uniparc['uniparc_id']
        #
        data = self._parse_match(match)
        assert all(c in self._uniparc_sequence_columns for c in data)
        self._uniparc_sequence.update(data)
        #
        line = ''
        sequence = ''
        for line in self._file_iterator:
            if line != '</sequence>':
                sequence += line.strip()
            else:
                break
        assert sequence.isupper()
        self._uniparc_sequence['sequence'] = sequence
        self._flush_uniparc_sequence()

    def _parse_uniparc_end(self, match):
        assert self._uniparc['uniparc_id'] is not None
        self._flush_uniparc()

    # === UniParc XRef ===
    def _parse_uniparc_xref(self, match):
        do_flush = False
        if match.endswith('/'):
            do_flush = True
            match = match.strip('/')
        #
        assert ('uniparc_id' not in self._uniparc_xref and
                self._uniparc['uniparc_id'] is not None)
        self._uniparc_xref['uniparc_id'] = self._uniparc['uniparc_id']
        #
        assert self._uniparc_xref['uniparc_xref_id'] is None
        self._uniparc_xref_id += 1
        self._uniparc_xref['uniparc_xref_id'] = self._uniparc_xref_id
        #
        data = self._parse_match(match)
        assert all(c in self._uniparc_xref_columns for c in data)
        self._uniparc_xref.update(data)
        if do_flush:
            self._flush_uniparc_xref()

    def _parse_uniparc_xref_end(self, match):
        assert 'uniparc_xref_id' in self._uniparc_xref
        self._flush_uniparc_xref()

    # === UniParc XRef Prop ===
    def _parse_uniparc_xref_prop(self, match):
        assert ('uniparc_xref_id' not in self._uniparc_xref_prop and
                self._uniparc_xref['uniparc_xref_id'] is not None)
        self._uniparc_xref_prop['uniparc_xref_id'] = self._uniparc_xref['uniparc_xref_id']
        #
        assert self._uniparc_xref_prop['uniparc_xref_prop_id'] is None
        self._uniparc_xref_prop_id += 1
        self._uniparc_xref_prop['uniparc_xref_prop_id'] = self._uniparc_xref_prop_id
        #
        data = self._parse_match(match)
        assert all(c in self._uniparc_xref_prop_columns for c in data)
        self._uniparc_xref_prop.update(data)
        self._flush_uniparc_xref_prop()

    # === Output ===
    def _flush_uniparc(self):
        if len(self._uniparc_cache) > 100:
            self._flush_cache()
        self._uniparc_cache.append(self._uniparc)
        self._uniparc = {'uniparc_id': None}

    def _flush_uniparc_sequence(self):
        self._uniparc_sequence_cache.append(self._uniparc_sequence)
        self._uniparc_sequence = {'uniparc_id': None}

    def _flush_uniparc_xref(self):
        self._uniparc_xref_cache.append(self._uniparc_xref)
        self._uniparc_xref = {'uniparc_xref_id': None}

    def _flush_uniparc_xref_prop(self):
        self._uniparc_xref_prop_cache.append(self._uniparc_xref_prop)
        self._uniparc_xref_prop = {'uniparc_xref_prop_id': None}

    def _flush_cache(self):
        """Flush cached data to files."""
        caches = [
            ('uniparc.tsv', self._uniparc_cache, self._uniparc_columns),
            ('uniparc_sequence.tsv', self._uniparc_sequence_cache, self._uniparc_sequence_columns),
            ('uniparc_xref.tsv', self._uniparc_xref_cache, self._uniparc_xref_columns),
            ('uniparc_xref_prop.tsv', self._uniparc_xref_prop_cache,
             self._uniparc_xref_prop_columns),
        ]
        for filename, cache, columns in caches:
            self._append_to_file(filename, cache, columns)
            del cache[:]

    def _append_to_file(self, filename, data, columns):
        if self.writer == 'pandas':
            self._append_to_file_pandas(filename, data, columns)
        elif self.writer == 'csv':
            self._append_to_file_csv(filename, data, columns)
        else:
            raise Exception

    def _append_to_file_pandas(self, filename, data, columns):
        import pandas as pd
        file_path = op.join(self.output_dir, filename)
        df = pd.DataFrame(data, columns=columns)
        write_header = False
        if not op.isfile(file_path):
            write_header = True
        with open(file_path, 'a+') as ofh:
            df.to_csv(ofh, sep='\t', na_rep='\\N', index=False, header=write_header)

    def _append_to_file_csv(self, filename, data, columns):
        import csv
        file_path = op.join(self.output_dir, filename)
        csv_writer_kwargs = {}
        if sys.version_info >= (3, 0):
            csv_writer_kwargs['newline'] = ''
        if not op.isfile(file_path):
            # Add header
            data.insert(0, {c: c for c in columns})
        with open(file_path, 'a+', **csv_writer_kwargs) as ofh:
            writer = csv.writer(
                ofh, delimiter='\t', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
            for row_dict in data:
                assert not set(row_dict) - set(columns), (
                    set(row_dict) - set(columns), filename)
                row = [row_dict.get(c, '\\N') for c in columns]
                writer.writerow(row)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--file_path', type=str)
    parser.add_argument('--output_dir', type=str)
    args = parser.parse_args()
    hacky_xml_parser = HackyXMLParser(args.file_path, args.output_dir)
    hacky_xml_parser.parse()
