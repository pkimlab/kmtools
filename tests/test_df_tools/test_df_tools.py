import numpy as np
import pytest
import tempfile
import kmtools.df_tools


class TestCVS:

    @classmethod
    def setup_class(cls):
        cls.df = kmtools.df_tools.random_df()

    @pytest.mark.parametrize(
        'sep, compression',
        [(sep, compression) for sep in (',', '\t') for compression in [None, '.gz', 'bz2', 'xz']])
    def test(self, sep, compression):
        df = self.df
        with tempfile.NamedTemporaryFile() as fh:
            file = fh.name
            # sep
            if sep == ',':
                file += '.csv'
            elif sep == '\t':
                file += '.tsv'
            else:
                file += '.txt'
            # compression
            if compression is not None:
                file += compression
            kmtools.df_tools.dump_csv(df, file)
            df_out = kmtools.df_tools.load_csv(file)
        assert (
            df.select_dtypes(include=[object]) == df_out.select_dtypes(include=[object])
        ).all().all()
        assert np.allclose(
            df.select_dtypes(exclude=[object]), df_out.select_dtypes(exclude=[object]))
