import tempfile
from pathlib import Path

import numpy as np
import pytest

from kmtools.sci_tools.weblogo import make_weblogo

try:
    import weblogo
except ImportError:
    weblogo = None


@pytest.mark.skipif(weblogo is None, reason="weblogo is not available")
def test_make_weblogo():
    n_seqs = 121
    seq_length = 18
    amino_acids = np.array(list("GVALICMFWPDESTYQNKRH"))
    assert len(amino_acids) == 20

    seqs = [
        "".join(seqs) for seqs in np.random.choice(amino_acids, size=(n_seqs, seq_length)).tolist()
    ]

    with tempfile.TemporaryDirectory() as tmp_dir:
        for format_ in ["png", "pdf", "svg"]:
            output_file = Path(tmp_dir, f"output.{format_}")
            assert not output_file.is_file()
            make_weblogo(seqs, format=format_, output_file=output_file)
            assert output_file.is_file()
