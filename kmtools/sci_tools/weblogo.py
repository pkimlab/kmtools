import io
import logging
from pathlib import Path
from typing import Any, List, Mapping
from unittest.mock import patch

import weblogo._cli
from PIL import Image

logger = logging.getLogger(__name__)


class _BytesIO(io.BytesIO):
    def __init__(self):
        super().__init__()
        self.buffer = self


def make_weblogo(seqs: List[str], **kwargs: Mapping[str, Any]):
    """Generate a logo for provided sequences.

    Args:
        seqs: List of sequences for which to generate a logo.
        units:
        color_scheme:
        stacks_per_line:
        format:
        output_file:
    """
    default_kwargs: Mapping[str, Any] = {
        "units": "bits",
        "sequence_type": "protein",
        "color_scheme": "charge",
        "stacks_per_line": 60,
        "format": "svg",
    }
    kwargs = {**default_kwargs, **kwargs}

    output_file = kwargs.pop("output_file")
    if not isinstance(output_file, (Path, str)):
        raise ValueError("output_file should be a Path or a string.")

    fin = io.StringIO()
    _write_sequences(seqs, fin)
    fin.seek(0)

    weblogo_args = ["weblogo"] + [
        f"--{k.replace('_', '-')}={str(v)}" for k, v in kwargs.items() if v is not None
    ]

    with patch("sys.stdin", fin), patch("weblogo._cli.sys.argv", weblogo_args), patch(
        "sys.stdout", new_callable=_BytesIO
    ) as patch_out:
        try:
            weblogo._cli.main()
        except RuntimeError as e:
            logger.error("Failed to create WebLogo image because of error: '%s'.", str(e))
            return None
        finally:
            patch_out.seek(0)
            img_data = patch_out.read()

    if output_file:
        with open(output_file, "wb") as fout:
            fout.write(img_data)

    if kwargs["format"] in ["eps", "png", "png_print", "jpeg"]:
        img = Image.open(io.BytesIO(img_data))
    else:
        img = None

    return img


def _write_sequences(seqs, fh):
    for i in range(len(seqs)):
        fh.write(f"> seq_{i}\n")
        fh.write(seqs[i] + "\n")
