from pathlib import Path
from typing import List

import pytest
import yaml

from kmtools import sequence_tools

# --- Load test data ---


def make_kwargs_constructor(obj):
    def kwargs_constructor(loader, node):
        values = loader.construct_mapping(node)
        return obj(**values)

    return kwargs_constructor


yaml.add_constructor("!HHRAlignment", make_kwargs_constructor(sequence_tools.HHRAlignment))

with Path(__file__).with_suffix(".yaml").open("rt") as fin:
    TEST_DATA = yaml.load(fin)


# --- Tests ---


@pytest.mark.parametrize("data, results", TEST_DATA["test_parse_hhr_data"], ids=["0"])
def test_parse_hhr_data(data: str, results: List[sequence_tools.HHRAlignment]):
    results_ = sequence_tools.parse_hhr_data(data.strip().split("\n"))
    for r, r_ in zip(results, results_):
        assert r == r_, (r, r_)
