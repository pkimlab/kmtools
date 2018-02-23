import json

import numpy as np
from hypothesis import given
from hypothesis.extra.numpy import arrays
from hypothesis.strategies import floats, integers, lists, sets, text

from kmtools import py_tools


def serialization(obj):
    serialized = py_tools.serialize_object_to_string(obj)
    assert isinstance(serialized, str)
    return py_tools.deserialize_object_from_string(serialized)


@given(integers())
def test_serialization_integers(obj):
    deserialized = serialization(obj)
    assert deserialized == obj


@given(floats())
def test_serialization_floats(obj):
    deserialized = serialization(obj)
    assert deserialized == obj or (np.isnan(deserialized) and np.isnan(obj))


@given(text())
def test_serialization_text(obj):
    deserialized = serialization(obj)
    assert deserialized == obj


@given(lists(integers()))
def test_serialization_lists(obj):
    deserialized = serialization(obj)
    assert deserialized == obj


@given(sets(integers()))
def test_serialization_sets(obj):
    deserialized = serialization(obj)
    assert deserialized == obj


@given(arrays(np.float, (6, 8)))
def test_serialization_arrays(obj):
    deserialized = serialization(obj)
    assert np.allclose(deserialized, obj, equal_nan=True)


@given(arrays(np.float, (6, 8)))
def test_json_encoder_numpy(obj):
    serialized = json.dumps(obj, cls=py_tools.JSONEncoderNumPy)
    assert isinstance(serialized, str)
    deserialized = json.loads(serialized)
    assert np.allclose(deserialized, obj, equal_nan=True)
