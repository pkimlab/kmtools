import pytest
from hypothesis import given
from hypothesis.strategies import characters, data, integers, lists, one_of, sets, text

from kmtools import py_tools


@given(data())
def test_struct(data):
    _elements = one_of(characters(), integers(), text())
    _sets = sets(elements=_elements, min_size=1, average_size=5)
    allowed_keys = data.draw(_sets)

    keys = list(allowed_keys)
    values = data.draw(lists(elements=_elements, min_size=len(keys), max_size=len(keys)))
    s = py_tools.Struct(allowed_keys)

    for key, value in zip(keys, values):
        s[key] = value
        assert s[key] == value

    assert s == {k: v for k, v in zip(keys, values)}

    prohibited_keys = {k for k in data.draw(_sets) if k not in allowed_keys}
    for key in prohibited_keys:
        with pytest.raises(KeyError):
            s[key] = data.draw(_elements)


@pytest.mark.parametrize("obj", [{'a': 10, 'b': 20}, 10, [1, 2, 3], range(100), "hello world"])
def test_object_serialization(obj):
    serialized = py_tools.serialize_object_to_string(obj)
    assert isinstance(serialized, str)
    deserialized = py_tools.deserialize_object_from_string(serialized)
    assert deserialized == obj
