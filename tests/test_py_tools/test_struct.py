import pytest

from kmtools import py_tools


class TestStructFactory:

    def setup_method(self, method):
        Point = py_tools.struct_factory('Point', ('a', 'b'))
        self.point = Point()

    def test_set_attr(self):
        self.point.a = 10
        assert self.point['a'] == 10

    def test_set_key(self):
        self.point['b'] = 'hello world'
        assert self.point.b == 'hello world'

    def test_set_attr_missing(self):
        with pytest.raises(AttributeError):
            self.point.c = 30

    def test_set_key_missing(self):
        with pytest.raises(AttributeError):
            self.point['c'] = 30

    def test_contains(self):
        assert 'a' in self.point
        assert 'c' not in self.point

    def test_iter(self):
        assert list(self.point) == ['a', 'b']
        assert tuple(self.point) == ('a', 'b')

    @pytest.mark.xfail(reason="May be desirable, but don't know how to implement.")
    def test_dict(self):
        assert dict(self.point) == {'a': None, 'b': None}


@pytest.mark.parametrize("obj", [{'a': 10, 'b': 20}, 10, [1, 2, 3], range(100), "hello world"])
def test_object_serialization(obj):
    serialized = py_tools.serialize_object_to_string(obj)
    assert isinstance(serialized, str)
    deserialized = py_tools.deserialize_object_from_string(serialized)
    assert deserialized == obj
