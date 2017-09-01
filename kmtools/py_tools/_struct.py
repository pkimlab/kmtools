import logging

logger = logging.getLogger(__name__)


def struct_factory(name, slots):
    def __setitem__(self, key, value):
        if isinstance(key, str):
            return setattr(self, key, value)
        else:
            return setattr(self, self.__slots__[key], value)

    def __getitem__(self, key):
        if isinstance(key, str):
            return getattr(self, key)
        else:
            return self.__slots__[key]

    variables = {'__slots__': tuple(slots), '__setitem__': __setitem__, '__getitem__': __getitem__}
    return type(name, (), variables)
