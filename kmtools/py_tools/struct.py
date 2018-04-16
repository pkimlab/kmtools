import logging

logger = logging.getLogger(__name__)


def struct_factory(name, slots):
    """
    Notes:
        - This is deprecated. Used `typing.NamedTuple` / `attr.s` instead.
    """

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

    def update(self, kwargs):
        for key, value in kwargs.items():
            self[key] = value

    def values(self):
        for key in self:
            yield self[key]

    def items(self):
        for key in self:
            yield key, self[key]

    methods = {
        '__slots__': tuple(slots),
        '__setitem__': __setitem__,
        '__getitem__': __getitem__,
        'update': update,
        'values': values,
        'items': items
    }
    return type(name, (), methods)
