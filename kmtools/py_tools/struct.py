import itertools
import logging

logger = logging.getLogger(__name__)


class Struct(dict):
    """A dictionary which allows only a restricted set of keys.

    Use this instead of a dictionary when you want to catch typos in names quickly!
    """

    # To avoid confusion, disable setting custom attributes
    __slots__ = ['_allowed_keys']

    def __init__(self, allowed_keys, *args, **kwargs):
        """
        Parameters
        ----------
        allowed_keys : set|list|tuple
            A collection of keys that the Struct can take.
        """
        self._allowed_keys = set(allowed_keys)
        self.update(*args, **kwargs)

    def __setitem__(self, key, value):
        """
        Examples
        --------
        >>> c = Struct({'a', 'b', 'c'})
        >>> c['a'] = 100
        >>> c['d'] = 100
        Traceback (most recent call last):
        KeyError:
        """
        if key not in self._allowed_keys:
            raise KeyError("The following key is not allowed: {}".format(repr(key)))
        super().__setitem__(key, value)

    def __getitem__(self, key):
        """
        Examples
        --------
        >>> c = Struct({'a', 'b', 'c'}, {'a': 100, 'b': 200})
        >>> c['a']
        100
        >>> c['c']
        >>> c['d']
        Traceback (most recent call last):
        KeyError:
        """
        if key in self:
            return super().__getitem__(key)
        elif key in self._allowed_keys:
            return None
        else:
            raise KeyError("The following key is not allowed: {}".format(repr(key)))

    def update(self, *args, **kwargs):
        """
        Examples
        --------
        >>> c = Struct({'a', 'b', 'c'})
        >>> c.update({'a': 'aaa', 'b': 'bbb'}, c='ccc')
        >>> c == {'a': 'aaa', 'b': 'bbb', 'c': 'ccc'}
        True
        """
        other = itertools.chain(*(d.items() for d in args) if args else {}.items(), kwargs.items())
        for key, value in other:
            self[key] = value
