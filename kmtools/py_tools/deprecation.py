import functools
import inspect
import warnings
from typing import Callable, Optional, Union


def deprecated(fn_or_reason: Union[Callable, str]):
    """
    This is a decorator which can be used to mark functions as deprecated.
    It will result in a warning being emitted when the function is used.

    Source: https://github.com/tantale/deprecated

    Args:
        reason: Reason message (or function/class/method to decorate).

    Returns:
        Decorated function that will emit a warning message when called.

    Examples:
        **Classic usage:**

        To use this, decorate your deprecated function with **@deprecated** decorator:

        >>> from deprecated import deprecated
        >>> @deprecated
        >>> def some_old_function(x, y):
            return x + y

        You can also decorate a class or a method:

        >>> from deprecated import deprecated
        >>> class SomeClass(object):
        >>> @deprecated
        >>> def some_old_method(self, x, y):
        ...     return x + y
        >>> @deprecated
        >>> class SomeOldClass(object):
        ...     pass

        You can give a "reason" message to help the developer to choose another function/class:

        >>> from deprecated import deprecated
        >>> @deprecated(reason="use another function")
        >>> def some_old_function(x, y):
        ...     return x + y
    """
    if isinstance(fn_or_reason, str):
        return functools.partial(_deprecate, reason=fn_or_reason)
    elif inspect.isfunction(fn_or_reason) or inspect.isclass(fn_or_reason):
        return _deprecate(fn_or_reason)
    else:
        raise TypeError(repr(type(fn_or_reason)))


def _deprecate(fn: Callable, reason: Optional[str] = None):
    if inspect.isclass(fn):
        message = f"Call to deprecated class {fn.__name__}."
    else:
        message = f"Call to deprecated function {fn.__name__}."

    if reason is not None:
        message += f' ({reason}).'

    @functools.wraps(fn)
    def new_fn(*args, **kwargs):
        warnings.simplefilter('always', DeprecationWarning)
        warnings.warn(
            message,
            category=DeprecationWarning,
            stacklevel=2
        )
        warnings.simplefilter('default', DeprecationWarning)
        return fn(*args, **kwargs)

    return new_fn
