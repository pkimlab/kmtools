import codecs
import logging

import dill

logger = logging.getLogger(__name__)


def serialize_object_to_string(obj):
    """Serialize a Python object into a string.

    Useful for when you want to pass an object as a command-line
    argument to a file.

    Parameters
    ----------
    obj : object
        The object that you want to serialize.

    Returns
    -------
    data_hex_string : str
        String representation of the object.
    """
    data_bytes = dill.dumps(obj)
    data_hex = codecs.encode(data_bytes, "hex")
    data_hex_string = data_hex.decode('utf-8')
    return data_hex_string


def deserialize_object_from_string(data_hex_string):
    """Deserialize a Python object from a string.

    Parameters
    ----------
    data_hex_string : str
        String representation of the Python object.

    Returns
    -------
    obj : object
        The object that was encoded inside the string argument.
    """
    data_hex = data_hex_string.encode('utf-8')
    data_bytes = codecs.decode(data_hex, "hex")
    obj = dill.loads(data_bytes)
    return obj
