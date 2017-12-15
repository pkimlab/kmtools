import codecs
import json
import logging

import dill
import numpy as np

logger = logging.getLogger(__name__)


class JSONEncoderNumPy(json.JSONEncoder):
    """Class for encoding numpy arrays as json documents

    Useful when you want to store numpy arrays inside the `results` attribute
    of a given class.
    """

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super().default(obj)


def serialize_object_to_string(obj: object) -> str:
    """Serialize a Python object into a string.

    Useful for when you want to pass an object as a command-line
    argument to a file.

    Args:
        obj: The object that you want to serialize.

    Returns:
        String representation of the object.
    """
    data_bytes = dill.dumps(obj)
    data_hex = codecs.encode(data_bytes, "hex")
    data_hex_string = data_hex.decode('utf-8')
    return data_hex_string


def deserialize_object_from_string(data_hex_string: str) -> object:
    """Deserialize a Python object from a string.

    Args:
        data_hex_string: String representation of the Python object.

    Returns:
        The object that was encoded inside the string argument.
    """
    data_hex = data_hex_string.encode('utf-8')
    data_bytes = codecs.decode(data_hex, "hex")
    obj = dill.loads(data_bytes)
    return obj
