import json
import logging
import pickle
import zlib

import cloudpickle
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
    return zlib.compress(cloudpickle.dumps(obj), level=9).hex()


def deserialize_object_from_string(data_hex_string: str) -> object:
    """Deserialize a Python object from a string.

    Args:
        data_hex_string: String representation of the Python object.

    Returns:
        The object that was encoded inside the string argument.
    """
    return pickle.loads(zlib.decompress(bytes.fromhex(data_hex_string)))
