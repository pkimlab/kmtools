import numpy as np
import numpy.linalg as la


def compute_angles(u0, u1, u2) -> np.ndarray:
    """Calculate the angles defined by three points.
    """
    v1 = u0 - u1
    # v1 /= la.norm(v1, axis=-1, keepdims=True)

    v2 = u2 - u1
    # v2 /= la.norm(v2, axis=-1, keepdims=True)

    x = (v1 * v2).sum(-1)
    y = np.cross(v1, v2)
    if y.ndim > 1:
        y = la.norm(y, axis=-1)

    return np.arctan2(y, x)


def compute_dihedrals(u0, u1, u2, u3) -> np.ndarray:
    """Calculate the dihedral angles defined by four points.
    """
    # https://math.stackexchange.com/a/47084/191500
    b1 = u1 - u0
    b2 = u2 - u1
    b3 = u3 - u2

    n1 = np.cross(b1, b2)
    # n1 /= la.norm(n1, axis=-1, keepdims=True)

    n2 = np.cross(b2, b3)
    # n2 /= la.norm(n2, axis=-1, keepdims=True)

    y = (b1 * n2).sum(-1) * la.norm(b2, axis=-1)
    x = (n1 * n2).sum(-1)

    return np.arctan2(y, x)
