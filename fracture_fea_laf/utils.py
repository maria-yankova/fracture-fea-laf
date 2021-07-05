import numpy as np

def polar2cart_2D(r, θ):
    """
    Convert 2D polar coordinates to Cartesian coordinates.

    """

    x = r * np.cos(θ)
    y = r * np.sin(θ)

    return x, y

def order_coplanar_points(points, normal, anticlockwise=True):
    """
    Find the clockwise or anticlockwise ordering of a set of coplanar 3D points.

    Parameters
    ----------
    points : ndarray of shape (3, N)
        The set of coplanar points (three-vector columns) whose ordering is to be found.
    normal : ndarray of shape (3, 1)
        Column three-vector representing the normal of the plane on which all points lie.

    Returns
    -------
    Ordered indices of points according to a clockwise or anticlockwise direction when
    looking in the opposite direction to `normal`.

    """

    # Normalise `normal` to a unit vector:
    normal = normal / np.linalg.norm(normal)

    # Compute the centroid:
    centroid = np.mean(points, axis=1)

    # Get the direction vectors from each point to the centroid
    p = points - centroid[:, np.newaxis]

    # Use the first point as the reference point
    # Find cross product of each point with the reference point
    crs = np.cross(p[:, 0:1], p, axis=0)

    # Find the scalar triple product of the point pairs and the normal vector
    stp = np.einsum('ij,ik->k', normal, crs)

    # Find the dot product of the point pairs
    dot = np.einsum('ij,ik->k', p[:, 0:1], p)

    # Find signed angles from reference point to each other point
    ang = np.arctan2(stp, dot)
    ang_order = np.argsort(ang)

    if not anticlockwise:
        ang_order = ang_order[::-1]

    return ang_order

def circle_points(r, n, centre):
    t = np.linspace(0, 2*np.pi, n)
    x = r * np.cos(t) + centre[0]
    y = r * np.sin(t) + centre[1]

    return np.c_[x, y]

def search_keep_order(A, B):
    sort_idx = A.argsort()
    out = sort_idx[np.searchsorted(A, B, sorter = sort_idx)]
    idx_layer_0 = np.nonzero(B[:,None] == A)[1]
    
    return idx_layer_0