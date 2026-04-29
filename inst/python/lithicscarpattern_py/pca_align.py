"""
pca_align.py
============
Area-weighted PCA alignment with energy-based sign convention.

Trimmed: removed clean_mesh, decimate_mesh, ICP, and all visualization.
Only robust_pca_alignment() remains - this is the core alignment routine.
"""

import numpy as np


def robust_pca_alignment(points: np.ndarray,
                         faces=None,
                         enforce_direction: bool = True,
                         verbose: bool = False):
    """
    Align a 3D point cloud so its three principal axes correspond to X, Y, Z.

    Sign convention (enforce_direction=True):
        For each axis, the positive direction is defined as the side carrying
        more squared projection energy. More stable for near-symmetric shapes
        than median-based rules.

    Area weighting (faces is not None):
        Each vertex is weighted by 1/3 of the total area of its adjacent
        faces, so high-curvature regions don't bias the principal axes.

    Returns
    -------
    aligned_points  : np.ndarray, shape (N, 3)
    rotation_matrix : np.ndarray, shape (3, 3)
    """
    if not isinstance(points, np.ndarray) or points.shape[1] != 3:
        raise ValueError("Input must be an Nx3 NumPy array")
    if len(points) < 3:
        raise ValueError("At least 3 points are required")

    # ---- Vertex weights ---------------------------------------------------
    if faces is not None:
        v0 = points[faces[:, 0]]
        v1 = points[faces[:, 1]]
        v2 = points[faces[:, 2]]
        face_areas     = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1)
        vertex_weights = np.zeros(len(points))
        for k in range(3):
            np.add.at(vertex_weights, faces[:, k], face_areas / 3.0)
        total = vertex_weights.sum()
        if total < 1e-12:
            vertex_weights = np.ones(len(points))
        else:
            vertex_weights /= total
    else:
        vertex_weights = np.ones(len(points)) / len(points)

    # ---- Weighted PCA ----------------------------------------------------
    centroid = (points * vertex_weights[:, None]).sum(axis=0)
    centered = points - centroid
    cov      = (centered * vertex_weights[:, None]).T @ centered

    _, _, Vt        = np.linalg.svd(cov)
    rotation_matrix = Vt.T

    if np.linalg.det(rotation_matrix) < 0:
        rotation_matrix[:, 2] *= -1

    # ---- Energy-based sign convention ------------------------------------
    if enforce_direction:
        projected = centered @ rotation_matrix
        for axis in range(3):
            proj  = projected[:, axis]
            w     = vertex_weights
            pos_e = np.dot(w[proj > 0], proj[proj > 0] ** 2) if np.any(proj > 0) else 0.0
            neg_e = np.dot(w[proj < 0], proj[proj < 0] ** 2) if np.any(proj < 0) else 0.0
            if neg_e > pos_e:
                rotation_matrix[:, axis] *= -1
        if np.linalg.det(rotation_matrix) < 0:
            rotation_matrix[:, 2] *= -1

    aligned_points = centered @ rotation_matrix

    if verbose:
        print(f"    PCA: det={np.linalg.det(rotation_matrix):.6f}, "
              f"area_weighted={faces is not None}")

    return aligned_points, rotation_matrix
