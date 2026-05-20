"""
mesh_processing.py
==================
Mesh I/O, decimation, smoothing, and normalization.

Rewritten to use open3d + trimesh only (igl dependency removed).
"""

import os
import struct
import tempfile
import gc
import numpy as np


# =============================================================================
# I/O + decimation
# =============================================================================

def load_and_decimate(stl_path: str,
                      target_faces: int = 20000,
                      pre_decimate_threshold: int = 3_000_000,
                      pre_decimate_target: int = 500_000,
                      verbose: bool = True):
    """
    Load an STL file and decimate to `target_faces` using open3d's quadric
    error metric. Files exceeding `pre_decimate_threshold` faces are first
    pre-decimated by streaming face sub-sampling (binary STLs) or one-shot
    open3d decimation (ASCII STLs) to keep memory under control.

    Returns
    -------
    vertices : np.ndarray, shape (N, 3)
    faces    : np.ndarray, shape (M, 3)
    n_raw    : int   original face count of the source file
    """
    import open3d as o3d

    # ---- Inspect header to detect format and face count ------------------
    with open(stl_path, 'rb') as f:
        header_bytes = f.read(80)
        is_ascii = header_bytes.lstrip().startswith(b'solid')
        if is_ascii:
            n_raw = None
        else:
            n_raw = struct.unpack('<I', f.read(4))[0]

    tmp_path  = None
    load_path = stl_path
    o3d_hold  = None  # for ASCII fast path

    try:
        # ---- Optional pre-decimation -----------------------------------
        if is_ascii:
            if verbose:
                print(f"    ASCII STL detected, loading via open3d...")
            o3d_mesh = o3d.io.read_triangle_mesh(stl_path)
            n_raw = len(o3d_mesh.triangles)
            if verbose:
                print(f"    {n_raw:,} faces")
            if n_raw > pre_decimate_threshold:
                if verbose:
                    print(f"    pre-decimating to {pre_decimate_target:,}...")
                o3d_mesh = o3d_mesh.simplify_quadric_decimation(pre_decimate_target)
                o3d_mesh.compute_vertex_normals()
                tmp_path = tempfile.mktemp(suffix='.stl')
                o3d.io.write_triangle_mesh(tmp_path, o3d_mesh)
                del o3d_mesh
                gc.collect()
                load_path = tmp_path
            else:
                o3d_hold  = o3d_mesh
                load_path = None
        else:
            if verbose:
                print(f"    {n_raw:,} faces")
            if n_raw > pre_decimate_threshold:
                if verbose:
                    print(f"    pre-decimating to {pre_decimate_target:,} via streaming...")
                tmp_path = tempfile.mktemp(suffix='.stl')
                step     = max(1, n_raw // pre_decimate_target)
                keep_set = set(range(0, n_raw, step))
                n_keep   = len(keep_set)
                with open(stl_path, 'rb') as fin, open(tmp_path, 'wb') as fout:
                    header = fin.read(80)
                    fin.read(4)
                    fout.write(header)
                    fout.write(struct.pack('<I', n_keep))
                    for i in range(n_raw):
                        face_data = fin.read(50)
                        if i in keep_set:
                            fout.write(face_data)
                gc.collect()
                load_path = tmp_path

        # ---- Final decimation to target_faces ----------------------------
        if o3d_hold is not None:
            mesh = o3d_hold
        else:
            mesh = o3d.io.read_triangle_mesh(load_path)

        mesh = mesh.simplify_quadric_decimation(target_faces)
        vertices = np.asarray(mesh.vertices)
        faces    = np.asarray(mesh.triangles)
        if len(faces) == 0:
            raise ValueError("open3d.simplify_quadric_decimation returned empty mesh")

        if verbose:
            print(f"    decimation done: {len(faces):,} faces")

        return vertices, faces, n_raw

    finally:
        if tmp_path and os.path.exists(tmp_path):
            os.remove(tmp_path)


def smooth_and_clean(vertices: np.ndarray,
                     faces: np.ndarray,
                     iterations: int = 3):
    """
    Drop face indices that exceed the vertex array, then run trimesh
    Laplacian smoothing for `iterations` passes (skip if iterations <= 0).
    """
    import trimesh
    from trimesh.smoothing import filter_laplacian

    valid_mask = np.all(faces < len(vertices), axis=1)
    faces      = faces[valid_mask]

    mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=True)
    mesh.remove_unreferenced_vertices()
    if len(mesh.vertices) == 0 or len(mesh.faces) == 0:
        raise ValueError("Mesh empty after cleaning unreferenced vertices")

    if iterations and iterations > 0:
        filter_laplacian(mesh, iterations=int(iterations),
                         volume_constraint=False)

    return mesh.vertices, mesh.faces


# =============================================================================
# Normalization
# =============================================================================

def normalize_mesh(vertices: np.ndarray, faces=None) -> np.ndarray:
    """
    Center the mesh at its volume centroid and scale to fit in the unit sphere.

    Parameters
    ----------
    vertices : (N, 3) array
    faces : (M, 3) array, optional. If given, uses volume centroid (mass);
        otherwise falls back to vertex mean.
    """
    if faces is not None:
        import trimesh
        mesh     = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
        centroid = mesh.center_mass
    else:
        centroid = np.mean(vertices, axis=0)

    centered    = vertices - centroid
    max_radius  = np.max(np.linalg.norm(centered, axis=1))
    return centered / max_radius
