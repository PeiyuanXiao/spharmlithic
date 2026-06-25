"""
api.py
======
Top-level entry points called by the R package via reticulate.

Two pipelines:
    spharm_from_directions(df, ...)   - Track B (vMF KDE -> SPHARM)
    spharm_from_meshes(stl_dir, ...)  - Track A (mesh -> SPHARM)

Both return a dict { specimen_id: { 'coefficients': ndarray,
                                    'power_spectrum': ndarray } }.
"""

import os
import traceback
import numpy as np
import pandas as pd
import pyshtools as pysh

from . import spherical_kde
from . import kde_to_spharm
from . import spherical_harmonics


# =============================================================================
# Track B: directions -> SPHARM
# =============================================================================

def spharm_from_directions(df: pd.DataFrame,
                           lmax: int       = 20,
                           bandwidth: float = 0.35,
                           n_bearing: int  = 72,
                           n_plunge: int   = 36,
                           dh_size: int    = 64,
                           verbose: bool   = True) -> dict:
    """
    Pipeline:
        directions -> vMF KDE -> DH grid interpolation -> SPHARM expansion

    Parameters
    ----------
    df : DataFrame with columns ID, ux, uy, uz
        Already filtered/cleaned by the R side.
    lmax, bandwidth, n_bearing, n_plunge, dh_size : pipeline parameters

    Returns
    -------
    dict { id (str) : { 'coefficients': ndarray (2, lmax+1, lmax+1),
                        'power_spectrum': ndarray (lmax+1,) } }
    """
    if verbose:
        print(f"[spharm_from_directions] {df['ID'].nunique()} specimens, "
              f"{len(df)} direction vectors")

    # Step 1: KDE
    kde_result = spherical_kde.batch_spherical_kde(
        df,
        bandwidth = bandwidth,
        n_bearing = n_bearing,
        n_plunge  = n_plunge,
        verbose   = verbose,
    )

    # Build sphere_grid frame once (kde_to_spharm needs bearing/plunge cols)
    G = kde_result["G"]
    sphere_grid = pd.DataFrame({
        "x":       G[:, 0],
        "y":       G[:, 1],
        "z":       G[:, 2],
        "bearing": np.arctan2(G[:, 1], G[:, 0]),
        "plunge":  np.arcsin(np.clip(G[:, 2], -1, 1)),
    })

    # Step 2 & 3: DH interpolation + SH expansion per specimen
    if verbose:
        print(f"  SPHARM: lmax={lmax}, DH grid={dh_size}x{dh_size*2}")

    out = {}
    for i, sid in enumerate(kde_result["ids"]):
        try:
            grid_2d = kde_to_spharm.kde_vector_to_dh_grid(
                kde_result["kde_matrix"][i],
                sphere_grid,
                dh_size = dh_size,
            )
            feats = kde_to_spharm.compute_spharm_features(grid_2d, lmax=lmax)
            out[sid] = {
                "coefficients":   feats["coefficients"],
                "power_spectrum": feats["power_spectrum"],
            }
            if verbose:
                print(f"    [{i+1:>3}/{len(kde_result['ids'])}] {sid}  OK")
        except Exception as e:
            if verbose:
                print(f"    [{i+1:>3}/{len(kde_result['ids'])}] {sid}  FAILED: {e}")
            out[sid] = None

    return out


# =============================================================================
# Track A: meshes -> SPHARM
# =============================================================================

def spharm_from_meshes(stl_dir: str,
                       lmax: int                   = 20,
                       target_faces: int           = 20000,
                       grid_size: int              = 256,
                       smooth_iterations: int      = 3,
                       pre_decimate_threshold: int = 3_000_000,
                       pre_decimate_target: int    = 500_000,
                       verbose: bool               = True) -> dict:
    """
    Pipeline (per STL file):
        load -> decimate (open3d) -> smooth (trimesh Laplacian) ->
        normalize -> PCA align -> spherical interp -> SPHARM expansion

    Returns
    -------
    dict { specimen_id : { 'coefficients': ..., 'power_spectrum': ... } or None }
    """
    # Local imports so the core install (without trimesh/open3d) still works.
    from . import mesh_processing, pca_align

    stl_files = sorted(
        os.path.join(stl_dir, f)
        for f in os.listdir(stl_dir)
        if f.lower().endswith('.stl')
    )

    if verbose:
        print(f"[spharm_from_meshes] {len(stl_files)} STL files in {stl_dir}")
        print(f"  lmax={lmax}, target_faces={target_faces}, "
              f"grid_size={grid_size}, smooth={smooth_iterations}")

    out = {}
    for i, stl_path in enumerate(stl_files):
        sid = os.path.splitext(os.path.basename(stl_path))[0]
        if verbose:
            print(f"\n  [{i+1}/{len(stl_files)}] {sid}")

        try:
            # 1-3: load, decimate, smooth
            v, f, n_raw = mesh_processing.load_and_decimate(
                stl_path,
                target_faces           = target_faces,
                pre_decimate_threshold = pre_decimate_threshold,
                pre_decimate_target    = pre_decimate_target,
                verbose                = verbose,
            )
            v, f = mesh_processing.smooth_and_clean(
                v, f, iterations=smooth_iterations
            )

            # 4: normalize
            v_norm = mesh_processing.normalize_mesh(v, faces=f)

            # 5: PCA align (internal step, not exposed to R users)
            v_aligned, _ = pca_align.robust_pca_alignment(
                v_norm, faces=f, enforce_direction=True
            )

            # 6: spherical conversion + interpolation
            sph    = spherical_harmonics.cartesian_to_spherical(v_aligned)
            R, theta, phi = sph.T
            grid_r = spherical_harmonics.spherical_interpolate(
                R, theta, phi, grid_size
            )
            if grid_r is None:
                raise ValueError("spherical_interpolate returned None "
                                 "(too few points)")

            # 7: SH expansion (zero-component normalized)
            clm = spherical_harmonics.compute_spherical_harmonics(grid_r)

            # Pad/clip to lmax and compute power spectrum via SHCoeffs
            clm_sh = pysh.SHCoeffs.from_array(
                clm, normalization='4pi', csphase=1, lmax=lmax
            ).pad(lmax=lmax)

            out[sid] = {
                "coefficients":   np.asarray(clm_sh.coeffs),
                "power_spectrum": np.asarray(clm_sh.spectrum())[:lmax + 1],
            }
            if verbose:
                print(f"    OK")

        except Exception as e:
            if verbose:
                print(f"    FAILED: {e}")
                traceback.print_exc()
            out[sid] = None

    return out
