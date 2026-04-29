"""
kde_to_spharm.py
================
Interpolate KDE values onto a Driscoll-Healy grid and run spherical
harmonic expansion. Outputs raw coefficients and per-degree power.

Trimmed: no spectral entropy, no SHE, no variance analysis.
"""

import numpy as np
import pandas as pd
import pyshtools as pysh


# =============================================================================
# KDE vector -> Driscoll-Healy 2D grid
# =============================================================================

def kde_vector_to_dh_grid(kde_vec: np.ndarray,
                          sphere_grid: pd.DataFrame,
                          dh_size: int = 64) -> np.ndarray:
    """
    Interpolate a KDE probability vector onto a Driscoll-Healy regular grid.

    Parameters
    ----------
    kde_vec     : ndarray, shape (n_grid,)   single specimen's KDE values
    sphere_grid : DataFrame with columns 'bearing', 'plunge' (radians)
    dh_size     : DH latitude size; longitude has 2*dh_size points

    Returns
    -------
    grid_2d : ndarray, shape (dh_size, 2*dh_size), area-normalised
    """
    plunge  = sphere_grid["plunge"].values
    bearing = sphere_grid["bearing"].values

    colat_src = np.pi / 2 - plunge
    lon_src   = bearing

    n_lat    = dh_size
    n_lon    = 2 * dh_size
    colat_dh = np.linspace(0, np.pi,    n_lat, endpoint=False)
    lon_dh   = np.linspace(0, 2 * np.pi, n_lon, endpoint=False)
    TH, PH   = np.meshgrid(colat_dh, lon_dh, indexing='ij')

    tx = np.sin(TH) * np.cos(PH)
    ty = np.sin(TH) * np.sin(PH)
    tz = np.cos(TH)

    sx = np.sin(colat_src) * np.cos(lon_src)
    sy = np.sin(colat_src) * np.sin(lon_src)
    sz = np.cos(colat_src)

    dot = np.clip(
        tx[:, :, None] * sx + ty[:, :, None] * sy + tz[:, :, None] * sz,
        -1, 1,
    )
    weights = np.exp(50 * dot)
    grid_2d = (np.sum(weights * kde_vec, axis=2) /
               np.sum(weights, axis=2))

    grid_2d = np.clip(grid_2d, 0, None)

    sin_weights = np.sin(colat_dh)[:, None]
    area_sum    = (grid_2d * sin_weights).sum()
    if area_sum > 0:
        grid_2d /= area_sum

    return grid_2d


# =============================================================================
# DH grid -> SH coefficients + power spectrum
# =============================================================================

def compute_spharm_features(grid_2d: np.ndarray,
                            lmax: int = 20) -> dict:
    """
    Spherical harmonic expansion of a DH grid.

    Returns
    -------
    dict with:
        coefficients   : ndarray (2, lmax+1, lmax+1) - 4pi-normalized SH coeffs
        power_spectrum : ndarray (lmax+1,)           - raw power per degree
    """
    sh_grid = pysh.SHGrid.from_array(grid_2d, grid='DH')
    clm     = sh_grid.expand(lmax_calc=lmax)

    return {
        "coefficients":   np.asarray(clm.coeffs),
        "power_spectrum": np.asarray(clm.spectrum())[:lmax + 1],
    }
