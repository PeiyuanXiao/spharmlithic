"""
kde_to_spharm.py
================
Interpolate KDE values onto a Driscoll-Healy grid, run spherical harmonic
expansion (forward), and reconstruct from coefficients (inverse).

Outputs raw coefficients and per-degree power. No spectral entropy / SHE /
variance analysis - those are left to the R side.
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
# DH grid -> SH coefficients + power spectrum (forward transform)
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


# =============================================================================
# SH coefficients -> DH grid (inverse transform)
# =============================================================================

def reconstruct_from_coeffs(coeffs: np.ndarray,
                            grid_size: int = 64) -> np.ndarray:
    """
    Inverse spherical harmonic transform.

    Given a coefficient array of shape (2, lmax+1, lmax+1), reconstruct
    the corresponding scalar field on a Driscoll-Healy grid. Negative
    values introduced by finite-degree truncation are clipped to zero.

    Parameters
    ----------
    coeffs    : ndarray, shape (2, lmax+1, lmax+1)  4pi-normalized SH coeffs.
                Real-valued. (R side strips imaginary parts before calling.)
    grid_size : int. Latitude resolution of the output grid; longitude
                resolution is 2 * grid_size.

    Returns
    -------
    grid_2d : ndarray, shape (grid_size, 2*grid_size).
              Non-negative reconstructed density.

    Notes
    -----
    The output grid size is determined by `grid_size`, NOT by the lmax
    of the input coefficients. pyshtools internally upsamples the SH
    expansion onto the requested grid.
    """
    coeffs = np.asarray(coeffs, dtype=np.float64)
    if coeffs.ndim != 3 or coeffs.shape[0] != 2:
        raise ValueError(
            f"coeffs must have shape (2, lmax+1, lmax+1); got {coeffs.shape}"
        )

    clm = pysh.SHCoeffs.from_array(coeffs, normalization='4pi', csphase=1)

    # extend='False' so output shape is exactly (grid_size, 2*grid_size).
    sh_grid = clm.expand(grid='DH', lmax=grid_size - 1, extend=False)
    grid_2d = np.asarray(sh_grid.to_array(), dtype=np.float64)

    # Clip negative values from truncation noise.
    grid_2d = np.clip(grid_2d, 0, None)

    return grid_2d
