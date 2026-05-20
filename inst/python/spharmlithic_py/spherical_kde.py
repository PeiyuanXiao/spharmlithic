"""
spherical_kde.py
================
von Mises-Fisher kernel density estimation on the unit sphere.

Trimmed for the lithicscarpattern R package backend.
"""

import numpy as np
import pandas as pd


# =============================================================================
# Sphere grid
# =============================================================================

def make_sphere_grid(n_bearing: int = 72, n_plunge: int = 36) -> np.ndarray:
    """
    Build a uniform evaluation grid on the unit sphere.

    Parameters
    ----------
    n_bearing : int
        Number of azimuthal divisions (longitude). Default 72 -> 5 deg steps.
    n_plunge : int
        Number of elevation divisions (latitude). Default 36 -> 5 deg steps.

    Returns
    -------
    G : np.ndarray, shape (n_bearing * n_plunge, 3)
        Unit vectors for each grid point (x, y, z).
    """
    # Tiny epsilon at the poles to avoid sin(colat)=0 singularities downstream.
    eps     = np.deg2rad(0.01)
    bearing = np.linspace(0, 2 * np.pi, n_bearing, endpoint=False)
    plunge  = np.linspace(-np.pi / 2 + eps, np.pi / 2 - eps, n_plunge)

    b, p = np.meshgrid(bearing, plunge)
    x = np.cos(p) * np.cos(b)
    y = np.cos(p) * np.sin(b)
    z = np.sin(p)
    return np.column_stack([x.ravel(), y.ravel(), z.ravel()])


# =============================================================================
# von Mises-Fisher KDE
# =============================================================================

def fit_vmf_kde(
    ux: np.ndarray,
    uy: np.ndarray,
    uz: np.ndarray,
    G: np.ndarray,
    bandwidth: float = 0.35,
) -> np.ndarray:
    """
    Fit a vMF kernel density estimate for one specimen.

    Density at each grid point g:
        density(g) = mean_i [ exp(kappa * dot(g, x_i)) ]
    where kappa = 1 / bandwidth^2.

    Returns
    -------
    density : np.ndarray, shape (n_grid,)
        Normalised density values summing to 1.
    """
    X     = np.column_stack([ux, uy, uz]).astype(np.float64)
    kappa = 1.0 / bandwidth ** 2

    dot_mat = G @ X.T                                  # (n_grid, n_scars)
    density = np.mean(np.exp(kappa * dot_mat), axis=1)

    total = density.sum()
    if total == 0:
        raise ValueError("KDE density sums to zero - check input vectors.")
    return density / total


# =============================================================================
# Batch
# =============================================================================

def batch_spherical_kde(
    directions_df: pd.DataFrame,
    bandwidth: float = 0.35,
    n_bearing: int   = 72,
    n_plunge:  int   = 36,
    id_col: str = "ID",
    ux_col: str = "ux",
    uy_col: str = "uy",
    uz_col: str = "uz",
    verbose: bool = True,
) -> dict:
    """
    Run spherical KDE for every specimen in a DataFrame.

    Returns
    -------
    dict with keys:
        kde_matrix : np.ndarray (n_specimens, n_grid)
        G          : np.ndarray (n_grid, 3)
        ids        : list of specimen IDs (as strings)
    """
    G       = make_sphere_grid(n_bearing, n_plunge)
    n_grid  = len(G)
    all_ids = directions_df[id_col].astype(str).unique()
    n_cores = len(all_ids)

    kde_matrix = np.full((n_cores, n_grid), np.nan)

    if verbose:
        print(f"  vMF KDE: {n_cores} specimens, "
              f"bandwidth={bandwidth}, grid={n_bearing}x{n_plunge}")

    for i, id_i in enumerate(all_ids):
        df_i = directions_df[directions_df[id_col].astype(str) == id_i]
        kde_matrix[i] = fit_vmf_kde(
            df_i[ux_col].values,
            df_i[uy_col].values,
            df_i[uz_col].values,
            G,
            bandwidth,
        )
        if verbose:
            print(f"    [{i+1:>3}/{n_cores}] {id_i}  (n_scars={len(df_i)})")

    return {
        "kde_matrix": kde_matrix,
        "G":          G,
        "ids":        list(all_ids),
    }
