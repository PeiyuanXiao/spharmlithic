"""
spherical_harmonics.py
======================
Cartesian -> spherical conversion, scattered -> regular grid interpolation,
and spherical harmonic expansion (DH sampling, 4pi normalization,
zero-component normalized).

Trimmed: removed clm_to_1d_standard, process_spherical_harmonics,
visualize_*, and the 'mean-radius' normalization branch.
"""

import numpy as np
import pyshtools.expand as shtools
from scipy.interpolate import griddata


def cartesian_to_spherical(vertices: np.ndarray) -> np.ndarray:
    """
    Cartesian -> spherical coordinates.

    Returns
    -------
    spherical : np.ndarray, shape (N, 3) with columns [r, theta, phi]
        theta in [0, pi]    colatitude
        phi   in [0, 2*pi)  longitude
    """
    x, y, z = vertices.T
    r       = np.sqrt(x**2 + y**2 + z**2)
    r       = np.where(r == 0, 1e-5, r)
    theta   = np.arccos(np.clip(z / r, -1.0, 1.0))
    phi     = np.arctan2(y, x) % (2 * np.pi)
    return np.column_stack([r, theta, phi])


def spherical_interpolate(R, theta, phi, grid_size: int):
    """
    Interpolate scattered (theta, phi, R) values onto a regular
    (grid_size, grid_size) grid.

    Returns
    -------
    grid : np.ndarray, shape (grid_size, grid_size), or None if too few points.
    """
    if len(R) < 4:
        return None

    I = np.linspace(0, np.pi,     grid_size, endpoint=False)
    J = np.linspace(0, 2 * np.pi, grid_size, endpoint=False)
    J, I = np.meshgrid(J, I)

    values = R
    points = np.array([theta, phi]).T

    # Polar completion
    points = np.concatenate((
        points,
        np.array([[0, 0], [0, 2 * np.pi], [np.pi, 0], [np.pi, 2 * np.pi]])
    ), axis=0)
    rmin = np.mean(R[theta == theta.min()])
    rmax = np.mean(R[theta == theta.max()])
    values = np.concatenate((values, [rmin, rmin, rmax, rmax]))

    # Longitude periodicity
    points = np.concatenate((points,
                             points - [0, 2 * np.pi],
                             points + [0, 2 * np.pi]), axis=0)
    values = np.concatenate((values, values, values))

    xi   = np.array([[I[i, j], J[i, j]]
                     for i in range(grid_size)
                     for j in range(grid_size)])
    grid = griddata(points, values, xi, method='linear')
    grid = grid.reshape((grid_size, grid_size))
    grid[:, -1] = grid[:, 0]
    return grid


def compute_spherical_harmonics(surface: np.ndarray) -> np.ndarray:
    """
    SH expansion of a 2D surface grid using DH sampling.

    Always uses zero-component normalization: all coefficients are divided
    by c(l=0, m=0), making them dimensionless relative to the mean radius.

    Parameters
    ----------
    surface : np.ndarray, shape (N, N) or (N, 2N). Both dims must be even.

    Returns
    -------
    harmonics : np.ndarray, shape (2, lmax+1, lmax+1)
        Complex coefficients in pyshtools '4pi' normalization, scaled so
        c(0, 0) == 1.
    """
    if surface.shape[1] % 2 or surface.shape[0] % 2:
        raise ValueError("Grid dimensions must be even")

    if surface.shape[1] == surface.shape[0]:
        sampling = 1
    elif surface.shape[1] == 2 * surface.shape[0]:
        sampling = 2
    else:
        raise ValueError("Grid must be (N, N) or (N, 2N)")

    harmonics = shtools.SHExpandDHC(surface.copy(), sampling=sampling)

    c00 = harmonics[0, 0, 0]
    if np.abs(c00) < 1e-10:
        raise ValueError("c(l=0, m=0) is near zero, cannot normalize")
    return harmonics / c00
