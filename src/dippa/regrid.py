"""Regridding utilities: coordinate conversion and PCHIP interpolation.

Ported from BIGdippaFunctions/tsinterpl.m. Converts X-ray diffraction
patterns from 2θ (scattering angle) coordinates to g = 2sinθ/λ (reciprocal
space) with constant-step interpolation using PCHIP (Piecewise Cubic Hermite
Interpolation). Used in the original tool for data preprocessing before fitting.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]


def _regrid_pchip(x_data: FloatArray, y_data: FloatArray, x_grid: FloatArray) -> FloatArray:
    """Interpolate y-values to a regular grid using PCHIP.

    Parameters
    ----------
    x_data : ndarray
        Irregular x coordinates (must be monotonic)
    y_data : ndarray
        Y values corresponding to x_data
    x_grid : ndarray
        Regular target grid (sorted, ascending)

    Returns
    -------
    ndarray
        Interpolated y values at x_grid positions
    """
    if len(x_data) < 2:
        raise ValueError("need at least 2 data points for interpolation")

    if not np.all(np.diff(x_data) > 0) and not np.all(np.diff(x_data) < 0):
        raise ValueError("x_data must be monotonic")

    # Use monotone cubic hermite interpolation (scipy's PchipInterpolator)
    # This is equivalent to MATLAB's 'PCHIP' option
    from scipy.interpolate import PchipInterpolator

    pchip = PchipInterpolator(x_data, y_data)
    return pchip(x_grid)


def regrid_to_g(
    data: FloatArray,
    wavelength: float,
    g_step: float = 5e-5,
    sort_if_reversed: bool = True,
    clip_negative: bool = False,
) -> FloatArray:
    """Regrid diffraction pattern from 2θ to g = 2sinθ/λ with constant step.

    Converts X-ray diffraction data from scattering angle (2θ, in degrees)
    to reciprocal space (g = 2sinθ/λ, in Ų⁻¹) on a constant-step grid,
    using PCHIP interpolation. Matches the original tool's `tsinterpl.m`
    behavior exactly (see AUDIT.md §16).

    Parameters
    ----------
    data : ndarray, shape (n_points, 2)
        Input diffraction pattern: columns [2θ (degrees), intensity (counts)].
    wavelength : float
        X-ray wavelength in Ų. Standard values:
        - Cu Kα: ~1.541 Ų
        - Co Kα: ~1.790 Ų
    g_step : float, optional
        Target grid step size in Ų⁻¹. Default 5e-5 matches the reference
        fixture (ni_combo.mat).
    sort_if_reversed : bool, optional
        If True (default), sort by g if data is in descending order.
        Matches MATLAB behavior (line 15–20 in tsinterpl.m).
    clip_negative : bool, optional
        If True, clip negative interpolated values to zero. PCHIP can overshoot
        on low-intensity data near noise. Default False (matches MATLAB).

    Returns
    -------
    ndarray, shape (n_grid, 2)
        Regridded pattern: columns [g (Ų⁻¹), interpolated intensity].
        The g-grid is uniform with spacing g_step, starting from the
        lowest g value (floor) to the highest (no rounding to grid).

    Raises
    ------
    ValueError
        If wavelength is non-positive or data has < 2 points.

    Notes
    -----
    The conversion formula used is:
        g = 2 * sin(θ) / λ
    where θ = 2θ_degrees / 2 * π/180 is the half-angle in radians.

    Int16 overflow (MATLAB line 26) is automatically handled: Python ints
    are arbitrary precision, so no risk of overflow in the index calculation.
    """
    if wavelength <= 0:
        raise ValueError(f"wavelength must be positive, got {wavelength}")

    if data.shape[0] < 2:
        raise ValueError(f"need at least 2 data points, got {data.shape[0]}")

    # Extract columns
    theta_deg = data[:, 0].astype(np.float64)  # 2θ in degrees
    intensity = data[:, 1].astype(np.float64)

    # Convert 2θ (degrees) to θ (radians), then to g
    # g = 2 * sin(θ) / λ where θ = 2θ_degrees / 2 * π/180
    theta_rad = theta_deg * np.pi / 360.0  # 2θ/2 * π/180
    g = 2.0 * np.sin(theta_rad) / wavelength

    # Handle reversed/descending g (MATLAB lines 15–20)
    if sort_if_reversed and g[0] > g[-1]:
        # Sort by g (ascending)
        sort_idx = np.argsort(g)
        g = g[sort_idx]
        intensity = intensity[sort_idx]

    # Create regular grid (MATLAB line 22–31)
    # Round g_start down to nearest multiple of g_step
    g_min = g[0]
    remainder = g_min % g_step
    g_start = g_min - remainder  # Floor to grid
    g_end = g[-1]  # No rounding at end (MATLAB line 28)

    # Create grid: start:step:end (inclusive)
    n_points = int(np.floor((g_end - g_start) / g_step)) + 1
    g_grid = g_start + np.arange(n_points) * g_step

    # PCHIP interpolation (MATLAB line 35)
    intensity_interp = _regrid_pchip(g, intensity, g_grid)

    # Optionally clip negative values (PCHIP can overshoot on noisy data)
    if clip_negative:
        intensity_interp = np.maximum(intensity_interp, 0.0)

    # Return [g, interpolated_intensity]
    return np.column_stack([g_grid, intensity_interp])


def theta_to_g(theta_deg: FloatArray, wavelength: float) -> FloatArray:
    """Convert 2θ (degrees) to g = 2sinθ/λ (Ų⁻¹).

    Utility function for single-step conversion without interpolation.

    Parameters
    ----------
    theta_deg : ndarray
        Scattering angles in degrees (2θ)
    wavelength : float
        X-ray wavelength in Ų

    Returns
    -------
    ndarray
        Reciprocal space coordinates g in Ų⁻¹
    """
    theta_rad = theta_deg * np.pi / 360.0
    return 2.0 * np.sin(theta_rad) / wavelength


def g_to_theta(g: FloatArray, wavelength: float) -> FloatArray:
    """Convert g = 2sinθ/λ (Ų⁻¹) to 2θ (degrees).

    Inverse of theta_to_g. Utility for coordinate conversion.

    Parameters
    ----------
    g : ndarray
        Reciprocal space coordinates in Ų⁻¹
    wavelength : float
        X-ray wavelength in Ų

    Returns
    -------
    ndarray
        Scattering angles in degrees (2θ)
    """
    # g = 2*sin(θ)/λ  =>  sin(θ) = g*λ/2  =>  θ = arcsin(g*λ/2)
    # Then 2θ = 2*θ
    sin_theta = g * wavelength / 2.0
    # Clip to [-1, 1] to handle numeric errors near Bragg limit
    sin_theta = np.clip(sin_theta, -1.0, 1.0)
    theta_rad = np.arcsin(sin_theta)
    return 2.0 * theta_rad * 180.0 / np.pi
