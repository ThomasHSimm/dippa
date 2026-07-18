"""Background estimation.

Ported in spirit, not in mechanism, from `steel_bcg.m`/`bcgminus.m`: the
original identifies data points outside a window around every approximate
peak position and fits a quadratic to just those background-only points.
This module does the same point-selection, but solves the quadratic fit in
closed form (ordinary least squares) rather than replicating the original's
`fminsearch` started from `rand(1, 3)`.

Note on that original choice, corrected from an earlier draft of this
project's audit notes: fitting a polynomial by least squares is convex in
its parameters (a unique global minimum, no local-minima risk from a random
start), so the original approach isn't unreliable — this module is simpler
and avoids a solver-tolerance dependency, not a correctness fix. See
`AUDIT.md` §10.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]


def background_mask(x: FloatArray, peak_positions: FloatArray, half_width: float) -> FloatArray:
    """Boolean mask selecting points at least ``half_width`` from every peak position.

    Ported from the point-selection logic in ``bcgminus.m``: for each
    approximate peak position, exclude a window of ``half_width`` either
    side. What's left is treated as background-only signal.
    """
    mask = np.ones_like(x, dtype=bool)
    for pos in peak_positions:
        mask &= np.abs(x - pos) > half_width
    return mask


def fit_background_quadratic(
    x: FloatArray, y: FloatArray, peak_positions: FloatArray, half_width: float
) -> tuple[float, float, float]:
    """Fit ``c0 + c1*x + c2*x**2`` to the points outside every peak window.

    Closed-form ordinary least squares — see module docstring for why this
    doesn't need (and isn't helped by) the original's iterative nonlinear
    solver. Returns ``(c0, c1, c2)`` matching the ``aa`` background-column
    convention used throughout ``profiles.py``.

    Raises ``ValueError`` if fewer than 3 background points remain after
    excluding peak windows — the fit is underdetermined below that, and a
    silently-wrong quadratic is worse than an explicit failure here.
    """
    mask = background_mask(x, peak_positions, half_width)
    n_background_points = int(mask.sum())
    if n_background_points < 3:
        raise ValueError(
            f"only {n_background_points} background points remain after excluding peak "
            f"windows (need >= 3) — half_width={half_width} is likely too large relative "
            "to the spacing between peaks"
        )

    x_bg, y_bg = x[mask], y[mask]
    design = np.column_stack([np.ones_like(x_bg), x_bg, x_bg**2])
    coeffs, *_ = np.linalg.lstsq(design, y_bg, rcond=None)
    c0, c1, c2 = coeffs
    return float(c0), float(c1), float(c2)
