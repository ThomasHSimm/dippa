"""Per-peak fitting.

Ported in spirit from `onepeak.m` and its caller `indivfit_GUI.m`. The key
idea, worth stating plainly because it's not obvious from `onepeak.m` in
isolation (see AUDIT.md for the full trace): each peak is fit against a
*local, decontaminated* window of data, not the whole pattern with other
peaks frozen. "Decontaminated" means the current best-known contribution of
every other peak, plus the full background, has already been subtracted —
what's left in the window is (approximately) just this peak's own signal,
which is why the original only needs a small local *linear* correction on
top of it rather than re-fitting the whole quadratic background each time.

Documented simplification versus the original: this module fits one peak
at a time, even when its window overlaps a neighbour, and relies on
multiple passes (each pass sees the latest fit of every other peak) to let
close peaks refine each other iteratively. The original instead fits
overlapping peaks jointly in a shared window. For well-separated peaks
these converge to the same place; for closely overlapping clusters this is
a real, currently-unaddressed gap — see TODO.md.
"""

from __future__ import annotations

import numpy as np
import scipy.optimize
from numpy.typing import NDArray

from dippa.profiles import evaluate_pattern, evaluate_peak

FloatArray = NDArray[np.float64]


def _local_window(x: FloatArray, x0: float, half_width: float) -> FloatArray:
    """Boolean mask selecting points within ``half_width`` of ``x0``."""
    return np.abs(x - x0) <= half_width


def _decontaminate(x_window: FloatArray, aa: FloatArray, peak_index: int, tube: str | None) -> FloatArray:
    """The signal to subtract from raw data: every other peak, plus the background.

    ``evaluate_pattern(...) - evaluate_peak(..., aa[:, peak_index], ...)``
    is exactly "everything except this one peak" — all other peaks plus the
    shared quadratic background, using the *current* (not yet re-fit)
    values for everything else.
    """
    total = evaluate_pattern(x_window, aa, tube=tube)
    this_peak = evaluate_peak(x_window, aa[:, peak_index], tube=tube)
    return total - this_peak


def _bounds_for_single_peak(
    peak_params: FloatArray, half_width: float, data_range: float
) -> tuple[FloatArray, FloatArray]:
    """Bounds for one peak's parameters plus a local linear background.

    Loosely follows `onepeak.m`'s bounds for the peak actually being fit
    (see AUDIT.md §10) — position within `half_width`, amplitude within an
    order of magnitude, width capped relative to the window, eta in
    [0, 1.3]. Simplified relative to the original in that there's no
    "other peaks in this window, freeze them" case: this function only
    ever bounds the one peak being fit, because decontamination has
    already removed everyone else's contribution from the target signal.
    """
    x0, amplitude = peak_params[0], peak_params[1]
    amplitude = max(amplitude, 1e-7)  # guard against a zero/negative starting amplitude

    n_shape_params = len(peak_params) - 2  # fwhm/eta parameters (2 for symmetric, 4 for asymmetric)
    lb_shape = np.zeros(n_shape_params)
    ub_shape = np.full(n_shape_params, 1.3)
    # widths sit at even indices (0, 2, ...) within the shape parameters; eta at odd indices.
    for i in range(0, n_shape_params, 2):
        lb_shape[i] = 1e-6
        ub_shape[i] = half_width / 1.5

    lb = np.concatenate([[x0 - half_width, 0.1 * amplitude], lb_shape, [-data_range, -data_range]])
    ub = np.concatenate([[x0 + half_width, 10.0 * amplitude], ub_shape, [data_range, data_range]])
    return lb, ub


def fit_one_peak(
    x: FloatArray,
    y: FloatArray,
    aa: FloatArray,
    peak_index: int,
    half_width: float,
    tube: str | None = None,
) -> FloatArray:
    """Fit one peak against its local, decontaminated window. Returns an updated ``aa``.

    Everything except column ``peak_index`` is left untouched — this
    function only ever changes one peak's parameters (plus that peak's own
    local linear background stand-in, folded back into the peak's fit
    rather than the shared background column, matching the original's
    per-peak `aabcg`-style local correction).
    """
    window = _local_window(x, aa[0, peak_index], half_width)
    if window.sum() < len(aa[:, peak_index]) + 2:
        raise ValueError(
            f"only {int(window.sum())} data points fall within half_width={half_width} of "
            f"peak {peak_index} at x0={aa[0, peak_index]:.6g} — not enough to fit"
        )

    x_window = x[window]
    y_window = y[window]
    y_adjusted = y_window - _decontaminate(x_window, aa, peak_index, tube)

    peak_params = aa[:, peak_index].copy()
    data_range = float(y_window.max() - y_window.min()) or 1.0
    lb, ub = _bounds_for_single_peak(peak_params, half_width, data_range)

    p0 = np.concatenate([peak_params, [0.0, 0.0]])  # start local background at zero, per decontamination
    p0 = np.clip(p0, lb, ub)

    n_peak_params = len(peak_params)

    def residuals(params: FloatArray) -> FloatArray:
        peak = params[:n_peak_params]
        c0_local, c1_local = params[n_peak_params], params[n_peak_params + 1]
        model = evaluate_peak(x_window, peak, tube=tube) + c0_local + c1_local * x_window
        return model - y_adjusted

    result = scipy.optimize.least_squares(residuals, p0, bounds=(lb, ub))

    aa_new = aa.copy()
    aa_new[:, peak_index] = result.x[:n_peak_params]
    return aa_new


def fit_pattern(
    x: FloatArray,
    y: FloatArray,
    aa_initial: FloatArray,
    half_width: float,
    tube: str | None = None,
    n_passes: int = 2,
) -> FloatArray:
    """Fit every peak in ``aa_initial``, refining the shared quadratic background first.

    Stage 1 (background) is deliberately not repeated here — the quadratic
    term is estimated once from the initial peak positions (see
    `background.py`) and held fixed through per-peak fitting, matching the
    original's design (see AUDIT.md §10). Stage 2 sweeps every peak
    `n_passes` times; each pass sees the latest fit of every other peak,
    which is how overlapping peaks refine each other despite being fit one
    at a time (see module docstring for the documented gap this leaves for
    tightly overlapping clusters).
    """
    aa = aa_initial.copy()
    n_peaks = aa.shape[1] - 1
    for _ in range(n_passes):
        for peak_index in range(n_peaks):
            aa = fit_one_peak(x, y, aa, peak_index, half_width, tube=tube)
    return aa
