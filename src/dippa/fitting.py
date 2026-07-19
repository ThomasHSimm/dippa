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

The local linear correction is fit jointly with the peak and *returned*, not
discarded — matching the original's ``aabcg`` design (`indivfit_GUI.m` seeds
each peak's fit from the stored ``aabcg`` column and stores the fitted values
back). Discarding them would let the local line absorb model mismatch during
optimisation and then vanish from the reported model; keeping them makes the
absorption inspectable. Note the same caveat as the original: ``aabcg`` is
*not* part of the global pattern model (``evaluate_pattern`` ignores it, in
both tools) — it exists only inside each peak's local fit.

Documented simplifications versus the original (see AUDIT.md §10–§11):

- The original fits peaks whose positions fall within the data window
  *jointly* (target peak loosely bounded; neighbours position-frozen but
  amplitude/width/eta loosely bounded). This module fits one peak at a time
  always, relying on multiple passes to let close peaks refine each other
  iteratively. Fine for well-separated peaks, an open gap for tightly
  overlapping clusters — see TODO.md.
- The original uses *two* preferences: ``sizfit`` for the data-window
  half-width and ``bcg2peak`` for the position bound / width cap inside
  ``onepeak.m``. This module uses a single ``half_width`` for both. The
  reference fixture's ``bcg2peak`` is 0.02; its ``sizfit`` is unknown (read
  live from the GUI, not stored — the same lesson as the ``alpha2``
  settings-file finding).

Every fit returns a structured result carrying optimiser success, cost,
bound-hit flags and the local background — see ``PeakFitResult`` /
``PatternFitResult``. A converged-looking array with a parameter silently
pinned at a bound is exactly the failure mode this exists to surface
(observed in practice: the reference pattern's weakest-amplitude peak pegs
``eta_right`` at the 1.3 bound — see AUDIT.md §15).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.optimize
from numpy.typing import NDArray

from dippa.profiles import evaluate_pattern, evaluate_peak

FloatArray = NDArray[np.float64]

_SYMMETRIC_NAMES = ("x0", "amplitude", "fwhm", "eta")
_ASYMMETRIC_NAMES = ("x0", "amplitude", "fwhm_right", "eta_left", "fwhm_left", "eta_right")
_LOCAL_BG_NAMES = ("c0_local", "c1_local")


def _param_names(n_peak_params: int) -> tuple[str, ...]:
    if n_peak_params == 4:
        return _SYMMETRIC_NAMES + _LOCAL_BG_NAMES
    if n_peak_params == 6:
        return _ASYMMETRIC_NAMES + _LOCAL_BG_NAMES
    raise ValueError(f"expected 4 or 6 peak parameters, got {n_peak_params}")


@dataclass(frozen=True)
class PeakFitResult:
    """Everything ``scipy.optimize.least_squares`` knew that a bare array would hide.

    ``success`` and the bound flags are the load-bearing fields: a fit that
    "converged" with a parameter pinned at its bound is not a trustworthy
    parameter estimate, even when the pattern-level R² looks fine (widths
    and mixing parameters can trade off against each other and against the
    local background while barely moving the residual — see AUDIT.md §15).
    """

    peak_index: int
    params: FloatArray
    local_background: tuple[float, float]
    success: bool
    message: str
    cost: float  # scipy's 0.5 * sum(residuals**2) over the local window
    nfev: int
    n_points: int
    at_lower_bound: tuple[str, ...] = field(default_factory=tuple)
    at_upper_bound: tuple[str, ...] = field(default_factory=tuple)

    @property
    def hit_bounds(self) -> bool:
        return bool(self.at_lower_bound or self.at_upper_bound)

    @property
    def clean(self) -> bool:
        """Converged and no parameter pinned at a bound."""
        return self.success and not self.hit_bounds


@dataclass(frozen=True)
class PatternFitResult:
    """Result of a full multi-pass pattern fit.

    ``aa`` follows the original array convention (peaks + shared quadratic
    background column). ``aabcg`` is the per-peak local linear background,
    shape ``(2, n_peaks)``, matching the original tool's ``aabcg`` variable
    exactly (see ``legacy_io.py``) — these are *not* included when
    evaluating ``aa`` with ``evaluate_pattern`` (true of the original tool
    too), so a large ``aabcg`` entry means the local fit absorbed signal
    the global model doesn't represent. Check ``warnings`` before using
    any width/shape parameter downstream.
    """

    aa: FloatArray
    aabcg: FloatArray
    peak_results: tuple[PeakFitResult, ...]  # final pass, one per peak
    n_passes: int

    @property
    def warnings(self) -> tuple[str, ...]:
        out: list[str] = []
        for r in self.peak_results:
            if not r.success:
                out.append(f"peak {r.peak_index}: optimiser did not converge ({r.message})")
            if r.at_lower_bound:
                out.append(f"peak {r.peak_index}: at lower bound: {', '.join(r.at_lower_bound)}")
            if r.at_upper_bound:
                out.append(f"peak {r.peak_index}: at upper bound: {', '.join(r.at_upper_bound)}")
        return tuple(out)

    @property
    def all_clean(self) -> bool:
        return all(r.clean for r in self.peak_results)


def _local_window(x: FloatArray, x0: float, half_width: float) -> FloatArray:
    """Boolean mask selecting points within ``half_width`` of ``x0``."""
    return np.abs(x - x0) <= half_width


def _decontaminate(
    x_window: FloatArray, aa: FloatArray, peak_index: int, tube: str | None
) -> FloatArray:
    """The signal to subtract from raw data: every other peak, plus the background.

    ``evaluate_pattern(...) - evaluate_peak(..., aa[:, peak_index], ...)``
    is exactly "everything except this one peak" — all other peaks plus the
    shared quadratic background, using the *current* (not yet re-fit)
    values for everything else. Matches `indivfit_GUI.m`'s
    ``I + pv_tv_aa(nearby, ...) - pv_tv_aa(everything, ...)`` identity.
    Other peaks' ``aabcg`` corrections are not subtracted — the original
    doesn't either (``pv_tv_aa`` never sees ``aabcg``).
    """
    total = evaluate_pattern(x_window, aa, tube=tube)
    this_peak = evaluate_peak(x_window, aa[:, peak_index], tube=tube)
    return total - this_peak


def _bounds_for_single_peak(
    peak_params: FloatArray, half_width: float, data_range: float
) -> tuple[FloatArray, FloatArray]:
    """Bounds for one peak's parameters plus a local linear background.

    Follows `onepeak.m`'s bounds for the *target* peak (see AUDIT.md §10):
    position within ``half_width`` (the original uses ``bcg2peak``),
    amplitude within an order of magnitude, widths capped at
    ``half_width / 1.5`` (original: ``bcg2peak/1.5``), eta in [0, 1.3]
    (original literal). The local background is effectively unbounded in
    the original (``-inf``/``inf``); here it's bounded by ``±data_range``,
    which is loose enough not to bind in practice but keeps the solver's
    trust region sane.
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


def _bound_hits(
    x: FloatArray, lb: FloatArray, ub: FloatArray, names: tuple[str, ...]
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Which parameters landed at (or numerically against) their bounds."""
    span = ub - lb
    tol = 1e-6 * np.where(np.isfinite(span), span, 1.0)
    at_lb = tuple(n for n, hit in zip(names, x - lb <= tol) if hit)
    at_ub = tuple(n for n, hit in zip(names, ub - x <= tol) if hit)
    return at_lb, at_ub


def fit_one_peak(
    x: FloatArray,
    y: FloatArray,
    aa: FloatArray,
    peak_index: int,
    half_width: float,
    tube: str | None = None,
    local_background_start: tuple[float, float] = (0.0, 0.0),
) -> tuple[FloatArray, PeakFitResult]:
    """Fit one peak against its local, decontaminated window.

    Returns ``(aa_new, result)``: the updated parameter array (only column
    ``peak_index`` changes) and a `PeakFitResult` carrying the fitted local
    linear background plus optimiser diagnostics. The local background is
    seeded from ``local_background_start`` — pass the previous fit's values
    (the ``aabcg`` convention) when iterating, as `indivfit_GUI.m` does.
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

    p0 = np.concatenate([peak_params, list(local_background_start)])
    p0 = np.clip(p0, lb, ub)

    n_peak_params = len(peak_params)
    names = _param_names(n_peak_params)

    def residuals(params: FloatArray) -> FloatArray:
        peak = params[:n_peak_params]
        c0_local, c1_local = params[n_peak_params], params[n_peak_params + 1]
        model = evaluate_peak(x_window, peak, tube=tube) + c0_local + c1_local * x_window
        return model - y_adjusted

    result = scipy.optimize.least_squares(residuals, p0, bounds=(lb, ub))
    at_lb, at_ub = _bound_hits(result.x, lb, ub, names)

    aa_new = aa.copy()
    aa_new[:, peak_index] = result.x[:n_peak_params]

    fit_result = PeakFitResult(
        peak_index=peak_index,
        params=result.x[:n_peak_params].copy(),
        local_background=(float(result.x[n_peak_params]), float(result.x[n_peak_params + 1])),
        success=bool(result.success),
        message=str(result.message),
        cost=float(result.cost),
        nfev=int(result.nfev),
        n_points=int(window.sum()),
        at_lower_bound=at_lb,
        at_upper_bound=at_ub,
    )
    return aa_new, fit_result


def fit_pattern(
    x: FloatArray,
    y: FloatArray,
    aa_initial: FloatArray,
    half_width: float,
    tube: str | None = None,
    n_passes: int = 2,
) -> PatternFitResult:
    """Fit every peak in ``aa_initial``, refining the shared quadratic background first.

    Stage 1 (background) is deliberately not repeated here — the quadratic
    term is estimated once from the initial peak positions (see
    `background.py`) and held fixed through per-peak fitting, matching the
    original's design (see AUDIT.md §10). Stage 2 sweeps every peak
    `n_passes` times; each pass sees the latest fit of every other peak,
    which is how overlapping peaks refine each other despite being fit one
    at a time (see module docstring for the documented gap this leaves for
    tightly overlapping clusters). Per-peak local backgrounds (``aabcg``)
    persist across passes, seeding each peak's next fit — as
    `indivfit_GUI.m` does.

    Check ``result.warnings`` before using width/shape parameters
    downstream (Williamson-Hall consumes breadth, and breadth is exactly
    what a bound-pegged fit gets wrong quietly).
    """
    aa = aa_initial.copy()
    n_peaks = aa.shape[1] - 1
    aabcg = np.zeros((2, n_peaks))
    last_results: list[PeakFitResult | None] = [None] * n_peaks

    for _ in range(n_passes):
        for peak_index in range(n_peaks):
            aa, peak_result = fit_one_peak(
                x,
                y,
                aa,
                peak_index,
                half_width,
                tube=tube,
                local_background_start=(aabcg[0, peak_index], aabcg[1, peak_index]),
            )
            aabcg[:, peak_index] = peak_result.local_background
            last_results[peak_index] = peak_result

    return PatternFitResult(
        aa=aa,
        aabcg=aabcg,
        peak_results=tuple(r for r in last_results if r is not None),
        n_passes=n_passes,
    )
