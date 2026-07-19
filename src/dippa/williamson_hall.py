"""Classical and modified Williamson-Hall fitting.

The ``mwhA``, ``mwhB`` and ``mwhC`` forms are ported verbatim from
``BIGdippaFunctions/getWH.m`` and ``dippa_fitWH.m`` in the original DPPA
tool. That implementation exposes four optimisation backends across its
branches (``lsqcurvefit``, pattern search, genetic algorithm and simulated
annealing). This port deliberately implements one inspectable backend,
``scipy.optimize.least_squares``, and reports its Jacobian-based uncertainty.

``classical`` uses the same mwhA path with q fixed to zero, making C constant
and X proportional to g. The original ``getWH.m`` has only the three modified
branches; this classical variant is an explicit extension in the Python API.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Sequence

import numpy as np
import scipy.optimize
import scipy.stats
from numpy.typing import NDArray

from dippa.breadth import DeltaKResult, ExcludedPeak
from dippa.contrast import Reflection, contrast_cubic, h_squared

FloatArray = NDArray[np.float64]
Variant = Literal["classical", "mwhA", "mwhB", "mwhC"]


@dataclass(frozen=True, slots=True)
class WHParameters:
    q: float
    size: float
    strain: float


@dataclass(frozen=True, slots=True)
class WHConfidenceIntervals:
    q: tuple[float, float]
    size: tuple[float, float]
    strain: tuple[float, float]


@dataclass(frozen=True, slots=True)
class WilliamsonHallResult:
    """Fit parameters, uncertainty, diagnostics and upstream exclusions."""

    variant: Variant
    parameters: WHParameters
    covariance: FloatArray
    confidence_intervals: WHConfidenceIntervals
    residuals: FloatArray
    fitted_delta_k: FloatArray
    x: FloatArray
    contrast: FloatArray
    n_points: int
    dof: int
    warnings: tuple[str, ...]
    excluded_peaks: tuple[ExcludedPeak, ...]
    success: bool
    message: str


def _model_from_x(x: FloatArray, size: float, strain: float, variant: Variant) -> FloatArray:
    if variant in {"classical", "mwhA"}:
        return size + strain * x
    if variant == "mwhB":
        return np.sqrt(size**2 + strain**2 * x**2)
    if variant == "mwhC":
        return size + strain * x**2
    raise ValueError(f"unknown Williamson-Hall variant: {variant!r}")


def williamson_hall_model(
    g: FloatArray,
    hkl: Sequence[Reflection] | NDArray[np.integer],
    ch00: float,
    parameters: WHParameters,
    variant: Variant,
) -> FloatArray:
    """Evaluate one getWH.m model with ``X = g * sqrt(C)``."""
    q = 0.0 if variant == "classical" else parameters.q
    contrast = contrast_cubic(hkl, ch00=ch00, q=q)
    if np.any(contrast <= 0):
        raise ValueError("contrast factors must be positive")
    x = np.asarray(g, dtype=np.float64) * np.sqrt(contrast)
    return _model_from_x(x, parameters.size, parameters.strain, variant)


def _confidence_intervals(
    values: FloatArray, covariance: FloatArray, dof: int, classical: bool
) -> WHConfidenceIntervals:
    if dof <= 0 or not np.all(np.isfinite(covariance)):
        intervals = np.full((3, 2), np.nan)
        if classical:
            intervals[0] = 0.0
    else:
        critical = float(scipy.stats.t.ppf(0.975, dof))
        half_width = critical * np.sqrt(np.maximum(np.diag(covariance), 0.0))
        intervals = np.column_stack((values - half_width, values + half_width))
        if classical:
            intervals[0] = 0.0
    return WHConfidenceIntervals(
        q=(float(intervals[0, 0]), float(intervals[0, 1])),
        size=(float(intervals[1, 0]), float(intervals[1, 1])),
        strain=(float(intervals[2, 0]), float(intervals[2, 1])),
    )


def fit_williamson_hall(
    breadths: DeltaKResult,
    hkl: Sequence[Reflection] | NDArray[np.integer],
    ch00: float,
    *,
    variant: Variant = "mwhA",
    initial: tuple[float, float, float] | None = None,
    q_bounds: tuple[float, float] | None = None,
) -> WilliamsonHallResult:
    """Fit q, size and strain jointly, or size/strain with q=0 classically.

    The HKL sequence is alongside the original peak array: retained peaks are
    selected using ``breadths.sample_peak_indices``. Reflections are never
    inferred from peak order, position or lattice symmetry.
    """
    if variant not in {"classical", "mwhA", "mwhB", "mwhC"}:
        raise ValueError(f"unknown Williamson-Hall variant: {variant!r}")
    indices = breadths.sample_peak_indices
    if len(hkl) <= int(indices.max(initial=-1)):
        raise ValueError("hkl must include an explicit reflection for every input peak")
    selected_hkl = [hkl[int(index)] for index in indices]
    h2 = h_squared(selected_hkl)
    n_parameters = 2 if variant == "classical" else 3
    n_points = len(breadths.delta_k)
    if n_points < n_parameters:
        raise ValueError(f"need at least {n_parameters} retained peaks, got {n_points}")

    if initial is None:
        initial = (min(2.0, 0.5 / max(float(h2.max()), 1e-12)), float(breadths.delta_k.min()), 0.01)
    q_initial, size_initial, strain_initial = (float(value) for value in initial)
    if q_bounds is None:
        positive_h2 = h2[h2 > 0]
        upper = 0.999999 / float(positive_h2.max()) if positive_h2.size else 100.0
        q_bounds = (-100.0, upper)

    if variant == "classical":
        p0 = np.maximum([size_initial, strain_initial], 1e-12)
        lower = np.array([0.0, 0.0])
        upper = np.array([np.inf, np.inf])

        def unpack(values: FloatArray) -> WHParameters:
            return WHParameters(0.0, float(values[0]), float(values[1]))

    else:
        p0 = np.array([q_initial, max(size_initial, 1e-12), max(strain_initial, 1e-12)])
        lower = np.array([q_bounds[0], 0.0, 0.0])
        upper = np.array([q_bounds[1], np.inf, np.inf])
        p0 = np.clip(p0, lower, upper)

        def unpack(values: FloatArray) -> WHParameters:
            return WHParameters(float(values[0]), float(values[1]), float(values[2]))

    def residuals(values: FloatArray) -> FloatArray:
        parameters = unpack(values)
        return (
            williamson_hall_model(
                breadths.positions, selected_hkl, ch00, parameters, variant
            )
            - breadths.delta_k
        )

    fit = scipy.optimize.least_squares(residuals, p0, bounds=(lower, upper))
    parameters = unpack(fit.x)
    fitted = williamson_hall_model(
        breadths.positions, selected_hkl, ch00, parameters, variant
    )
    residual = fitted - breadths.delta_k
    dof = n_points - n_parameters
    messages: list[str] = list(breadths.policy_messages)
    if dof <= 2:
        messages.append(f"low degrees of freedom: {dof}; confidence intervals are fragile")
    if not fit.success:
        messages.append(f"optimiser did not converge: {fit.message}")
    parameter_names = ("size", "strain") if variant == "classical" else ("q", "size", "strain")
    span = upper - lower
    tolerance = 1e-6 * np.where(np.isfinite(span), span, np.maximum(np.abs(fit.x), 1.0))
    at_lower = [
        name for name, value, bound, tol in zip(parameter_names, fit.x, lower, tolerance)
        if value - bound <= tol
    ]
    at_upper = [
        name for name, value, bound, tol in zip(parameter_names, fit.x, upper, tolerance)
        if bound - value <= tol
    ]
    if at_lower:
        messages.append(f"parameters at lower bound: {', '.join(at_lower)}")
    if at_upper:
        messages.append(f"parameters at upper bound: {', '.join(at_upper)}")

    covariance_fit = np.full((n_parameters, n_parameters), np.nan)
    if dof > 0:
        rank = int(np.linalg.matrix_rank(fit.jac))
        if rank < n_parameters:
            messages.append(
                f"rank-deficient Jacobian ({rank}/{n_parameters}); covariance is not identifiable"
            )
        else:
            variance = float(np.sum(residual**2) / dof)
            covariance_fit = variance * np.linalg.inv(fit.jac.T @ fit.jac)
    covariance = np.zeros((3, 3)) if variant == "classical" else covariance_fit
    if variant == "classical":
        covariance[1:, 1:] = covariance_fit
    values = np.array([parameters.q, parameters.size, parameters.strain])
    intervals = _confidence_intervals(values, covariance, dof, variant == "classical")
    contrast = contrast_cubic(
        selected_hkl, ch00=ch00, q=0.0 if variant == "classical" else parameters.q
    )
    x = breadths.positions * np.sqrt(contrast)
    return WilliamsonHallResult(
        variant=variant,
        parameters=parameters,
        covariance=covariance,
        confidence_intervals=intervals,
        residuals=residual,
        fitted_delta_k=fitted,
        x=x,
        contrast=contrast,
        n_points=n_points,
        dof=dof,
        warnings=tuple(messages),
        excluded_peaks=breadths.excluded_peaks,
        success=bool(fit.success),
        message=str(fit.message),
    )
