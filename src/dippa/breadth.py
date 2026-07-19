"""Peak breadths and instrumental subtraction.

Breadth definitions and nearest-position instrumental matching port
``BIGdippaFunctions/getFW_IB.m`` from the original DPPA tool. Instrumental
breadths are subtracted linearly (``K_R - K_I``), which assumes Lorentzian
addition of breadths. No quadratic/Gaussian subtraction alternative is
implemented yet.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal
import warnings as warnings_module

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]
BreadthKind = Literal["IB", "FW"]


@dataclass(frozen=True, slots=True)
class ExcludedPeak:
    """A sample peak omitted from instrumental subtraction or WH fitting."""

    peak_index: int
    reasons: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class DeltaKResult:
    """Instrument-corrected breadths plus every policy decision."""

    delta_k: FloatArray
    positions: FloatArray
    sample_peak_indices: NDArray[np.int64]
    instrument_peak_indices: NDArray[np.int64]
    excluded_peaks: tuple[ExcludedPeak, ...]
    policy_messages: tuple[str, ...]
    kind: BreadthKind


def _parameter_array(value: object, label: str) -> FloatArray:
    array = np.asarray(getattr(value, "aa", value), dtype=np.float64)
    if array.ndim != 2 or array.shape[0] not in {4, 6} or array.shape[1] < 2:
        raise ValueError(f"{label} aa must have shape (4 or 6, n_peaks + 1)")
    return array


def peak_breadths(aa: FloatArray, kind: BreadthKind = "IB") -> FloatArray:
    """Return per-peak FW or IB, ported from ``getFW_IB.m``.

    For asymmetric peaks, FWHM and eta are each side-averaged before the
    integral-breadth expression is evaluated. ``IB`` is the default because
    it is the better-identified downstream quantity in this port.
    """
    parameters = _parameter_array(aa, "peak")
    if kind not in {"IB", "FW"}:
        raise ValueError(f"kind must be 'IB' or 'FW', got {kind!r}")
    if parameters.shape[0] == 4:
        fwhm = parameters[2, :-1]
        eta = parameters[3, :-1]
    else:
        fwhm = 0.5 * (parameters[2, :-1] + parameters[4, :-1])
        eta = 0.5 * (parameters[3, :-1] + parameters[5, :-1])
    if kind == "FW":
        return fwhm.copy()
    return 0.5 * fwhm * (np.pi * eta + (1.0 - eta) * np.sqrt(np.pi / np.log(2.0)))


def _flagged_peaks(value: object) -> dict[int, tuple[str, ...]]:
    results = getattr(value, "peak_results", ())
    flagged: dict[int, tuple[str, ...]] = {}
    for result in results:
        reasons: list[str] = []
        if not result.success:
            reasons.append("optimiser did not converge")
        if result.at_lower_bound:
            reasons.append(f"at lower bound: {', '.join(result.at_lower_bound)}")
        if result.at_upper_bound:
            reasons.append(f"at upper bound: {', '.join(result.at_upper_bound)}")
        if reasons:
            flagged[int(result.peak_index)] = tuple(reasons)
    return flagged


def delta_k(
    sample: object,
    instrument: object,
    *,
    kind: BreadthKind = "IB",
    include_flagged: bool = False,
    drop_nonpositive: bool = False,
) -> DeltaKResult:
    """Subtract nearest-position instrumental breadths from sample breadths.

    Flagged sample or matched-instrument peaks are excluded by default;
    ``include_flagged=True`` explicitly overrides that safety policy. A
    nonpositive corrected breadth raises by default. With
    ``drop_nonpositive=True`` it is excluded and a Python warning is emitted.
    Both policies are recorded in ``excluded_peaks`` and ``policy_messages``.
    """
    sample_aa = _parameter_array(sample, "sample")
    instrument_aa = _parameter_array(instrument, "instrument")
    sample_widths = peak_breadths(sample_aa, kind)
    instrument_widths = peak_breadths(instrument_aa, kind)
    positions = sample_aa[0, :-1]
    instrument_positions = instrument_aa[0, :-1]
    nearest = np.argmin(np.abs(positions[:, np.newaxis] - instrument_positions), axis=1)

    sample_flagged = _flagged_peaks(sample)
    instrument_flagged = _flagged_peaks(instrument)
    exclusion_reasons: dict[int, list[str]] = {}
    messages: list[str] = []
    if not include_flagged:
        for peak_index, reasons in sample_flagged.items():
            exclusion_reasons.setdefault(peak_index, []).extend(reasons)
        for peak_index, instrument_index in enumerate(nearest):
            if int(instrument_index) in instrument_flagged:
                reasons = (
                    f"matched instrument peak {instrument_index}: {reason}"
                    for reason in instrument_flagged[int(instrument_index)]
                )
                exclusion_reasons.setdefault(peak_index, []).extend(reasons)
        if exclusion_reasons:
            messages.append(
                "excluded fit-flagged peaks: "
                + ", ".join(str(index) for index in sorted(exclusion_reasons))
            )

    corrected = sample_widths - instrument_widths[nearest]
    active = np.array(
        [index for index in range(len(positions)) if index not in exclusion_reasons], dtype=int
    )
    nonpositive = active[corrected[active] <= 0]
    if nonpositive.size and not drop_nonpositive:
        listed = ", ".join(str(index) for index in nonpositive)
        raise ValueError(f"nonpositive delta_k for sample peak(s): {listed}")
    if nonpositive.size:
        listed = ", ".join(str(index) for index in nonpositive)
        message = f"dropped nonpositive delta_k for sample peak(s): {listed}"
        warnings_module.warn(message, UserWarning, stacklevel=2)
        messages.append(message)
        for peak_index in nonpositive:
            exclusion_reasons.setdefault(int(peak_index), []).append("nonpositive delta_k")

    included = np.array(
        [index for index in range(len(positions)) if index not in exclusion_reasons], dtype=int
    )
    excluded = tuple(
        ExcludedPeak(index, tuple(exclusion_reasons[index])) for index in sorted(exclusion_reasons)
    )
    return DeltaKResult(
        delta_k=corrected[included],
        positions=positions[included],
        sample_peak_indices=included.astype(np.int64),
        instrument_peak_indices=nearest[included].astype(np.int64),
        excluded_peaks=excluded,
        policy_messages=tuple(messages),
        kind=kind,
    )
