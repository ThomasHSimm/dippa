"""Pseudo-Voigt peak profile functions.

Ported from ``pk_voigt2.m``, ``pk_voigt2asymm.m``, ``pk_alpha.m``,
``pk_alpha_asymm.m`` and ``pv_tv_aa.m`` in the original MATLAB tool
(github.com/ThomasHSimm/DPPA, ``dippa_v3``). Verified against the tool's own
bundled ``fit.mat``/``data.mat`` reference case: reproducing the saved fit as
a single asymmetric pseudo-Voigt per peak (no Kα doublet) gives R² = 0.994
against the real measured pattern. See ``AUDIT.md`` for the full write-up.

Peak parameter convention (one column per peak in the original ``aa`` array):

- 4 parameters: symmetric pseudo-Voigt   ``[x0, amplitude, fwhm, eta]``
- 6 parameters: asymmetric pseudo-Voigt  ``[x0, amplitude, fwhm_right,
  eta_left, fwhm_left, eta_right]``

The 6-parameter case cross-indexes FWHM and eta between the two sides — the
side with ``x < x0`` uses ``(fwhm_left, eta_left)`` and the side with
``x >= x0`` uses ``(fwhm_right, eta_right)``, but the parameter *positions* in
the array are ``[x0, amplitude, fwhm_right, eta_left, fwhm_left, eta_right]``,
not grouped by side. That pairing is preserved exactly as found in the
original source, not a design choice of this port — it is confirmed correct
against the reference fit, not merely faithful-by-assumption.

A separate, easy-to-get-backwards gotcha in the original code: the MATLAB
preference flag is named ``alpha2``, and ``alpha2 == 1`` means the doublet is
*not* fitted (single wavelength / no Kα2 component), while ``alpha2 == 0``
means the doublet *is* fitted. This port avoids that inversion by using an
explicit ``tube: str | None`` argument instead — ``None`` means no doublet,
a tube name means fit the Kα1/Kα2 doublet for that tube. Do not assume a
dataset's ``genset.mat``-style "alpha2" setting reflects how a *saved* fit
was actually produced — it reflects the GUI's state at whatever time that
file was last written, which is not necessarily fit time. The bundled
reference fit used here has ``alpha2 == 0`` in its ``genset.mat``, which
would imply "fit the doublet" — but the fit only reproduces the real pattern
(R² = 0.994) when no doublet is applied. Trust the reconstructed fit, not
the settings file, when the two disagree.

Coordinate convention (made explicit after review — see ``AUDIT.md`` §13):
``x`` is the diffraction-vector magnitude **g = 2 sin θ / λ = 1/d, in
Å⁻¹**, matching the original tool's reference data (whose peak positions
reproduce the fcc austenite 1/d sequence for a ≈ 3.6 Å). The Kα doublet
displacement below (``x0 * Δλ/λ̄``) is only valid in this coordinate —
passing degrees 2θ would silently produce a wrong doublet spacing, and
window widths like ``half_width=0.02`` are Å⁻¹ quantities. A typed pattern
object carrying coordinate/units/wavelength is planned (TODO.md); until
then this docstring is the contract.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]

# Kα1/Kα2 wavelengths in Å, as hardcoded in pk_alpha.m / pk_alpha_asymm.m.
# Only these three tubes are supported by the original tool; there is no
# path for an arbitrary user-supplied wavelength pair. Carried over as a
# known limitation (see AUDIT.md) rather than silently "fixed" here, so a
# port-parity check against the original stays meaningful.
TUBE_WAVELENGTHS: dict[str, tuple[float, float]] = {
    "Cu": (1.54056, 1.54439),
    "Co": (1.78897, 1.79285),
    "Fe": (1.93604, 1.93998),
}


def pseudo_voigt(x: FloatArray, x0: float, amplitude: float, fwhm: float, eta: float) -> FloatArray:
    """Symmetric pseudo-Voigt: ``amplitude * (eta * lorentzian + (1-eta) * gaussian)``.

    Ported from ``pk_voigt2.m``. ``eta`` is the Lorentzian fraction (0 =
    pure Gaussian, 1 = pure Lorentzian); the original GUI allows it up to
    1.2 ("super-Lorentzian"), which this function does not forbid either.
    """
    delx = x - x0
    gau = np.exp(-((delx / (0.600561 * fwhm)) ** 2))
    lor = 1.0 + (delx / (0.5 * fwhm)) ** 2
    return amplitude * ((eta / lor) + (1.0 - eta) * gau)


def asymmetric_pseudo_voigt(
    x: FloatArray,
    x0: float,
    amplitude: float,
    fwhm_right: float,
    eta_left: float,
    fwhm_left: float,
    eta_right: float,
) -> FloatArray:
    """Asymmetric (split) pseudo-Voigt: different width/mixing either side of x0.

    Ported from ``pk_voigt2asymm.m``. Parameter order matches the original
    ``aa`` row convention exactly (see module docstring) — the cross-indexing
    between which side uses which named parameter is intentional and
    confirmed against the reference fit, not a typo.
    """
    delx = x - x0
    out = np.empty_like(x, dtype=np.float64)
    left = delx < 0
    right = ~left

    dl = delx[left]
    gau_l = np.exp(-((dl / (0.600561 * fwhm_left)) ** 2))
    lor_l = 1.0 + (dl / (0.5 * fwhm_left)) ** 2
    out[left] = amplitude * ((eta_left / lor_l) + (1.0 - eta_left) * gau_l)

    dr = delx[right]
    gau_r = np.exp(-((dr / (0.600561 * fwhm_right)) ** 2))
    lor_r = 1.0 + (dr / (0.5 * fwhm_right)) ** 2
    out[right] = amplitude * ((eta_right / lor_r) + (1.0 - eta_right) * gau_r)

    return out


def _single_component(x: FloatArray, params: FloatArray) -> FloatArray:
    """Dispatch to symmetric or asymmetric pseudo-Voigt by parameter count."""
    if params.shape[0] == 4:
        return pseudo_voigt(x, *params)
    if params.shape[0] == 6:
        return asymmetric_pseudo_voigt(x, *params)
    raise ValueError(
        f"expected 4 (symmetric) or 6 (asymmetric) peak parameters, got {params.shape[0]}"
    )


def evaluate_peak(x: FloatArray, params: FloatArray, tube: str | None = None) -> FloatArray:
    """Evaluate one peak: a single pseudo-Voigt, or a Kα1/Kα2 doublet if ``tube`` is given.

    ``params`` is a length-4 (symmetric) or length-6 (asymmetric) array, using
    the same convention as one column of the original ``aa`` array. The
    second, weaker (half-intensity) component is placed at
    ``x0 * (lambda_alpha2 - lambda_alpha1) / mean(lambda_alpha1, lambda_alpha2)``
    beyond ``x0``, matching ``pk_alpha.m`` / ``pk_alpha_asymm.m`` exactly.
    """
    f1 = _single_component(x, params)
    if tube is None:
        return f1
    if tube not in TUBE_WAVELENGTHS:
        raise ValueError(f"unsupported tube {tube!r}; supported: {sorted(TUBE_WAVELENGTHS)}")
    lam1, lam2 = TUBE_WAVELENGTHS[tube]
    wavelen = (lam1 + lam2) / 2.0
    alphadiff = params[0] * (lam2 - lam1) / wavelen

    params2 = params.copy()
    params2[0] = params[0] + alphadiff
    params2[1] = 0.5 * params[1]
    f2 = _single_component(x, params2)
    return f1 + f2


def evaluate_pattern(x: FloatArray, aa: FloatArray, tube: str | None = None) -> FloatArray:
    """Evaluate a full multi-peak pattern plus shared quadratic background.

    Ported from ``pv_tv_aa.m``. ``aa`` has shape ``(rows, n_peaks + 1)``:
    columns 0..n_peaks-1 are one peak each (4 or 6 rows), and the last
    column holds the shared background ``[c0, c1, c2]`` such that
    ``background = c0 + c1*x + c2*x**2``. Rows beyond the third in the
    background column are unused padding in the original array and are
    ignored here too.
    """
    n_peaks = aa.shape[1] - 1
    total = np.zeros_like(x, dtype=np.float64)
    for i in range(n_peaks):
        total += evaluate_peak(x, aa[:, i], tube=tube)

    c0, c1, c2 = aa[0, -1], aa[1, -1], aa[2, -1]
    total += c0 + c1 * x + c2 * x**2
    return total


@dataclass(frozen=True)
class PeakFit:
    """A single fitted peak, in the parameter convention above."""

    x0: float
    amplitude: float
    params: FloatArray  # full 4- or 6-element array, for round-tripping to `aa`

    @property
    def is_asymmetric(self) -> bool:
        return self.params.shape[0] == 6
