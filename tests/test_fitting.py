"""Tests for dippa.fitting.

The important test here is `test_fit_pattern_recovers_real_pattern_from_rough_guess`
— it starts from a starting point that is *worse than predicting the mean*
(R^2 < 0) and checks the staged fitter finds its way back to the real
answer. Evaluating known-good parameters (as profiles.py's tests do) proves
the physics is right; this proves the *fitter* can find an answer.

Deliberate framing, after a review pass (see AUDIT.md §15): this is a
**single-fixture recovery test, not MATLAB parity**. The starting guess is
a perturbation of the known answer (positions jittered within the window,
amplitudes scaled by a factor drawn from [0.5, 1.5] of the true value,
widths/etas replaced by generic placeholders) — not an independent starting
point — and recovery is asserted on R², positions, amplitudes and
*integral breadth*, not on every parameter matching MATLAB's saved values.
FWHM and eta individually are not fully identified by this data (they trade
off against each other; one peak pins eta at its bound — pinned by
`test_fit_pattern_reports_bound_hits`, so the behaviour is visible rather
than silently assumed away). Genuine parity would need MATLAB and Python
runs from identical starts with per-parameter tolerances — see TODO.md.
"""

from pathlib import Path

import numpy as np
import pytest
import scipy.io

from dippa.background import fit_background_quadratic
from dippa.fitting import PatternFitResult, PeakFitResult, fit_one_peak, fit_pattern
from dippa.profiles import evaluate_pattern

FIXTURES = Path(__file__).parent / "fixtures"


def _load_reference():
    fit = scipy.io.loadmat(FIXTURES / "reference_fit.mat")
    data = scipy.io.loadmat(FIXTURES / "reference_data.mat")
    aa = np.asarray(fit["aa"])
    d = np.asarray(data["data"])
    return aa, d[:, 0], d[:, 1]


def _r_squared(y: np.ndarray, model: np.ndarray) -> float:
    ss_res = np.sum((y - model) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1 - ss_res / ss_tot


def _integral_breadth(aa: np.ndarray, peak: int) -> float:
    """Pseudo-Voigt integral breadth, matching the original `getFW_IB.m` ('IB' case).

    ``0.5 * fwhm * (pi*eta + (1-eta)*sqrt(pi/ln 2))`` with fwhm and eta each
    averaged over the two sides for the asymmetric 6-parameter case —
    this is the exact quantity the original tool can feed to
    (modified) Williamson-Hall, so it's the recovery target that matters
    downstream, more than either side's FWHM individually.
    """
    if aa.shape[0] == 4:
        fw, eta = aa[2, peak], aa[3, peak]
    else:
        fw = 0.5 * (aa[2, peak] + aa[4, peak])
        eta = 0.5 * (aa[3, peak] + aa[5, peak])
    return 0.5 * fw * (np.pi * eta + (1 - eta) * np.sqrt(np.pi / np.log(2)))


def _rough_guess(aa_real: np.ndarray, x: np.ndarray, y: np.ndarray, seed: int = 42):
    n_peaks = aa_real.shape[1] - 1
    rng = np.random.default_rng(seed)
    aa_guess = aa_real.copy()
    for i in range(n_peaks):
        aa_guess[0, i] += rng.uniform(-0.005, 0.005)  # position jitter, within the window
        aa_guess[1, i] *= rng.uniform(0.5, 1.5)  # scale the true amplitude by up to +/-50%
        aa_guess[2, i] = 0.002  # generic rough width, not the real value
        aa_guess[3, i] = 0.5  # generic eta
        aa_guess[4, i] = 0.002
        aa_guess[5, i] = 0.5

    c0, c1, c2 = fit_background_quadratic(x, y, aa_guess[0, :n_peaks], half_width=0.02)
    aa_guess[0, -1], aa_guess[1, -1], aa_guess[2, -1] = c0, c1, c2
    return aa_guess


# --- Synthetic recovery test ------------------------------------------------


def test_fit_one_peak_recovers_synthetic_peak():
    x = np.linspace(0.4, 0.6, 4000)
    true_params = np.array([0.5, 3.0, 0.003, 0.5, 0.0025, 0.4])
    aa = np.zeros((6, 2))  # one peak + background
    aa[:, 0] = true_params
    aa[:, 1] = [0.01, -0.005, 0.0, 0.0, 0.0, 0.0]  # small linear-ish background

    rng = np.random.default_rng(1)
    y = evaluate_pattern(x, aa) + rng.normal(scale=1e-3, size=x.shape)

    guess = aa.copy()
    guess[0, 0] += 0.003  # wrong position, within the fitting window
    guess[1, 0] *= 0.6  # wrong amplitude
    guess[2, 0] = 0.001  # wrong width
    guess[4, 0] = 0.001

    fitted, result = fit_one_peak(x, y, guess, peak_index=0, half_width=0.02)

    assert isinstance(result, PeakFitResult)
    assert result.success
    assert fitted[0, 0] == pytest.approx(true_params[0], abs=1e-4)
    assert fitted[1, 0] == pytest.approx(true_params[1], rel=0.05)
    # The result object and the array must agree — one source of truth.
    np.testing.assert_array_equal(result.params, fitted[:, 0])


def test_fit_one_peak_raises_with_no_nearby_data():
    x = np.linspace(0.4, 0.6, 100)
    y = np.ones_like(x)
    aa = np.zeros((4, 2))
    aa[:, 0] = [50.0, 1.0, 0.01, 0.5]  # x0 way outside the data range
    with pytest.raises(ValueError, match="not enough to fit"):
        fit_one_peak(x, y, aa, peak_index=0, half_width=0.02)


def test_fit_one_peak_returns_local_background():
    """The local linear correction must be reported, not fitted-then-discarded.

    Synthetic data with a deliberate linear offset the shared background
    doesn't know about: the fitted local background should pick it up and
    the returned values should reflect it (the aabcg convention).
    """
    x = np.linspace(0.4, 0.6, 4000)
    aa = np.zeros((6, 2))
    aa[:, 0] = [0.5, 3.0, 0.003, 0.5, 0.0025, 0.4]
    aa[:, 1] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # shared background: none

    # Real signal includes a local linear term the shared background omits.
    y = evaluate_pattern(x, aa) + (0.05 - 0.06 * x)

    _, result = fit_one_peak(x, y, aa, peak_index=0, half_width=0.02)
    c0, c1 = result.local_background
    # Over the window the recovered line should reproduce the injected one
    # (c0 and c1 individually are strongly correlated over a narrow window,
    # so assert on the evaluated line, not the coefficients).
    x_mid = 0.5
    assert c0 + c1 * x_mid == pytest.approx(0.05 - 0.06 * x_mid, abs=5e-3)


# --- Real-pattern recovery tests --------------------------------------------


def test_fit_pattern_recovers_real_pattern_from_rough_guess():
    aa_real, x, y = _load_reference()
    n_peaks = aa_real.shape[1] - 1
    aa_guess = _rough_guess(aa_real, x, y)

    starting_r2 = _r_squared(y, evaluate_pattern(x, aa_guess))
    assert starting_r2 < 0.5, "the test fixture should start from a genuinely bad guess"

    result = fit_pattern(x, y, aa_guess, half_width=0.02, n_passes=3)
    assert isinstance(result, PatternFitResult)
    final_r2 = _r_squared(y, evaluate_pattern(x, result.aa))

    assert final_r2 > 0.99, f"expected the fitter to recover R^2 > 0.99, got {final_r2:.4f}"

    # Position recovery should be tight — this is the parameter a real user
    # supplies only approximately (by clicking), so it's the one that
    # matters most for "did the fitter actually find the right peak".
    position_errors = np.abs(result.aa[0, :n_peaks] - aa_real[0, :n_peaks])
    assert np.all(position_errors < 1e-3), f"largest position error: {position_errors.max():.6f}"

    # Amplitude recovery, since the guess deliberately scaled it by +/-50%.
    amp_rel_errors = np.abs(result.aa[1, :n_peaks] / aa_real[1, :n_peaks] - 1)
    assert np.all(amp_rel_errors < 0.10), f"largest amplitude error: {amp_rel_errors.max():.3f}"


def test_fit_pattern_recovers_integral_breadth():
    """Breadth is what (modified) Williamson-Hall consumes, so test it directly.

    FWHM and eta individually are *not* asserted here: on this data they
    trade off against each other (and against the local background) without
    moving the residual much, so per-parameter agreement with the MATLAB
    values is not something this fitter achieves or currently claims — see
    AUDIT.md §15. Integral breadth is the better-identified combination and
    the one that flows into mWH; 15% is its observed recovery envelope on
    this fixture, documented rather than aspirational.
    """
    aa_real, x, y = _load_reference()
    n_peaks = aa_real.shape[1] - 1
    aa_guess = _rough_guess(aa_real, x, y)
    result = fit_pattern(x, y, aa_guess, half_width=0.02, n_passes=3)

    for peak in range(n_peaks):
        ib_fit = _integral_breadth(result.aa, peak)
        ib_real = _integral_breadth(aa_real, peak)
        assert ib_fit == pytest.approx(ib_real, rel=0.15), (
            f"peak {peak}: integral breadth {ib_fit:.5f} vs reference {ib_real:.5f}"
        )


def test_fit_pattern_reports_bound_hits():
    """The known bound-pegging case must be *reported*, not silent.

    On this fixture, at least one peak's eta lands on the 1.3 upper bound
    (observed for the weakest, last peak — its right tail prefers a
    super-Lorentzian in this objective; see AUDIT.md §15). This test pins
    the diagnostic, not the pathology: if a future change makes every fit
    clean, flip this test to assert `result.all_clean` instead — that's an
    improvement, and this test failing in that direction is good news.
    """
    aa_real, x, y = _load_reference()
    aa_guess = _rough_guess(aa_real, x, y)
    result = fit_pattern(x, y, aa_guess, half_width=0.02, n_passes=3)

    assert not result.all_clean
    assert any("upper bound" in w for w in result.warnings)
    # Every result must carry its diagnostics regardless of outcome.
    for r in result.peak_results:
        assert r.n_points > 0
        assert r.nfev > 0
        assert np.isfinite(r.cost)


def test_fit_pattern_returns_aabcg():
    """Local backgrounds persist across passes and come back in aabcg form.

    Shape and convention match the original tool's ``aabcg`` variable
    ((2, n_peaks): offset row then slope row) as read by ``legacy_io``.
    """
    aa_real, x, y = _load_reference()
    n_peaks = aa_real.shape[1] - 1
    aa_guess = _rough_guess(aa_real, x, y)
    result = fit_pattern(x, y, aa_guess, half_width=0.02, n_passes=3)

    assert result.aabcg.shape == (2, n_peaks)
    assert np.all(np.isfinite(result.aabcg))
    # The local corrections should stay small relative to the peaks they sit
    # under — a large value here means the local line is absorbing model
    # mismatch, which is exactly what reporting it is meant to catch.
    local_at_peak = np.abs(result.aabcg[0] + result.aabcg[1] * result.aa[0, :n_peaks])
    assert np.all(local_at_peak < 0.1 * result.aa[1, :n_peaks]), (
        f"local background magnitude vs amplitude: {np.round(local_at_peak / result.aa[1, :n_peaks], 3)}"
    )
