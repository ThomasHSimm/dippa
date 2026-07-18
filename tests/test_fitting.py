"""Tests for dippa.fitting.

The important test here is `test_fit_pattern_recovers_real_peaks_from_rough_guess`
— it starts from a starting point that is *worse than predicting the mean*
(R^2 < 0) and checks the staged fitter finds its way back to the real
answer. Evaluating known-good parameters (as profiles.py's tests do) proves
the physics is right; this proves the *fitter* actually works, which is a
different and stronger claim.
"""

from pathlib import Path

import numpy as np
import pytest
import scipy.io

from dippa.background import fit_background_quadratic
from dippa.fitting import fit_one_peak, fit_pattern
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

    fitted = fit_one_peak(x, y, guess, peak_index=0, half_width=0.02)

    assert fitted[0, 0] == pytest.approx(true_params[0], abs=1e-4)
    assert fitted[1, 0] == pytest.approx(true_params[1], rel=0.05)


def test_fit_one_peak_raises_with_no_nearby_data():
    x = np.linspace(0.4, 0.6, 100)
    y = np.ones_like(x)
    aa = np.zeros((4, 2))
    aa[:, 0] = [50.0, 1.0, 0.01, 0.5]  # x0 way outside the data range
    with pytest.raises(ValueError, match="not enough to fit"):
        fit_one_peak(x, y, aa, peak_index=0, half_width=0.02)


# --- Real parity test --------------------------------------------------------


def test_fit_pattern_recovers_real_peaks_from_rough_guess():
    aa_real, x, y = _load_reference()
    n_peaks = aa_real.shape[1] - 1

    rng = np.random.default_rng(42)
    aa_guess = aa_real.copy()
    for i in range(n_peaks):
        aa_guess[0, i] += rng.uniform(-0.005, 0.005)  # position jitter, within the window
        aa_guess[1, i] *= rng.uniform(0.5, 1.5)  # +/-50% amplitude
        aa_guess[2, i] = 0.002  # generic rough width, not the real value
        aa_guess[3, i] = 0.5  # generic eta
        aa_guess[4, i] = 0.002
        aa_guess[5, i] = 0.5

    c0, c1, c2 = fit_background_quadratic(x, y, aa_guess[0, :n_peaks], half_width=0.02)
    aa_guess[0, -1], aa_guess[1, -1], aa_guess[2, -1] = c0, c1, c2

    starting_r2 = _r_squared(y, evaluate_pattern(x, aa_guess))
    assert starting_r2 < 0.5, "the test fixture should start from a genuinely bad guess"

    aa_fitted = fit_pattern(x, y, aa_guess, half_width=0.02, n_passes=3)
    final_r2 = _r_squared(y, evaluate_pattern(x, aa_fitted))

    assert final_r2 > 0.99, f"expected the fitter to recover R^2 > 0.99, got {final_r2:.4f}"

    # Position recovery should be tight — this is the parameter a real user
    # supplies only approximately (by clicking), so it's the one that
    # matters most for "did the fitter actually find the right peak".
    position_errors = np.abs(aa_fitted[0, :n_peaks] - aa_real[0, :n_peaks])
    assert np.all(position_errors < 1e-3), f"largest position error: {position_errors.max():.6f}"
