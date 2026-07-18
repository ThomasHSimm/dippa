"""Tests for dippa.background."""

from pathlib import Path

import numpy as np
import pytest
import scipy.io

from dippa.background import background_mask, fit_background_quadratic

FIXTURES = Path(__file__).parent / "fixtures"


def test_background_mask_excludes_peak_windows():
    x = np.linspace(0, 10, 1001)
    mask = background_mask(x, peak_positions=np.array([3.0, 7.0]), half_width=1.0)
    assert not mask[np.argmin(np.abs(x - 3.0))]
    assert not mask[np.argmin(np.abs(x - 7.0))]
    assert mask[np.argmin(np.abs(x - 0.0))]
    assert mask[np.argmin(np.abs(x - 5.0))]


def test_fit_background_quadratic_recovers_known_coefficients():
    rng = np.random.default_rng(0)
    x = np.linspace(0, 10, 2000)
    true_c0, true_c1, true_c2 = 1.0, -0.3, 0.05
    y_background = true_c0 + true_c1 * x + true_c2 * x**2
    # add a fake peak that should be excluded, and small noise everywhere
    y = y_background + 5.0 * np.exp(-((x - 5.0) ** 2) / (2 * 0.1**2))
    y += rng.normal(scale=1e-4, size=x.shape)

    c0, c1, c2 = fit_background_quadratic(x, y, peak_positions=np.array([5.0]), half_width=0.5)

    assert c0 == pytest.approx(true_c0, abs=1e-2)
    assert c1 == pytest.approx(true_c1, abs=1e-2)
    assert c2 == pytest.approx(true_c2, abs=1e-2)


def test_fit_background_quadratic_raises_with_too_few_points():
    x = np.linspace(0, 1, 10)
    y = np.ones_like(x)
    # half_width covering the whole range leaves ~nothing outside the peak window
    with pytest.raises(ValueError, match="background points remain"):
        fit_background_quadratic(x, y, peak_positions=np.array([0.5]), half_width=0.6)


def test_background_close_to_real_reference_fit():
    """Not a strict parity test (the original's point-selection window isn't
    known exactly) — just a sanity check that the closed-form fit lands in
    the right neighbourhood of the real saved background coefficients."""
    fit = scipy.io.loadmat(FIXTURES / "reference_fit.mat")
    data = scipy.io.loadmat(FIXTURES / "reference_data.mat")
    aa = np.asarray(fit["aa"])
    d = np.asarray(data["data"])
    x, y = d[:, 0], d[:, 1]

    c0, c1, c2 = fit_background_quadratic(x, y, peak_positions=aa[0, :10], half_width=0.02)

    assert c0 == pytest.approx(aa[0, -1], abs=0.01)
    assert c1 == pytest.approx(aa[1, -1], abs=0.03)
    assert c2 == pytest.approx(aa[2, -1], abs=0.02)
