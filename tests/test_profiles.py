"""Tests for dippa.profiles.

Two kinds of test here, deliberately kept separate:

1. Synthetic self-consistency tests — no external data, check the profile
   functions behave sanely (peak at x0, symmetric case matches the
   asymmetric case when both sides are equal, doublet adds a second
   component, etc).
2. A real parity test against ``fit.mat``/``data.mat`` bundled from the
   original MATLAB tool's own example dataset (MIT licensed,
   github.com/ThomasHSimm/DPPA) — this is the actual "does the port
   reproduce the original tool" check, not just an internal consistency
   check. See AUDIT.md for how this fixture was chosen and validated.
"""

from pathlib import Path

import numpy as np
import pytest
import scipy.io

from dippa.profiles import (
    TUBE_WAVELENGTHS,
    asymmetric_pseudo_voigt,
    evaluate_pattern,
    evaluate_peak,
    pseudo_voigt,
)

FIXTURES = Path(__file__).parent / "fixtures"


# --- Synthetic self-consistency tests --------------------------------------


def test_pseudo_voigt_peaks_at_x0():
    x = np.linspace(-1, 1, 2001)
    y = pseudo_voigt(x, x0=0.0, amplitude=5.0, fwhm=0.2, eta=0.5)
    assert y.argmax() == 1000  # x=0 is the centre index
    assert y[1000] == pytest.approx(5.0)


def test_pseudo_voigt_pure_lorentzian_and_gaussian_bounds():
    x = np.linspace(-1, 1, 501)
    lorentzian = pseudo_voigt(x, 0.0, 1.0, 0.2, eta=1.0)
    gaussian = pseudo_voigt(x, 0.0, 1.0, 0.2, eta=0.0)
    # Lorentzian has heavier tails than Gaussian for the same FWHM.
    assert lorentzian[0] > gaussian[0]


def test_asymmetric_matches_symmetric_when_both_sides_equal():
    x = np.linspace(-1, 1, 1001)
    sym = pseudo_voigt(x, 0.0, 3.0, 0.15, 0.4)
    asym = asymmetric_pseudo_voigt(
        x, x0=0.0, amplitude=3.0, fwhm_right=0.15, eta_left=0.4, fwhm_left=0.15, eta_right=0.4
    )
    np.testing.assert_allclose(sym, asym, atol=1e-12)


def test_asymmetric_sides_differ_when_parameters_differ():
    x = np.linspace(-1, 1, 1001)
    asym = asymmetric_pseudo_voigt(
        x, x0=0.0, amplitude=1.0, fwhm_right=0.3, eta_left=0.5, fwhm_left=0.05, eta_right=0.5
    )
    # Narrower left FWHM => faster falloff just left of centre than just right.
    i_left = np.searchsorted(x, -0.1)
    i_right = np.searchsorted(x, 0.1)
    assert asym[i_left] < asym[i_right]


def test_doublet_adds_a_second_smaller_component():
    x = np.linspace(0.4, 0.6, 5001)
    params = np.array([0.5, 4.0, 0.002, 0.5, 0.002, 0.5])
    single = evaluate_peak(x, params, tube=None)
    doublet = evaluate_peak(x, params, tube="Co")
    # The doublet's alpha2 component sits to the right of x0 and adds area,
    # so the doublet curve must be >= the single curve everywhere and
    # strictly greater somewhere to the right of x0.
    assert np.all(doublet >= single - 1e-12)
    assert doublet[x > 0.502].sum() > single[x > 0.502].sum()


def test_tube_wavelengths_match_original_hardcoded_values():
    # Pinned against pk_alpha.m / pk_alpha_asymm.m — these are exact literals
    # from the original source, not computed, so a straight equality check
    # is the right test (any drift here silently changes doublet spacing).
    assert TUBE_WAVELENGTHS["Cu"] == (1.54056, 1.54439)
    assert TUBE_WAVELENGTHS["Co"] == (1.78897, 1.79285)
    assert TUBE_WAVELENGTHS["Fe"] == (1.93604, 1.93998)


def test_unsupported_tube_raises():
    x = np.linspace(0, 1, 10)
    params = np.array([0.5, 1.0, 0.01, 0.5])
    with pytest.raises(ValueError, match="unsupported tube"):
        evaluate_peak(x, params, tube="Mo")


def test_wrong_parameter_count_raises():
    x = np.linspace(0, 1, 10)
    with pytest.raises(ValueError, match="4 .* or 6"):
        evaluate_peak(x, np.array([0.5, 1.0, 0.01]), tube=None)


def test_evaluate_pattern_background_only_with_zero_peaks():
    x = np.linspace(0, 1, 100)
    aa = np.array([[0.1], [0.2], [-0.05], [0.0], [0.0], [0.0]])  # 0 peaks, background only
    y = evaluate_pattern(x, aa)
    expected = 0.1 + 0.2 * x - 0.05 * x**2
    np.testing.assert_allclose(y, expected)


# --- Real parity test against the original MATLAB tool ---------------------


def _load_reference():
    fit = scipy.io.loadmat(FIXTURES / "reference_fit.mat")
    data = scipy.io.loadmat(FIXTURES / "reference_data.mat")
    aa = np.asarray(fit["aa"])
    d = np.asarray(data["data"])
    return aa, d[:, 0], d[:, 1]


def test_parity_against_original_tool_no_doublet():
    """Reproduce the original tool's own saved fit on its own bundled data.

    This is the key finding from the audit: the bundled genset.mat for this
    example has alpha2=0 (which in the original tool's inverted convention
    means "fit the doublet"), but the saved fit only reproduces the real
    pattern when NO doublet is applied. See the profiles.py module
    docstring and AUDIT.md — this test pins that finding down so it can't
    silently regress if someone "fixes" it back to using the doublet.
    """
    aa, x, y = _load_reference()
    model = evaluate_pattern(x, aa, tube=None)

    ss_res = np.sum((y - model) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r_squared = 1 - ss_res / ss_tot

    assert r_squared > 0.99, f"expected R^2 > 0.99 against the real reference pattern, got {r_squared:.4f}"


def test_parity_with_doublet_is_much_worse():
    """Guard rail: confirm the doublet really is wrong for this fixture,
    not just "also fine". If this test starts failing, the no-doublet
    result above may need re-checking too."""
    aa, x, y = _load_reference()
    model_with_doublet = evaluate_pattern(x, aa, tube="Co")

    ss_res = np.sum((y - model_with_doublet) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r_squared = 1 - ss_res / ss_tot

    assert r_squared < 0.8, f"expected the doublet-applied fit to be clearly worse, got R^2={r_squared:.4f}"
