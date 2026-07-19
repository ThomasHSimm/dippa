"""Tests for breadth definitions and instrumental subtraction policies."""

from types import SimpleNamespace

import numpy as np
import pytest

from dippa.breadth import delta_k, peak_breadths


def _aa(positions, fwhm, eta=0.5):
    positions = np.asarray(positions, dtype=float)
    fwhm = np.broadcast_to(fwhm, positions.shape)
    eta = np.broadcast_to(eta, positions.shape)
    aa = np.zeros((6, len(positions) + 1))
    aa[0, :-1] = positions
    aa[1, :-1] = 1.0
    aa[2, :-1] = fwhm * 0.8
    aa[4, :-1] = fwhm * 1.2
    aa[3, :-1] = eta * 0.5
    aa[5, :-1] = eta * 1.5
    return aa


def _fit(aa, flagged=()):
    results = []
    for peak_index in range(aa.shape[1] - 1):
        results.append(
            SimpleNamespace(
                peak_index=peak_index,
                success=True,
                at_lower_bound=(),
                at_upper_bound=("eta_left",) if peak_index in flagged else (),
            )
        )
    return SimpleNamespace(aa=aa, peak_results=tuple(results))


def test_peak_breadths_side_average_before_ib():
    aa = _aa([0.5], fwhm=0.004, eta=0.6)
    expected = 0.5 * 0.004 * (
        np.pi * 0.6 + (1.0 - 0.6) * np.sqrt(np.pi / np.log(2.0))
    )
    np.testing.assert_allclose(peak_breadths(aa, "FW"), [0.004])
    np.testing.assert_allclose(peak_breadths(aa), [expected])


def test_delta_k_matches_instrument_by_nearest_position():
    sample = _aa([0.50, 0.80], [0.005, 0.008])
    instrument = _aa([0.79, 0.49], [0.002, 0.001])
    result = delta_k(sample, instrument, kind="FW")

    np.testing.assert_allclose(result.delta_k, [0.004, 0.006])
    np.testing.assert_array_equal(result.instrument_peak_indices, [1, 0])
    assert result.excluded_peaks == ()


def test_flagged_peaks_excluded_by_default_and_overridable():
    sample = _fit(_aa([0.5, 0.8], [0.005, 0.008]), flagged=(1,))
    instrument = _aa([0.5, 0.8], [0.001, 0.002])

    safe = delta_k(sample, instrument, kind="FW")
    np.testing.assert_array_equal(safe.sample_peak_indices, [0])
    assert safe.excluded_peaks[0].peak_index == 1
    assert "upper bound" in safe.excluded_peaks[0].reasons[0]
    assert safe.policy_messages

    overridden = delta_k(sample, instrument, kind="FW", include_flagged=True)
    np.testing.assert_array_equal(overridden.sample_peak_indices, [0, 1])
    assert overridden.excluded_peaks == ()


def test_nonpositive_delta_k_raises_or_warns_and_drops():
    sample = _aa([0.5, 0.8], [0.001, 0.008])
    instrument = _aa([0.5, 0.8], [0.002, 0.002])

    with pytest.raises(ValueError, match=r"peak\(s\): 0"):
        delta_k(sample, instrument, kind="FW")

    with pytest.warns(UserWarning, match=r"peak\(s\): 0"):
        dropped = delta_k(sample, instrument, kind="FW", drop_nonpositive=True)
    np.testing.assert_array_equal(dropped.sample_peak_indices, [1])
    assert dropped.excluded_peaks[0].reasons == ("nonpositive delta_k",)
    assert "nonpositive" in dropped.policy_messages[-1]
