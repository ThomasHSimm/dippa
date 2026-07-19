"""Tests for classical and modified Williamson-Hall fitting."""

import numpy as np
import pytest

from dippa.breadth import DeltaKResult, ExcludedPeak
from dippa.structure import Phase, bind_reflections, generate_reflections, reciprocal_g
from dippa.williamson_hall import (
    WHParameters,
    fit_williamson_hall,
    williamson_hall_model,
)

PHASE = Phase("FCC", 3.52)
HKL = generate_reflections(PHASE, max_index=8, limit=15)
G = np.array([reciprocal_g(PHASE, reflection) for reflection in HKL])
BINDING = bind_reflections(G, PHASE, HKL, tol=1e-12)


def _breadths(delta_k, excluded=()):
    count = len(delta_k)
    return DeltaKResult(
        delta_k=np.asarray(delta_k),
        positions=G[:count],
        sample_peak_indices=np.arange(count, dtype=np.int64),
        instrument_peak_indices=np.arange(count, dtype=np.int64),
        excluded_peaks=excluded,
        policy_messages=(),
        kind="IB",
    )


@pytest.mark.parametrize("variant", ["classical", "mwhA", "mwhB", "mwhC"])
def test_synthetic_round_trip_recovers_parameters_within_ci(variant):
    true = WHParameters(q=0.0 if variant == "classical" else 1.7, size=0.0012, strain=0.004)
    exact = williamson_hall_model(G, BINDING, 0.317, true, variant)
    noise = 2e-5 * np.array([0.2, -0.8, 0.5, 0.9, -0.4] * 3)
    result = fit_williamson_hall(
        _breadths(exact + noise),
        BINDING,
        0.317,
        variant=variant,
        initial=(1.2, 0.001, 0.003),
    )

    assert result.success
    assert result.n_points == 15
    assert result.dof == (13 if variant == "classical" else 12)
    for value, interval in [
        (true.q, result.confidence_intervals.q),
        (true.size, result.confidence_intervals.size),
        (true.strain, result.confidence_intervals.strain),
    ]:
        assert interval[0] <= value <= interval[1]
    assert result.covariance.shape == (3, 3)
    assert np.all(np.isfinite(result.residuals))


def test_classical_is_mwha_with_q_fixed_zero():
    parameters = WHParameters(q=2.5, size=0.001, strain=0.004)
    classical = williamson_hall_model(G, BINDING, 0.317, parameters, "classical")
    mwha = williamson_hall_model(
        G, BINDING, 0.317, WHParameters(q=0.0, size=0.001, strain=0.004), "mwhA"
    )
    np.testing.assert_allclose(classical, mwha)


def test_model_rejects_coordinates_that_do_not_match_binding():
    parameters = WHParameters(q=1.5, size=0.001, strain=0.004)
    with pytest.raises(ValueError, match="does not match"):
        williamson_hall_model(G + 0.01, BINDING, 0.317, parameters, "mwhA")


def test_low_dof_and_exclusions_are_reported():
    parameters = WHParameters(q=1.5, size=0.001, strain=0.004)
    binding = bind_reflections(G[:5], PHASE, HKL[:5], tol=1e-12)
    delta = williamson_hall_model(G[:5], binding, 0.317, parameters, "mwhA")
    excluded = (ExcludedPeak(5, ("nonpositive delta_k",)),)
    result = fit_williamson_hall(_breadths(delta, excluded), BINDING, 0.317, variant="mwhA")

    assert result.dof == 2
    assert any("low degrees of freedom" in warning for warning in result.warnings)
    assert result.excluded_peaks == excluded


def test_too_few_retained_peaks_raises():
    with pytest.raises(ValueError, match="at least 3"):
        fit_williamson_hall(_breadths([0.001, 0.002]), BINDING, 0.317, variant="mwhA")
