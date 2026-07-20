"""Tests for the object-oriented shell; numerical fitting remains tested in test_fitting.py."""

import matplotlib
import numpy as np
import pytest

from dippa.analysis import Analysis, FitConfig, Pattern
from dippa.fitting import PatternFitResult
from dippa.profiles import PeakParams, evaluate_pattern
from dippa.structure import Phase, bind_reflections, reciprocal_g

matplotlib.use("Agg")


def _analysis() -> Analysis:
    phase = Phase("FCC", a=3.6)
    hkl = [(1, 1, 1)]
    x0 = reciprocal_g(phase, hkl[0])
    x = np.linspace(x0 - 0.04, x0 + 0.04, 301)
    aa = np.array(
        [
            [x0, 0.2],
            [5.0, 0.1],
            [0.006, 0.0],
            [0.4, 0.0],
            [0.006, 0.0],
            [0.4, 0.0],
        ]
    )
    pattern = Pattern(x, evaluate_pattern(x, aa))
    binding = bind_reflections(np.array([x0]), phase, hkl, tol=1e-12)
    return Analysis(pattern, phase, binding, FitConfig(half_width=0.02, n_passes=1), aa)


def test_peak_params_round_trip_and_evaluate_like_array():
    array = np.array([0.5, 4.0, 0.01, 0.2, 0.02, 0.8])
    params = PeakParams.from_array(array)
    np.testing.assert_array_equal(params.to_array(), array)
    x = np.linspace(0.4, 0.6, 101)
    from dippa.profiles import evaluate_peak

    np.testing.assert_allclose(evaluate_peak(x, params), evaluate_peak(x, array))


@pytest.mark.parametrize(
    "array, message",
    [
        ([0.5, 1.0, 0.01], "expected 6"),
        ([0.5, -1.0, 0.01, 0.2, 0.01, 0.2], "amplitude"),
        ([0.5, 1.0, 0.0, 0.2, 0.01, 0.2], "fwhm"),
        ([0.5, 1.0, 0.01, 1.4, 0.01, 0.2], "eta"),
    ],
)
def test_peak_params_validation(array, message):
    with pytest.raises(ValueError, match=message):
        PeakParams.from_array(np.array(array))


def test_pattern_and_fit_config_validation():
    with pytest.raises(ValueError, match="matching shapes"):
        Pattern(np.arange(2.0), np.arange(3.0))
    with pytest.raises(ValueError, match="coordinate"):
        Pattern(np.arange(2.0), np.arange(2.0), coordinate="2theta")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="positive integer"):
        FitConfig(half_width=0.02, n_passes=0)


def test_analysis_fit_delegates_and_retains_result():
    analysis = _analysis()
    result = analysis.fit()
    assert isinstance(result, PatternFitResult)
    assert analysis.results is result
    assert result.n_passes == analysis.config.n_passes


def test_plot_peak_shades_exact_fit_window():
    analysis = _analysis()
    _, ax = analysis.plot_peak(0)
    display_vertices = ax.patches[0].get_transform().transform(ax.patches[0].get_path().vertices)
    vertices = ax.transData.inverted().transform(display_vertices)
    assert vertices[:, 0].min() == pytest.approx(analysis.aa[0, 0] - analysis.config.half_width)
    assert vertices[:, 0].max() == pytest.approx(analysis.aa[0, 0] + analysis.config.half_width)


def test_from_matlab_sample_uses_explicit_tube():
    from pathlib import Path

    from dippa.io.matlab_io import load_matlab_samples

    sample = load_matlab_samples(Path(__file__).parent / "fixtures" / "ni_combo_minimal.mat")[0]
    phase = Phase("FCC", a=3.52)
    binding = bind_reflections(
        sample.aa[0, :-1],
        phase,
        [(1, 1, 1), (2, 0, 0), (2, 2, 0), (3, 1, 1), (2, 2, 2)],
        tol=0.05,
    )
    analysis = Analysis.from_matlab_sample(sample, phase, binding=binding, tube=None)
    assert analysis.pattern.tube is None
    assert analysis.config.tube is None
    np.testing.assert_array_equal(analysis.aa, sample.aa)
