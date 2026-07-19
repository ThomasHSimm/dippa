"""Tests for legacy_io. Fixture shapes (6, 11) / (2, 10) match the real
0. variables/fit.mat bundled with the original dippa_v3.zip (10 peaks,
asymmetric 6-parameter profile), not an arbitrary example."""

from pathlib import Path

import numpy as np
import pytest
import scipy.io

from dippa import load_legacy_fit


@pytest.fixture
def fit_mat_path(tmp_path: Path) -> Path:
    rng = np.random.default_rng(0)
    aa = rng.normal(size=(6, 11))
    aabcg = -np.abs(rng.normal(scale=0.005, size=(2, 10)))
    path = tmp_path / "fit.mat"
    scipy.io.savemat(path, {"aa": aa, "aabcg": aabcg})
    return path


def test_load_legacy_fit_shapes(fit_mat_path: Path):
    fit = load_legacy_fit(fit_mat_path)
    assert fit.aa.shape == (6, 11)
    assert fit.aabcg.shape == (2, 10)
    assert fit.n_peaks == 10


def test_background_coeffs_is_last_column(fit_mat_path: Path):
    fit = load_legacy_fit(fit_mat_path)
    np.testing.assert_array_equal(fit.background_coeffs, fit.aa[:3, -1])


def test_peak_params_out_of_range_raises(fit_mat_path: Path):
    fit = load_legacy_fit(fit_mat_path)
    with pytest.raises(IndexError):
        fit.peak_params(10)


def test_missing_variables_raise_keyerror(tmp_path: Path):
    path = tmp_path / "empty.mat"
    scipy.io.savemat(path, {"something_else": np.array([1])})
    with pytest.raises(KeyError):
        load_legacy_fit(path)
