"""Tests for matlab_io module."""

from pathlib import Path

import numpy as np
import pytest

from dippa.matlab_io import (
    MatlabSample,
    extract_struct_array,
    load_matlab_samples,
    squeeze_value,
    struct_to_dict,
)

FIXTURES = Path(__file__).parent / "fixtures"


class TestSqueezeValue:
    def test_scalar_array(self):
        val = np.array([[42.0]])
        assert squeeze_value(val) == 42.0

    def test_singleton_array(self):
        val = np.array([42.0])
        assert squeeze_value(val) == 42.0

    def test_regular_array(self):
        val = np.array([[1, 2, 3], [4, 5, 6]])
        result = squeeze_value(val)
        np.testing.assert_array_equal(result, val)

    def test_scalar_types(self):
        assert squeeze_value(42) == 42
        assert squeeze_value("hello") == "hello"


class TestLoadMatlabSamples:
    def test_load_ni_combo(self):
        # Use minimal fixture (1 sample, decimated data) instead of full ni_combo.mat (2.3 MB)
        samples = load_matlab_samples(FIXTURES / "ni_combo_minimal.mat")

        assert len(samples) == 1
        assert all(isinstance(s, MatlabSample) for s in samples)

    def test_sample_structure(self):
        samples = load_matlab_samples(FIXTURES / "ni_combo_minimal.mat")
        first = samples[0]

        # Pin hardcoded values verified against full ni_combo.mat
        assert first.name == "niHmid_halfpc"
        assert first.n_peaks == 5
        # Minimal fixture decimates data/data_I to every 10th point
        assert first.data.shape == (1161, 2)
        assert first.data_I.shape == (1127, 2)
        assert first.aa.shape == (6, 6)

    def test_background_coeffs(self):
        samples = load_matlab_samples(FIXTURES / "ni_combo_minimal.mat")
        first = samples[0]

        bg = first.background_coeffs
        assert bg.shape == (3,)
        # Pin against full ni_combo.mat sample 0
        np.testing.assert_array_almost_equal(
            bg, [85.92835305, -208.24732128, 139.02646391]
        )

    def test_peak_params(self):
        samples = load_matlab_samples(FIXTURES / "ni_combo_minimal.mat")
        first = samples[0]

        peak_0 = first.peak_params(0)
        assert peak_0.shape == (6,)

    def test_peak_params_out_of_range(self):
        samples = load_matlab_samples(FIXTURES / "ni_combo_minimal.mat")
        first = samples[0]

        with pytest.raises(IndexError):
            first.peak_params(5)

    def test_missing_fita_raises(self):
        with pytest.raises(KeyError, match="expected variable 'fita'"):
            load_matlab_samples(FIXTURES / "reference_fit.mat")


class TestStructToDict:
    def test_unwraps_singleton_dimensions(self):
        import scipy.io

        data = scipy.io.loadmat(str(FIXTURES / "ni_combo_minimal.mat"))
        struct = data["fita"][0][0]

        d = struct_to_dict(struct)

        assert isinstance(d, dict)
        assert "instr" in d
        assert isinstance(d["instr"], str)
        assert isinstance(d["name"], str)


class TestExtractStructArray:
    def test_extracts_all_samples(self):
        import scipy.io

        data = scipy.io.loadmat(str(FIXTURES / "ni_combo_minimal.mat"))
        dicts = extract_struct_array(data["fita"])

        assert len(dicts) == 1  # Minimal fixture has 1 sample
        assert all(isinstance(d, dict) for d in dicts)
        assert all("name" in d for d in dicts)
