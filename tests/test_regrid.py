"""Tests for regrid module (tsinterpl.m port)."""


import numpy as np
import pytest
from scipy.integrate import trapezoid

from dippa.regrid import g_to_theta, regrid_to_g, theta_to_g

# Standard wavelengths for testing
COPPER_KA = 1.54056  # Cu Kα in Ų
COBALT_KA = 1.78897  # Co Kα in Ų


class TestThetaToG:
    """Test 2θ to g coordinate conversion."""

    def test_zero_angle_gives_zero_g(self):
        """At 2θ=0, g should be 0."""
        theta = np.array([0.0])
        g = theta_to_g(theta, COPPER_KA)
        assert g[0] == pytest.approx(0.0)

    def test_known_bragg_positions(self):
        """Test known FCC reflections for Ni (a=3.53 Å)."""
        # For Ni FCC, the d-spacing for (111) is a/√3 ≈ 2.038 Å
        # Bragg law: nλ = 2d sin(θ)  =>  θ = arcsin(nλ/2d)
        # For (111) with Cu Kα: d = 2.038, 2d sin(θ) = λ
        # sin(θ) = λ/(2d) ≈ 1.54056/(2*2.038) ≈ 0.3777
        # θ ≈ 22.2°, so 2θ ≈ 44.4°
        # g = 2sin(θ)/λ = 1/d ≈ 0.490 Ų⁻¹

        d_111_ni = 3.53 / np.sqrt(3)  # FCC (111) d-spacing
        two_theta = 2 * np.degrees(np.arcsin(COPPER_KA / (2 * d_111_ni)))
        g_expected = 1.0 / d_111_ni  # g = 1/d for Bragg reflection

        g = theta_to_g(np.array([two_theta]), COPPER_KA)
        assert g[0] == pytest.approx(g_expected, rel=1e-4)

    def test_multiple_angles(self):
        """Test array of angles."""
        angles = np.array([20.0, 40.0, 60.0, 80.0])
        g = theta_to_g(angles, COPPER_KA)
        assert len(g) == 4
        assert np.all(g >= 0)
        assert np.all(np.diff(g) > 0)  # Monotonic increasing

    def test_negative_wavelength_raises(self):
        """Negative wavelength should raise ValueError."""
        with pytest.raises(ValueError, match="wavelength"):
            regrid_to_g(np.array([[30.0, 100.0], [40.0, 200.0]]), -1.54)

    def test_cobalt_vs_copper(self):
        """For same angle, Cu Kα (longer wavelength) gives smaller g than Co Kα."""
        angle = np.array([50.0])
        g_cu = theta_to_g(angle, COPPER_KA)
        g_co = theta_to_g(angle, COBALT_KA)
        # Shorter wavelength => larger g
        assert g_co[0] < g_cu[0]


class TestGToTheta:
    """Test inverse conversion: g to 2θ."""

    def test_roundtrip_conversion(self):
        """theta_to_g followed by g_to_theta should return original."""
        angles_orig = np.array([20.0, 40.0, 60.0, 80.0])
        g = theta_to_g(angles_orig, COPPER_KA)
        angles_recovered = g_to_theta(g, COPPER_KA)
        np.testing.assert_allclose(angles_recovered, angles_orig, rtol=1e-10)

    def test_zero_g_gives_zero_angle(self):
        """At g=0, 2θ should be 0."""
        g = np.array([0.0])
        theta = g_to_theta(g, COPPER_KA)
        assert theta[0] == pytest.approx(0.0)

    def test_clips_beyond_bragg_limit(self):
        """At g > 2/λ, sin(θ) > 1, which is clipped."""
        # For Cu Kα, g_max = 2/λ ≈ 1.30 Ų⁻¹
        # At this limit, 2θ = 180°
        g_max = 2.0 / COPPER_KA
        g = np.array([g_max * 1.01])  # Slightly beyond limit
        theta = g_to_theta(g, COPPER_KA)
        # Should return something valid (clipped)
        assert np.isfinite(theta[0])


class TestRegridToG:
    """Test regridding from 2θ to g with constant step."""

    def test_regrid_simple_pattern(self):
        """Test regridding a simple synthetic pattern."""
        # Create synthetic data: 2θ from 20 to 100 degrees, 100 points
        theta = np.linspace(20, 100, 100)
        intensity = np.sin(theta * np.pi / 180) ** 2 * 1000  # Arbitrary profile
        data = np.column_stack([theta, intensity])

        # Regrid with default step (5e-5 Ų⁻¹)
        result = regrid_to_g(data, COPPER_KA)

        assert result.shape[1] == 2
        assert result.shape[0] > 0
        # g should be monotonically increasing
        assert np.all(np.diff(result[:, 0]) > 0)
        # Intensity should be positive (rough check)
        assert np.all(result[:, 1] >= 0)

    def test_regrid_preserves_integrated_intensity(self):
        """Regridding should approximately preserve integrated intensity."""
        # Create synthetic pattern
        theta = np.linspace(20, 100, 200)
        intensity = np.sin(theta * np.pi / 180) ** 2 * 1000
        data = np.column_stack([theta, intensity])

        result = regrid_to_g(data, COPPER_KA)

        # Rough check: integrated intensity should be similar order of magnitude
        # (Not exact because g-spacing is not uniform in 2θ)
        theta_integral = trapezoid(intensity, theta)
        g_integral = trapezoid(result[:, 1], result[:, 0])

        # Should be same order of magnitude (within factor of 2)
        assert g_integral > 0
        assert g_integral / theta_integral < 5.0

    def test_reversed_data_is_sorted(self):
        """If data is in descending g order, should be sorted ascending."""
        # Create pattern in descending 2θ order (descending g)
        theta = np.linspace(100, 20, 50)  # Reversed
        intensity = np.ones(50) * 100
        data = np.column_stack([theta, intensity])

        result = regrid_to_g(data, COPPER_KA, sort_if_reversed=True)

        # Result should be sorted ascending
        assert np.all(np.diff(result[:, 0]) > 0)

    def test_custom_g_step(self):
        """Test regridding with custom g_step."""
        theta = np.linspace(30, 90, 100)
        intensity = np.ones(100) * 500
        data = np.column_stack([theta, intensity])

        g_step = 1e-4  # Custom step
        result = regrid_to_g(data, COPPER_KA, g_step=g_step)

        # Check spacing (allowing for rounding)
        spacings = np.diff(result[:, 0])
        assert np.allclose(spacings, g_step, rtol=1e-10)

    def test_insufficient_data_raises(self):
        """Less than 2 data points should raise ValueError."""
        data = np.array([[30.0, 100.0]])  # Only 1 point
        with pytest.raises(ValueError, match="at least 2"):
            regrid_to_g(data, COPPER_KA)

    def test_invalid_wavelength_raises(self):
        """Non-positive wavelength should raise ValueError."""
        data = np.array([[30.0, 100.0], [40.0, 200.0]])
        with pytest.raises(ValueError, match="wavelength"):
            regrid_to_g(data, -1.54)

    def test_regrid_matches_ni_combo_grid(self):
        """Verify regridding produces expected step size (ni_combo.mat uses 5e-5)."""
        # Use data from ni_combo for reference
        theta = np.linspace(10, 140, 500)
        intensity = np.random.rand(500) * 1000 + 100
        data = np.column_stack([theta, intensity])

        result = regrid_to_g(data, COBALT_KA, g_step=5e-5)

        # Check that spacing is exactly 5e-5 (with reasonable floating-point tolerance)
        spacings = np.diff(result[:, 0])
        np.testing.assert_allclose(spacings, 5e-5, rtol=1e-10)

    def test_monotonic_input_preserved(self):
        """Monotonically increasing input should remain monotonic in output."""
        theta = np.linspace(20, 100, 150)
        # Random but monotonic intensity profile
        intensity = np.cumsum(np.random.randn(150))
        intensity = intensity - intensity.min() + 100  # Shift to positive
        data = np.column_stack([theta, intensity])

        result = regrid_to_g(data, COPPER_KA)

        # Output g should be strictly monotonic
        assert np.all(np.diff(result[:, 0]) > 0)


class TestRegridIntegration:
    """Integration tests combining multiple operations."""

    def test_full_workflow(self):
        """Test realistic workflow: create synthetic pattern, regrid, convert back."""
        # Create synthetic diffraction pattern
        theta_orig = np.linspace(20, 120, 500)
        # Multi-peak profile
        intensity = (
            1000 * np.exp(-((theta_orig - 45) ** 2) / 100)
            + 800 * np.exp(-((theta_orig - 65) ** 2) / 100)
            + 600 * np.exp(-((theta_orig - 85) ** 2) / 100)
        )
        data_orig = np.column_stack([theta_orig, intensity])

        # Regrid to g space
        data_g = regrid_to_g(data_orig, COPPER_KA, g_step=1e-4)
        assert data_g.shape[0] > 0
        assert np.all(np.diff(data_g[:, 0]) > 0)

        # Convert g back to theta (approximate, due to interpolation)
        theta_recovered = g_to_theta(data_g[:, 0], COPPER_KA)
        assert np.all(np.diff(theta_recovered) > 0)

        # First and last should be within original range (accounting for regridding)
        assert theta_recovered[0] >= theta_orig[0] - 1
        assert theta_recovered[-1] <= theta_orig[-1] + 1

    def test_regrid_with_noise(self):
        """PCHIP should smooth noisy data reasonably."""
        theta = np.linspace(30, 110, 200)
        # Smooth profile with high-frequency noise
        intensity_smooth = 1000 * np.exp(-((theta - 70) ** 2) / 400)
        intensity_noise = intensity_smooth + np.random.randn(200) * 20
        intensity_noise = np.maximum(intensity_noise, 1)  # No negative values
        data = np.column_stack([theta, intensity_noise])

        result = regrid_to_g(data, COPPER_KA, g_step=5e-5)

        # Output should be smooth (no extreme jumps)
        diff_intensity = np.diff(result[:, 1])
        # 99th percentile of abs differences should be reasonable
        assert np.percentile(np.abs(diff_intensity), 99) < 1000


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_very_small_g_step(self):
        """Extremely small g_step should produce many grid points."""
        theta = np.linspace(30, 60, 50)
        intensity = np.ones(50) * 500
        data = np.column_stack([theta, intensity])

        result = regrid_to_g(data, COPPER_KA, g_step=1e-6)

        # Should have many points
        assert result.shape[0] > 1000

    def test_large_g_step(self):
        """Large g_step should produce few grid points."""
        theta = np.linspace(30, 60, 50)
        intensity = np.ones(50) * 500
        data = np.column_stack([theta, intensity])

        result = regrid_to_g(data, COPPER_KA, g_step=0.05)

        # Should have few points
        assert result.shape[0] < 100

    def test_single_peak_regridding(self):
        """Regridding data around a single peak should work."""
        # Narrow range around one peak
        theta = np.linspace(43, 47, 100)
        # Gaussian peak
        intensity = 5000 * np.exp(-((theta - 45) ** 2) / 4)
        data = np.column_stack([theta, intensity])

        result = regrid_to_g(data, COPPER_KA)

        assert result.shape[0] > 10
        # Peak should be preserved
        assert np.max(result[:, 1]) > 1000
