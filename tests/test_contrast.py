"""Tests for cubic reflection and contrast-factor calculations."""

import numpy as np
import pytest

from dippa.contrast import contrast_cubic, h_squared
from dippa.structure import Phase, Reflection, bind_reflections, reciprocal_g

FCC_HKL = [
    Reflection(1, 1, 1),
    Reflection(2, 0, 0),
    Reflection(2, 2, 0),
    Reflection(3, 1, 1),
    Reflection(2, 2, 2),
]


def _fcc_binding():
    phase = Phase("FCC", 1.0)
    positions = np.array([reciprocal_g(phase, reflection) for reflection in FCC_HKL])
    return bind_reflections(positions, phase, FCC_HKL, tol=1e-12)


def test_h_squared_vectorised_fcc_reflections():
    np.testing.assert_allclose(h_squared(FCC_HKL), [1 / 3, 0, 1 / 4, 19 / 121, 1 / 3])


def test_h_squared_rejects_invalid_reflections():
    with pytest.raises(ValueError, match="shape"):
        h_squared([[1, 1]])
    with pytest.raises(ValueError, match="undefined"):
        h_squared([[0, 0, 0]])


def test_contrast_cubic_matches_ss316_fcres_set_zero():
    # Exact oracle from FCres_set(1).FCres.C in the original DPPA file
    # 4. RES/Fourier/SS316_logINDI_RES.mat. That saved WA-track result stores
    # chk0=0.317 and q=2.7131932196981099 for these five FCC reflections.
    stored = np.array(
        [
            0.030305916451899746,
            0.317,
            0.1019794373389248,
            0.18194576229552301,
            0.030305916451899746,
        ]
    )
    actual = contrast_cubic(_fcc_binding(), ch00=0.317, q=2.7131932196981099)
    np.testing.assert_allclose(actual, stored, rtol=0.0, atol=5e-15)


def test_contrast_cubic_requires_explicit_valid_constants():
    with pytest.raises(ValueError, match="ch00"):
        contrast_cubic(_fcc_binding(), ch00=0.0, q=2.0)
    with pytest.raises(ValueError, match="q"):
        contrast_cubic(_fcc_binding(), ch00=0.317, q=np.nan)
