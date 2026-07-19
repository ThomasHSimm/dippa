"""Cubic diffraction contrast factors.

``h_squared`` ports ``FourierFunctions/Hsq.m`` from the original DPPA tool.
The original obtains ``ch00`` from the GUI preference ``WHpref.chk0``;
neither that value nor ``q`` is inferred here, so both remain explicit.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatArray = NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class Reflection:
    """One explicitly assigned crystallographic reflection."""

    h: int
    k: int
    l: int  # noqa: E741 - conventional Miller index name

    def as_tuple(self) -> tuple[int, int, int]:
        return self.h, self.k, self.l


def _hkl_array(hkl: ArrayLike | Sequence[Reflection]) -> FloatArray:
    values = [item.as_tuple() if isinstance(item, Reflection) else item for item in hkl]
    array = np.asarray(values, dtype=np.float64)
    if array.ndim == 1:
        array = array[np.newaxis, :]
    if array.ndim != 2 or array.shape[1] != 3:
        raise ValueError(f"hkl must have shape (n, 3), got {array.shape}")
    if np.any(np.sum(array**2, axis=1) == 0):
        raise ValueError("the (0, 0, 0) reflection is undefined")
    return array


def h_squared(hkl: ArrayLike | Sequence[Reflection]) -> FloatArray:
    """Return vectorised cubic H² values, ported from ``Hsq.m``.

    ``H² = (h²k² + h²l² + k²l²) / (h² + k² + l²)²``.
    Reflections must be supplied explicitly; this module never assigns HKLs
    from peak order or position.
    """
    squared = _hkl_array(hkl) ** 2
    h2, k2, l2 = squared.T
    return (h2 * k2 + h2 * l2 + k2 * l2) / (h2 + k2 + l2) ** 2


def contrast_cubic(
    hkl: ArrayLike | Sequence[Reflection], ch00: float, q: float
) -> FloatArray:
    """Return ``C = ch00 * (1 - q * H²)`` for cubic reflections."""
    if not np.isfinite(ch00) or ch00 <= 0:
        raise ValueError(f"ch00 must be finite and positive, got {ch00}")
    if not np.isfinite(q):
        raise ValueError(f"q must be finite, got {q}")
    return ch00 * (1.0 - q * h_squared(hkl))
