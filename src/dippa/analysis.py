"""Object-oriented convenience shell around dippa's functional fitting core."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import NDArray

from dippa.fitting import PatternFitResult, fit_pattern
from dippa.io.matlab_io import MatlabSample
from dippa.profiles import TUBE_WAVELENGTHS, evaluate_peak
from dippa.structure import Phase, ReflectionBinding, assign_reflections

FloatArray = NDArray[np.float64]


def _validate_tube(tube: str | None) -> None:
    if tube is not None and tube not in TUBE_WAVELENGTHS:
        raise ValueError(f"unsupported tube {tube!r}; supported: {sorted(TUBE_WAVELENGTHS)}")


@dataclass(frozen=True, slots=True)
class FitConfig:
    """Settings passed unchanged to :func:`dippa.fitting.fit_pattern`."""

    half_width: float
    n_passes: int = 2
    tube: str | None = None

    def __post_init__(self) -> None:
        if not np.isfinite(self.half_width) or self.half_width <= 0:
            raise ValueError("half_width must be finite and positive")
        if isinstance(self.n_passes, bool) or not isinstance(self.n_passes, int) or self.n_passes < 1:
            raise ValueError("n_passes must be a positive integer")
        _validate_tube(self.tube)


@dataclass(frozen=True, slots=True)
class Pattern:
    """Measured intensity on dippa's reciprocal-space coordinate."""

    x: FloatArray
    I: FloatArray  # noqa: E741 - intensity follows the established diffraction notation
    coordinate: Literal["g_A_inv"] = "g_A_inv"
    tube: str | None = None

    def __post_init__(self) -> None:
        x = np.asarray(self.x, dtype=np.float64)
        intensity = np.asarray(self.I, dtype=np.float64)
        if x.ndim != 1 or intensity.ndim != 1 or x.shape != intensity.shape or not x.size:
            raise ValueError("x and I must be nonempty one-dimensional arrays with matching shapes")
        if not np.all(np.isfinite(x)) or not np.all(np.isfinite(intensity)):
            raise ValueError("x and I must contain only finite values")
        if self.coordinate != "g_A_inv":
            raise ValueError("coordinate must be 'g_A_inv'")
        _validate_tube(self.tube)
        object.__setattr__(self, "x", x.copy())
        object.__setattr__(self, "I", intensity.copy())


@dataclass(slots=True)
class Analysis:
    """Compose data and crystallography while delegating numerical work to the core API."""

    pattern: Pattern
    phase: Phase
    binding: ReflectionBinding
    config: FitConfig
    aa: FloatArray
    results: PatternFitResult | None = None

    def __post_init__(self) -> None:
        aa = np.asarray(self.aa, dtype=np.float64)
        if aa.ndim != 2 or aa.shape[0] not in (4, 6) or aa.shape[1] < 2:
            raise ValueError("aa must have 4 or 6 rows and at least one peak plus background")
        if len(self.binding.assignments) != aa.shape[1] - 1:
            raise ValueError("binding must contain one assignment per peak")
        if self.binding.phase != self.phase:
            raise ValueError("binding phase must match analysis phase")
        if self.pattern.tube != self.config.tube:
            raise ValueError("pattern.tube and config.tube must match")
        self.aa = aa.copy()

    def fit(self) -> PatternFitResult:
        """Fit through the existing functional core and retain its structured result."""
        self.results = fit_pattern(
            self.pattern.x,
            self.pattern.I,
            self.aa,
            half_width=self.config.half_width,
            tube=self.config.tube,
            n_passes=self.config.n_passes,
        )
        return self.results

    def plot_peak(self, i: int):
        """Plot one measured peak and shade its configured fitting window."""
        if not 0 <= i < self.aa.shape[1] - 1:
            raise IndexError(f"peak index {i} is out of range")
        import matplotlib.pyplot as plt

        aa = self.results.aa if self.results is not None else self.aa
        x0 = float(aa[0, i])
        mask = np.abs(self.pattern.x - x0) <= self.config.half_width
        fig, ax = plt.subplots()
        ax.plot(self.pattern.x[mask], self.pattern.I[mask], ".", label="measured")
        ax.plot(
            self.pattern.x[mask],
            evaluate_peak(self.pattern.x[mask], aa[:, i], tube=self.config.tube),
            label="peak model",
        )
        ax.axvspan(
            x0 - self.config.half_width,
            x0 + self.config.half_width,
            alpha=0.15,
            color="tab:blue",
            label="fit window",
        )
        ax.set(xlabel="g (Å⁻¹)", ylabel="intensity", title=f"Peak {i}")
        ax.legend()
        return fig, ax

    @classmethod
    def from_matlab_sample(
        cls,
        sample: MatlabSample,
        phase: Phase,
        *,
        binding: ReflectionBinding | None = None,
        config: FitConfig | None = None,
        binding_tolerance: float = 0.02,
        tube: str | None = None,
    ) -> Analysis:
        """Construct from a loaded sample; ``tube`` is explicit, never inferred from alpha2."""
        pattern = Pattern(sample.data[:, 0], sample.data[:, 1], tube=tube)
        config = config or FitConfig(half_width=0.02, tube=tube)
        if binding is None:
            binding = assign_reflections(sample.aa[0, :-1], phase, tol=binding_tolerance)
        return cls(pattern, phase, binding, config, sample.aa)
