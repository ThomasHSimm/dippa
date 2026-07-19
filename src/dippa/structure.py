"""Crystal phases, reflection assignment and lattice refinement.

The ``Phase`` schema mirrors the original MATLAB ``dsettings`` convention:
``lat1`` is the basal/cubic lattice parameter a and ``lat2`` is c for HCP.
The HCP metric ports ``BIGdippaFunctions/dspacing_dippa.m``. Reflections are
stored internally as three-index (h, k, l), matching ``dsettings.index``;
four-index Miller-Bravais input is validated at the API boundary.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Sequence

import numpy as np
from numpy.typing import NDArray

CrystalStructure = Literal["FCC", "BCC", "HCP"]
FloatArray = NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class Reflection:
    """One reflection stored in the original tool's three-index form."""

    h: int
    k: int
    l: int  # noqa: E741 - conventional Miller index name

    def as_tuple(self) -> tuple[int, int, int]:
        return self.h, self.k, self.l


ReflectionLike = Reflection | Sequence[int]


def _reflection(value: ReflectionLike) -> Reflection:
    if isinstance(value, Reflection):
        return value
    indices = tuple(int(index) for index in value)
    if len(indices) == 3:
        return Reflection(*indices)
    if len(indices) == 4:
        h, k, i, ell = indices
        expected = -(h + k)
        if i != expected:
            raise ValueError(
                f"invalid Miller-Bravais reflection {indices}: i must equal -(h+k)={expected}"
            )
        return Reflection(h, k, ell)
    raise ValueError(f"reflection must have 3 or 4 indices, got {len(indices)}")


@dataclass(frozen=True, slots=True)
class Phase:
    """Cubic or HCP phase using the original ``dsettings`` lattice fields."""

    cstruct: CrystalStructure
    a: float
    c: float | None = None

    def __post_init__(self) -> None:
        structure = self.cstruct.upper()
        if structure not in {"FCC", "BCC", "HCP"}:
            raise ValueError(f"unsupported crystal structure: {self.cstruct!r}")
        object.__setattr__(self, "cstruct", structure)
        if not np.isfinite(self.a) or self.a <= 0:
            raise ValueError(f"a must be finite and positive, got {self.a}")
        if structure == "HCP":
            if self.c is None or not np.isfinite(self.c) or self.c <= 0:
                raise ValueError("HCP phases require finite positive c")
        elif self.c is not None:
            raise ValueError(f"{structure} phases do not use c/lat2")

    @property
    def lat1(self) -> float:
        return self.a

    @property
    def lat2(self) -> float | None:
        return self.c

    def reflection(self, value: ReflectionLike) -> Reflection:
        """Validate notation and return the internal three-index form."""
        return _reflection(value)


@dataclass(frozen=True, slots=True)
class ReflectionAssignment:
    peak_index: int
    reflection: Reflection
    observed_g: float
    expected_g: float
    residual: float


@dataclass(frozen=True, slots=True)
class ReflectionBinding:
    """Validated one-to-one binding between fitted peaks and reflections."""

    phase: Phase
    assignments: tuple[ReflectionAssignment, ...]
    tolerance: float

    def __post_init__(self) -> None:
        if not np.isfinite(self.tolerance) or self.tolerance <= 0:
            raise ValueError("binding tolerance must be finite and positive")
        if not self.assignments:
            raise ValueError("a reflection binding cannot be empty")
        peak_indices = tuple(item.peak_index for item in self.assignments)
        if peak_indices != tuple(range(len(self.assignments))):
            raise ValueError("binding peak indices must be contiguous and ordered")
        reflections = tuple(item.reflection for item in self.assignments)
        if len(set(reflections)) != len(reflections):
            raise ValueError("each reflection may be bound to only one peak")
        if any(not reflection_allowed(item.reflection, self.phase.cstruct) for item in self.assignments):
            raise ValueError("binding contains a forbidden reflection")
        if any(abs(item.residual) > self.tolerance for item in self.assignments):
            raise ValueError("binding residual exceeds its validation tolerance")
        if any(
            not np.isfinite((item.observed_g, item.expected_g, item.residual)).all()
            for item in self.assignments
        ):
            raise ValueError("binding coordinates and residuals must be finite")

    @property
    def positions(self) -> FloatArray:
        return np.array([item.observed_g for item in self.assignments])

    @property
    def reflections(self) -> tuple[Reflection, ...]:
        return tuple(item.reflection for item in self.assignments)

    def select(self, peak_indices: Sequence[int]) -> ReflectionBinding:
        """Return a reindexed binding for a retained subset of peaks."""
        indices = tuple(int(index) for index in peak_indices)
        if len(set(indices)) != len(indices):
            raise ValueError("selected peak indices must be unique")
        if any(index < 0 or index >= len(self.assignments) for index in indices):
            raise ValueError("selected peak index is outside the binding")
        selected = tuple(self.assignments[index] for index in indices)
        return ReflectionBinding(
            phase=self.phase,
            assignments=tuple(
                ReflectionAssignment(
                    peak_index=new_index,
                    reflection=item.reflection,
                    observed_g=item.observed_g,
                    expected_g=item.expected_g,
                    residual=item.residual,
                )
                for new_index, item in enumerate(selected)
            ),
            tolerance=self.tolerance,
        )


@dataclass(frozen=True, slots=True)
class LatticeRefinementResult:
    cstruct: CrystalStructure
    a: float
    c: float | None
    c_over_a: float | None
    fitted_g: FloatArray
    residuals: FloatArray
    n_points: int


def format_hkil(hkl: ReflectionLike, cstruct: CrystalStructure = "HCP") -> str:
    """Format HCP with implicit i as ``hk.l`` and cubic as ``(hkl)``."""
    reflection = _reflection(hkl)
    structure = cstruct.upper()
    if structure == "HCP":
        return f"{reflection.h}{reflection.k}.{reflection.l}"
    if structure in {"FCC", "BCC"}:
        return f"({reflection.h}{reflection.k}{reflection.l})"
    raise ValueError(f"unsupported crystal structure: {cstruct!r}")


def reflection_allowed(hkl: ReflectionLike, cstruct: CrystalStructure) -> bool:
    """Return the FCC, BCC or HCP systematic-absence selection rule."""
    reflection = _reflection(hkl)
    h, k, ell = reflection.as_tuple()
    if h == k == ell == 0:
        return False
    structure = cstruct.upper()
    if structure == "FCC":
        parity = (abs(h) % 2, abs(k) % 2, abs(ell) % 2)
        return parity == (0, 0, 0) or parity == (1, 1, 1)
    if structure == "BCC":
        return (h + k + ell) % 2 == 0
    if structure == "HCP":
        return not ((h + 2 * k) % 3 == 0 and ell % 2 != 0)
    raise ValueError(f"unsupported crystal structure: {cstruct!r}")


def reciprocal_g(phase: Phase, hkl: ReflectionLike) -> float:
    """Return ``g = 1/d`` for a cubic or HCP reflection."""
    reflection = phase.reflection(hkl)
    h, k, ell = reflection.as_tuple()
    if not reflection_allowed(reflection, phase.cstruct):
        raise ValueError(f"forbidden {phase.cstruct} reflection: {format_hkil(reflection, phase.cstruct)}")
    if phase.cstruct in {"FCC", "BCC"}:
        g_squared = (h * h + k * k + ell * ell) / phase.a**2
    else:
        assert phase.c is not None
        g_squared = (4.0 / 3.0) * (h * h + h * k + k * k) / phase.a**2
        g_squared += ell * ell / phase.c**2
    return float(np.sqrt(g_squared))


def d_spacing(phase: Phase, hkl: ReflectionLike) -> float:
    """Return d for one allowed reflection."""
    return 1.0 / reciprocal_g(phase, hkl)


def generate_reflections(
    phase: Phase, *, max_index: int = 8, limit: int | None = None
) -> tuple[Reflection, ...]:
    """Generate unique positive-index families in increasing reciprocal g."""
    if max_index < 1:
        raise ValueError("max_index must be at least 1")
    candidates: list[Reflection] = []
    if phase.cstruct == "HCP":
        for h in range(max_index + 1):
            for k in range(h + 1):
                for ell in range(max_index + 1):
                    reflection = Reflection(h, k, ell)
                    if reflection_allowed(reflection, phase.cstruct):
                        candidates.append(reflection)
    else:
        for h in range(max_index + 1):
            for k in range(h + 1):
                for ell in range(k + 1):
                    reflection = Reflection(h, k, ell)
                    if reflection_allowed(reflection, phase.cstruct):
                        candidates.append(reflection)
    ordered = sorted(candidates, key=lambda reflection: reciprocal_g(phase, reflection))
    return tuple(ordered[:limit] if limit is not None else ordered)


def assign_reflections(
    positions: FloatArray,
    phase: Phase,
    tol: float,
    *,
    max_index: int = 8,
) -> ReflectionBinding:
    """Bind peaks only when exactly one unused candidate lies within ``tol``."""
    observed = np.asarray(positions, dtype=np.float64)
    if observed.ndim != 1 or not observed.size or not np.all(np.isfinite(observed)):
        raise ValueError("positions must be a nonempty finite one-dimensional array")
    if not np.isfinite(tol) or tol <= 0:
        raise ValueError("tol must be finite and positive")
    candidates = generate_reflections(phase, max_index=max_index)
    if len(candidates) < len(observed):
        raise ValueError("not enough candidate reflections for the observed peaks")
    expected = np.array([reciprocal_g(phase, reflection) for reflection in candidates])
    available = set(range(len(candidates)))
    assignments: list[ReflectionAssignment | None] = [None] * len(observed)
    for peak_index in np.argsort(observed):
        matches = [
            index for index in available if abs(expected[index] - observed[peak_index]) <= tol
        ]
        if not matches:
            raise ValueError(f"peak {peak_index} has no reflection match within tol={tol}")
        if len(matches) > 1:
            labels = ", ".join(format_hkil(candidates[index], phase.cstruct) for index in matches)
            raise ValueError(f"peak {peak_index} has ambiguous reflection matches: {labels}")
        candidate_index = matches[0]
        available.remove(candidate_index)
        assignments[int(peak_index)] = ReflectionAssignment(
            peak_index=int(peak_index),
            reflection=candidates[candidate_index],
            observed_g=float(observed[peak_index]),
            expected_g=float(expected[candidate_index]),
            residual=float(observed[peak_index] - expected[candidate_index]),
        )
    return ReflectionBinding(
        phase=phase,
        assignments=tuple(assignment for assignment in assignments if assignment is not None),
        tolerance=tol,
    )


def bind_reflections(
    positions: FloatArray,
    phase: Phase,
    hkl: Sequence[ReflectionLike],
    tol: float,
) -> ReflectionBinding:
    """Validate an explicit peak↔HKL assignment without auto-indexing."""
    observed = np.asarray(positions, dtype=np.float64)
    reflections = tuple(phase.reflection(value) for value in hkl)
    if observed.ndim != 1 or len(observed) != len(reflections):
        raise ValueError("positions and hkl must be one-dimensional with matching lengths")
    if len(set(reflections)) != len(reflections):
        raise ValueError("explicit hkl assignments must be unique")
    if not np.isfinite(tol) or tol <= 0:
        raise ValueError("tol must be finite and positive")
    assignments = tuple(
        ReflectionAssignment(
            peak_index=index,
            reflection=reflection,
            observed_g=float(position),
            expected_g=reciprocal_g(phase, reflection),
            residual=float(position - reciprocal_g(phase, reflection)),
        )
        for index, (position, reflection) in enumerate(zip(observed, reflections))
    )
    return ReflectionBinding(phase=phase, assignments=assignments, tolerance=tol)


def refine_lattice_parameter(
    positions: FloatArray,
    hkl: Sequence[ReflectionLike],
    cstruct: CrystalStructure = "FCC",
) -> LatticeRefinementResult:
    """Refine cubic a or HCP a/c from assigned reciprocal peak positions."""
    observed = np.asarray(positions, dtype=np.float64)
    reflections = tuple(_reflection(value) for value in hkl)
    if observed.ndim != 1 or len(observed) != len(reflections):
        raise ValueError("positions and hkl must be one-dimensional with matching lengths")
    if not len(observed) or not np.all(np.isfinite(observed)) or np.any(observed <= 0):
        raise ValueError("positions must be finite and positive")
    structure = cstruct.upper()
    if any(not reflection_allowed(reflection, structure) for reflection in reflections):
        raise ValueError("hkl contains a forbidden reflection")
    indices = np.asarray([reflection.as_tuple() for reflection in reflections], dtype=float)
    h, k, ell = indices.T
    target = observed**2
    if structure == "HCP":
        design = np.column_stack((h * h + h * k + k * k, ell * ell))
        if np.linalg.matrix_rank(design) < 2:
            raise ValueError("HCP refinement needs independent basal and c-axis reflections")
        coefficient_a, coefficient_c = np.linalg.lstsq(design, target, rcond=None)[0]
        if coefficient_a <= 0 or coefficient_c <= 0:
            raise ValueError("refined HCP reciprocal-metric coefficients must be positive")
        a = float(np.sqrt(4.0 / (3.0 * coefficient_a)))
        c = float(1.0 / np.sqrt(coefficient_c))
        fitted_squared = design @ np.array([coefficient_a, coefficient_c])
        ratio: float | None = c / a
    elif structure in {"FCC", "BCC"}:
        metric = h * h + k * k + ell * ell
        slope = float(np.linalg.lstsq(metric[:, np.newaxis], target, rcond=None)[0][0])
        if slope <= 0:
            raise ValueError("refined cubic reciprocal-metric coefficient must be positive")
        a = float(1.0 / np.sqrt(slope))
        c = None
        ratio = None
        fitted_squared = metric * slope
    else:
        raise ValueError(f"unsupported crystal structure: {cstruct!r}")
    fitted = np.sqrt(fitted_squared)
    return LatticeRefinementResult(
        cstruct=structure,
        a=a,
        c=c,
        c_over_a=ratio,
        fitted_g=fitted,
        residuals=observed - fitted,
        n_points=len(observed),
    )


refine_lattice = refine_lattice_parameter
