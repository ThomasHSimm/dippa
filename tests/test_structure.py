"""Tests for cubic/HCP phases, notation, assignment and lattice refinement."""

import numpy as np
import pytest

from dippa.structure import (
    Phase,
    Reflection,
    assign_reflections,
    bind_reflections,
    format_hkil,
    generate_reflections,
    reciprocal_g,
    refine_lattice_parameter,
    reflection_allowed,
)


def test_phase_mirrors_dsettings_lat1_lat2():
    titanium = Phase(cstruct="HCP", a=2.9506, c=4.6835)
    assert titanium.lat1 == 2.9506
    assert titanium.lat2 == 4.6835
    assert Phase(cstruct="fcc", a=3.52).cstruct == "FCC"
    with pytest.raises(ValueError, match="require.*c"):
        Phase(cstruct="HCP", a=2.9506)


def test_miller_bravais_input_is_validated_and_stored_three_index():
    phase = Phase(cstruct="HCP", a=2.9506, c=4.6835)
    assert phase.reflection((1, 0, -1, 2)) == Reflection(1, 0, 2)
    with pytest.raises(ValueError, match=r"i must equal -\(h\+k\)"):
        phase.reflection((1, 0, 0, 2))


def test_hcp_dot_notation_and_cubic_notation():
    assert format_hkil((1, 0, 0), "HCP") == "10.0"
    assert format_hkil((0, 0, 2), "HCP") == "00.2"
    assert format_hkil((1, 0, -1, 1), "HCP") == "10.1"
    assert format_hkil((1, 1, 1), "FCC") == "(111)"


def test_hcp_forbidden_reflections():
    assert not reflection_allowed((0, 0, 1), "HCP")
    assert not reflection_allowed((1, 1, 1), "HCP")
    assert reflection_allowed((1, 0, 1), "HCP")
    assert reflection_allowed((0, 0, 2), "HCP")


def test_titanium_reflections_have_canonical_order_and_absences():
    phase = Phase(cstruct="HCP", a=2.9506, c=4.6835)
    reflections = generate_reflections(phase, max_index=4, limit=10)
    labels = [format_hkil(reflection, "HCP") for reflection in reflections]
    assert labels[:6] == ["10.0", "00.2", "10.1", "10.2", "11.0", "10.3"]
    assert "00.1" not in labels
    assert "11.1" not in labels


def test_hcp_lattice_refinement_round_trip():
    phase = Phase(cstruct="HCP", a=2.9506, c=4.6835)
    reflections = generate_reflections(phase, max_index=4, limit=10)
    positions = np.array([reciprocal_g(phase, reflection) for reflection in reflections])
    result = refine_lattice_parameter(positions, reflections, "HCP")

    assert result.a == pytest.approx(phase.a, rel=1e-12)
    assert result.c == pytest.approx(phase.c, rel=1e-12)
    assert result.c_over_a == pytest.approx(phase.c / phase.a, rel=1e-12)
    assert result.c_over_a == pytest.approx(1.587, abs=5e-4)
    np.testing.assert_allclose(result.residuals, 0.0, atol=1e-15)


@pytest.mark.parametrize(
    ("phase", "expected"),
    [
        (Phase("FCC", 3.52), ["(111)", "(200)", "(220)", "(311)", "(222)"]),
        (Phase("BCC", 2.87), ["(110)", "(200)", "(211)", "(220)", "(310)"]),
    ],
)
def test_cubic_reflection_generation_and_refinement(phase, expected):
    reflections = generate_reflections(phase, max_index=5, limit=5)
    assert [format_hkil(reflection, phase.cstruct) for reflection in reflections] == expected
    positions = np.array([reciprocal_g(phase, reflection) for reflection in reflections])
    result = refine_lattice_parameter(positions, reflections, phase.cstruct)
    assert result.a == pytest.approx(phase.a, rel=1e-12)
    assert result.c is None
    assert result.c_over_a is None


@pytest.mark.parametrize(
    "phase",
    [Phase("FCC", 3.52), Phase("HCP", 2.9506, 4.6835)],
)
def test_assign_reflections_works_for_cubic_and_hcp(phase):
    expected = generate_reflections(phase, max_index=5, limit=6)
    exact = np.array([reciprocal_g(phase, reflection) for reflection in expected])
    observed = exact + np.linspace(-2e-5, 2e-5, len(exact))
    binding = assign_reflections(observed, phase, tol=3e-5, max_index=5)
    assert list(binding.reflections) == list(expected)
    np.testing.assert_allclose(
        [assignment.residual for assignment in binding.assignments], observed - exact
    )


def test_fcc_nickel_assignment_and_lattice_refinement_match_ni_combo():
    from dippa import load_matlab_samples

    sample = load_matlab_samples("tests/fixtures/ni_combo_minimal.mat")[0]
    positions = sample.aa[0, :-1]
    phase = Phase("FCC", 3.53)
    binding = assign_reflections(positions, phase, tol=2.5e-4, max_index=5)

    assert [format_hkil(item, "FCC") for item in binding.reflections] == [
        "(111)",
        "(200)",
        "(220)",
        "(311)",
        "(222)",
    ]
    result = refine_lattice_parameter(positions, binding.reflections, "FCC")
    assert result.a == pytest.approx(3.52982, abs=5e-5)
    assert 3.52 < result.a < 3.54


def test_assignment_errors_when_tolerance_makes_match_ambiguous():
    phase = Phase("FCC", 3.52)
    first_two = generate_reflections(phase, max_index=2, limit=2)
    positions = np.array([reciprocal_g(phase, item) for item in first_two])
    midpoint = np.array([positions.mean()])

    with pytest.raises(ValueError, match="ambiguous"):
        assign_reflections(midpoint, phase, tol=0.05, max_index=2)


def test_explicit_binding_rejects_duplicate_hkl():
    phase = Phase("FCC", 3.52)
    position = reciprocal_g(phase, (1, 1, 1))
    with pytest.raises(ValueError, match="unique"):
        bind_reflections([position, position], phase, [(1, 1, 1), (1, 1, 1)], tol=1e-6)


def test_binding_subset_is_reindexed_and_remains_validated():
    phase = Phase("FCC", 3.52)
    reflections = generate_reflections(phase, max_index=5, limit=5)
    positions = np.array([reciprocal_g(phase, item) for item in reflections])
    binding = bind_reflections(positions, phase, reflections, tol=1e-12)

    selected = binding.select([0, 2, 4])

    assert [item.peak_index for item in selected.assignments] == [0, 1, 2]
    assert selected.reflections == (reflections[0], reflections[2], reflections[4])
