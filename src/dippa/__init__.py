"""dippa: diffraction peak profile analysis in Python.

Pre-alpha. See AUDIT.md in the repository root for the state of the port
from the original MATLAB tool (github.com/ThomasHSimm/DPPA).
"""

from dippa.io import (
    LegacyFit,
    MatlabSample,
    RawScan,
    TUBE_WAVELENGTHS,
    extract_struct_array,
    load_legacy_fit,
    load_dat,
    load_matlab_samples,
    load_udf,
    merge_scans,
    squeeze_value,
    struct_to_dict,
)
from dippa.contrast import contrast_cubic, h_squared
from dippa.structure import (
    LatticeRefinementResult,
    Phase,
    Reflection,
    ReflectionAssignment,
    ReflectionBinding,
    assign_reflections,
    bind_reflections,
    d_spacing,
    format_hkil,
    generate_reflections,
    reciprocal_g,
    refine_lattice,
    refine_lattice_parameter,
    reflection_allowed,
)
from dippa.breadth import DeltaKResult, ExcludedPeak, delta_k, peak_breadths
from dippa.williamson_hall import (
    WHConfidenceIntervals,
    WHParameters,
    WilliamsonHallResult,
    fit_williamson_hall,
    williamson_hall_model,
)
from dippa.regrid import (
    g_to_theta,
    regrid_to_g,
    theta_to_g,
)
from dippa.analysis import Analysis, FitConfig, Pattern
from dippa.profiles import PeakParams

__version__ = "0.0.1.dev0"

__all__ = [
    "DeltaKResult",
    "Analysis",
    "FitConfig",
    "ExcludedPeak",
    "LatticeRefinementResult",
    "Phase",
    "Pattern",
    "PeakParams",
    "WHConfidenceIntervals",
    "WHParameters",
    "WilliamsonHallResult",
    "Reflection",
    "ReflectionAssignment",
    "ReflectionBinding",
    "LegacyFit",
    "MatlabSample",
    "RawScan",
    "TUBE_WAVELENGTHS",
    "extract_struct_array",
    "load_legacy_fit",
    "load_dat",
    "load_matlab_samples",
    "load_udf",
    "merge_scans",
    "squeeze_value",
    "struct_to_dict",
    "contrast_cubic",
    "assign_reflections",
    "bind_reflections",
    "d_spacing",
    "delta_k",
    "h_squared",
    "fit_williamson_hall",
    "format_hkil",
    "generate_reflections",
    "peak_breadths",
    "reciprocal_g",
    "refine_lattice",
    "refine_lattice_parameter",
    "reflection_allowed",
    "williamson_hall_model",
    "regrid_to_g",
    "theta_to_g",
    "g_to_theta",
]
