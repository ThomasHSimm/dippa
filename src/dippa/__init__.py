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
from dippa.contrast import Reflection, contrast_cubic, h_squared
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

__version__ = "0.0.1.dev0"

__all__ = [
    "DeltaKResult",
    "ExcludedPeak",
    "WHConfidenceIntervals",
    "WHParameters",
    "WilliamsonHallResult",
    "Reflection",
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
    "delta_k",
    "h_squared",
    "fit_williamson_hall",
    "peak_breadths",
    "williamson_hall_model",
    "regrid_to_g",
    "theta_to_g",
    "g_to_theta",
]
