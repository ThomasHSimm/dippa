"""dippa: diffraction peak profile analysis in Python.

Pre-alpha. See AUDIT.md in the repository root for the state of the port
from the original MATLAB tool (github.com/ThomasHSimm/DPPA).
"""

from dippa.io import (
    LegacyFit,
    MatlabSample,
    extract_struct_array,
    load_legacy_fit,
    load_matlab_samples,
    squeeze_value,
    struct_to_dict,
)
from dippa.regrid import (
    g_to_theta,
    regrid_to_g,
    theta_to_g,
)

__version__ = "0.0.1.dev0"

__all__ = [
    "LegacyFit",
    "MatlabSample",
    "extract_struct_array",
    "load_legacy_fit",
    "load_matlab_samples",
    "squeeze_value",
    "struct_to_dict",
    "regrid_to_g",
    "theta_to_g",
    "g_to_theta",
]
