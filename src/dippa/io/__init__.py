"""Input/output helpers for legacy DPPA and MATLAB data files."""

from dippa.io.legacy_io import LegacyFit, load_legacy_fit
from dippa.io.matlab_io import (
    MatlabSample,
    extract_struct_array,
    load_matlab_samples,
    squeeze_value,
    struct_to_dict,
)

__all__ = [
    "LegacyFit",
    "MatlabSample",
    "extract_struct_array",
    "load_legacy_fit",
    "load_matlab_samples",
    "squeeze_value",
    "struct_to_dict",
]
