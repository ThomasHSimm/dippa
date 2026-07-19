"""Input/output helpers for legacy DPPA and MATLAB data files."""

from dippa.io.legacy_io import LegacyFit, load_legacy_fit
from dippa.io.matlab_io import (
    MatlabSample,
    extract_struct_array,
    load_matlab_samples,
    squeeze_value,
    struct_to_dict,
)
from dippa.io.raw import RawScan, TUBE_WAVELENGTHS, load_dat, load_udf, merge_scans

__all__ = [
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
]
