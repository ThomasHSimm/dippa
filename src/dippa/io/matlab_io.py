"""Utilities for extracting data from MATLAB `.mat` files.

MATLAB structs loaded via scipy.io.loadmat have awkward structure:
- Struct arrays are nested in extra dimensions (shape (1,) or (1, N))
- Field values are wrapped in singleton dimensions
- Accessing nested fields requires brittle indexing like struct[0][0]

This module provides clean extraction functions that unwrap these
abstractions and return properly shaped arrays or dicts.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, TypeVar

import numpy as np
import scipy.io

T = TypeVar("T")


def squeeze_value(val: Any) -> Any:
    """Unwrap a MATLAB field value to remove singleton dimensions.

    MATLAB scalar fields come as shape (1,) or (1, 1) arrays; this
    returns the scalar. Arrays are unchanged.
    """
    if not isinstance(val, np.ndarray):
        return val
    if val.size == 1:
        return val.item()
    return val


def struct_to_dict(struct: np.void) -> dict[str, Any]:
    """Convert a numpy structured array element to a plain dict.

    Parameters
    ----------
    struct : np.void
        A single element from a structured array (e.g. fita[0][0]).

    Returns
    -------
    dict
        Fields unwrapped to remove singleton dimensions where present.
    """
    result = {}
    for field in struct.dtype.names:
        val = struct[field]
        result[field] = squeeze_value(val)
    return result


def extract_struct_array(
    arr: np.ndarray, flat: bool = True
) -> list[dict[str, Any]]:
    """Extract a MATLAB struct array to a list of dicts.

    Parameters
    ----------
    arr : np.ndarray
        A struct array from scipy.io.loadmat (shape is typically (1, N) or
        (1,) for single element; each element is a np.void with named fields).
    flat : bool
        If True (default), squeeze the array to 1-D before extracting
        (removes the (1, N) shape to (N,)).

    Returns
    -------
    list[dict]
        One dict per struct element, with fields unwrapped.
    """
    if flat:
        arr = arr.flat
    return [struct_to_dict(elem) for elem in arr]


@dataclass(frozen=True)
class MatlabSample:
    """A single measurement from a MATLAB combo fit file.

    Extracted from scipy.io.loadmat('.mat')['fita'] elements. Wraps the
    original tool's output struct; not all fields are populated in all files.

    **Field semantics**:

    - `name` (str): Sample identifier
    - `val`, `valstrain`, `valsstress` (all float): **Same value in all three fields**.
      No separate strain/stress tracking; all three equal.
    - `instr` (str): Tube type (e.g., 'Co' for Cobalt Kα radiation).
    - `alpha2` (int): **Inverted MATLAB convention**: 0 = Kα1/Kα2 doublet fitted (ON);
      1 = single wavelength (OFF).
    - `data` (ndarray, shape (n_meas, 2)): Measured diffraction pattern.
      Columns: [g coordinate, intensity (counts)].
    - `data_I` (ndarray, shape (n_meas, 2)): Instrumental standard pattern
      (same coordinate grid).
    - `aa`, `aabcg` (ndarrays): Fitted peak parameters and background
      (sample fit).
    - `aa_I`, `aabcg_I` (ndarrays): Fitted peak parameters and background
      (instrumental standard fit).
    - `IB` (ndarray, shape (1, 11)): Unidentified field. Treated as opaque payload.
    """

    instr: str
    alpha2: int
    data_I: np.ndarray
    aa_I: np.ndarray
    aa: np.ndarray
    data: np.ndarray
    aabcg_I: np.ndarray
    aabcg: np.ndarray
    IB: np.ndarray
    name: str
    val: float
    FileN: str
    PathN: str
    valstrain: float
    valsstress: float

    @property
    def n_peaks(self) -> int:
        """Number of peaks in the fit (aa has n_peaks + 1 columns for background)."""
        return self.aa.shape[1] - 1

    @property
    def background_coeffs(self) -> np.ndarray:
        """Shared quadratic background [c0, c1, c2] from aa's last column."""
        return self.aa[:3, -1]

    def peak_params(self, index: int) -> np.ndarray:
        """Parameters for one peak (0-indexed)."""
        if not 0 <= index < self.n_peaks:
            raise IndexError(f"index {index} out of range for {self.n_peaks} peaks")
        return self.aa[:, index]


def load_matlab_samples(path: str | Path) -> list[MatlabSample]:
    """Load a MATLAB combo file (e.g. ni_combo.mat) into structured samples.

    Parameters
    ----------
    path : str | Path
        Path to the `.mat` file.

    Returns
    -------
    list[MatlabSample]
        One sample per measurement in the file.

    Raises
    ------
    KeyError
        If the file does not contain the expected 'fita' variable.
    ValueError
        If a sample is missing required fields.
    """
    data = scipy.io.loadmat(str(path))
    if "fita" not in data:
        raise KeyError(
            f"expected variable 'fita' in {path}, found: "
            f"{sorted(k for k in data if not k.startswith('__'))}"
        )

    samples = []
    fita = data["fita"]
    for elem in fita.flat:
        d = struct_to_dict(elem)
        try:
            sample = MatlabSample(
                instr=str(d["instr"]),
                alpha2=int(d["alpha2"]),
                data_I=d["data_I"],
                aa_I=d["aa_I"],
                aa=d["aa"],
                data=d["data"],
                aabcg_I=d["aabcg_I"],
                aabcg=d["aabcg"],
                IB=d["IB"],
                name=str(d["name"]),
                val=float(d["val"]),
                FileN=str(d["FileN"]),
                PathN=str(d["PathN"]),
                valstrain=float(d["valstrain"]),
                valsstress=float(d["valsstress"]),
            )
            samples.append(sample)
        except KeyError as e:
            raise ValueError(f"sample missing required field: {e}") from e

    return samples
