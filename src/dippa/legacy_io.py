"""Read state files produced by the original MATLAB BIGdippa/dippaFC app.

This module is I/O only: it loads the legacy `.mat` fit-state files into
plain numpy arrays with documented shapes, so a saved MATLAB session can be
inspected or migrated without re-deriving the fit. It does not implement or
re-derive any fitting logic. Shapes and field names below were confirmed
against the actual `0. variables/fit.mat` / `fit_I.mat` files bundled with
dippa_v3.zip (github.com/ThomasHSimm/DPPA), not assumed from documentation.

Array conventions (from BIGdippaFunctions/pv_tv_aa.m and onepeak.m):
- `aa` has shape (rows, n_peaks + 1). rows is 4 for the symmetric
  pseudo-Voigt peak model (pk_alpha/pk_voigt2) or 6 for the asymmetric
  variant (pk_alpha_asymm). Columns 0..n_peaks-1 are one peak each; the
  last column holds the shared quadratic background coefficients
  [c0, c1, c2] used as backgd = c0 + c1*x + c2*x**2.
- `aabcg` has shape (2, n_peaks): a local linear background offset/slope
  per peak, separate from the shared background in `aa`'s last column.
- The `_I`-suffixed variants (`aa_I`, `aabcg_I`) hold the same structure
  for the instrumental-standard fit, stored in a separate file
  (`fit_I.mat`) in the original tool.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy.io


@dataclass(frozen=True)
class LegacyFit:
    """A single fitted pattern (sample OR instrument) from the legacy tool."""

    aa: np.ndarray  # shape (rows, n_peaks + 1)
    aabcg: np.ndarray  # shape (2, n_peaks)

    @property
    def n_peaks(self) -> int:
        return self.aa.shape[1] - 1

    @property
    def background_coeffs(self) -> np.ndarray:
        """Shared quadratic background [c0, c1, c2] from aa's last column."""
        return self.aa[:3, -1]

    def peak_params(self, index: int) -> np.ndarray:
        """Parameters for one peak (0-indexed), excluding the background column."""
        if not 0 <= index < self.n_peaks:
            raise IndexError(f"index {index} out of range for {self.n_peaks} peaks")
        return self.aa[:, index]


def load_legacy_fit(path: str | Path, instrument: bool = False) -> LegacyFit:
    """Load a `fit.mat` / `fit_I.mat` file saved by the original MATLAB app.

    Parameters
    ----------
    path:
        Path to the `.mat` file (e.g. "0. variables/fit.mat").
    instrument:
        If True, read the `aa_I` / `aabcg_I` variables (instrumental-standard
        fit) instead of `aa` / `aabcg` (sample fit). The original tool stores
        these in separate files by convention, but this flag is provided in
        case both are ever consolidated into one file.
    """
    data = scipy.io.loadmat(str(path))
    aa_key, aabcg_key = ("aa_I", "aabcg_I") if instrument else ("aa", "aabcg")
    if aa_key not in data or aabcg_key not in data:
        raise KeyError(
            f"expected variables '{aa_key}' and '{aabcg_key}' in {path}, "
            f"found: {sorted(k for k in data if not k.startswith('__'))}"
        )
    return LegacyFit(aa=np.asarray(data[aa_key]), aabcg=np.asarray(data[aabcg_key]))
