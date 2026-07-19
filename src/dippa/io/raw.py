"""Read raw diffraction scans produced by the original DPPA workflow.

The UDF angle reconstruction follows ``udfcon.pl`` from
github.com/ThomasHSimm/DPPA: each count is assigned
``start_angle + index * step_size``. Tube wavelengths follow
``BIGdippaFunctions/loadb_GUI.m``; in particular, Co is loaded using
``0.5 * (1.78897 + 1.79285) = 1.79091`` Angstrom.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Sequence

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]

TUBE_WAVELENGTHS: dict[str, float] = {
    "Cu": 0.5 * (1.54056 + 1.54439),
    "Co": 0.5 * (1.78897 + 1.79285),
}


@dataclass(frozen=True)
class RawScan:
    """A raw two-theta diffraction scan and its available metadata."""

    data: FloatArray
    sample_id: str | None = None
    tube: str | None = None
    wavelength: float | None = None


def load_dat(path: str | Path) -> FloatArray:
    """Load a two-column angle/count ``.dat`` file."""
    data = np.loadtxt(path, dtype=np.float64)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"expected at least two columns in {path}")
    data = data[:, :2]
    if not np.all(np.diff(data[:, 0]) > 0):
        raise ValueError(f"angles must be strictly increasing in {path}")
    return data


def load_udf(path: str | Path) -> RawScan:
    """Load a Philips/PANalytical UDF scan, including tube and sample ID.

    Counts are converted to angle/count pairs with the same formula as the
    original DPPA ``udfcon.pl`` converter. The header's wavelength fields are
    retained only as source metadata: for supported tubes, ``wavelength`` uses
    the load-time constants selected by ``loadb_GUI.m``.
    """
    header: dict[str, list[str]] = {}
    counts: list[float] = []
    in_scan = False

    with Path(path).open(encoding="ascii") as udf:
        for raw_line in udf:
            line = raw_line.strip()
            if not line:
                continue
            if line == "RawScan":
                in_scan = True
                continue
            if in_scan:
                values = (value.strip() for value in line.rstrip("/").split(","))
                counts.extend(float(value) for value in values if value)
                continue
            fields = [field.strip() for field in line.rstrip("/").split(",")]
            header[fields[0]] = fields[1:]

    try:
        start, end = (float(value) for value in header["DataAngleRange"][:2])
        step = float(header["ScanStepSize"][0])
    except (KeyError, IndexError, ValueError) as error:
        raise ValueError(f"missing or invalid UDF scan range in {path}") from error
    if not counts:
        raise ValueError(f"no RawScan counts found in {path}")

    angles = start + np.arange(len(counts), dtype=np.float64) * step
    if not np.isclose(angles[-1], end, rtol=0.0, atol=step * 1e-6):
        raise ValueError(
            f"UDF count range ends at {angles[-1]:g}, expected {end:g}; "
            "the number of counts does not match the detector positions"
        )

    tube = header.get("Anode", [None])[0]
    sample_id = header.get("SampleIdent", [None])[0]
    wavelength = TUBE_WAVELENGTHS.get(tube) if tube is not None else None
    return RawScan(
        data=np.column_stack((angles, np.asarray(counts, dtype=np.float64))),
        sample_id=sample_id,
        tube=tube,
        wavelength=wavelength,
    )


def merge_scans(
    segments: Sequence[RawScan | FloatArray],
    overlap: Literal["finer", "first", "last"] = "finer",
) -> FloatArray:
    """Merge angle/count scan segments into one strictly increasing array.

    Segments may have different step sizes and may arrive out of order. Within
    an overlapping angle range, ``overlap="finer"`` (the default) keeps points
    from the segment with the smaller median step and discards the coarser
    segment throughout that range. ``"first"`` or ``"last"`` instead gives
    priority by input order. The result is suitable for :func:`dippa.regrid_to_g`.
    """
    if not segments:
        raise ValueError("need at least one scan segment")
    if overlap not in {"finer", "first", "last"}:
        raise ValueError(f"unsupported overlap policy: {overlap!r}")

    arrays = [np.asarray(segment.data if isinstance(segment, RawScan) else segment) for segment in segments]
    for data in arrays:
        if data.ndim != 2 or data.shape[1] < 2 or data.shape[0] < 2:
            raise ValueError("each segment must have shape (n >= 2, 2)")
        if not np.all(np.diff(data[:, 0]) > 0):
            raise ValueError("each segment must have strictly increasing angles")

    steps = [float(np.median(np.diff(data[:, 0]))) for data in arrays]
    if overlap == "finer":
        priority = sorted(range(len(arrays)), key=lambda index: (steps[index], index))
    elif overlap == "first":
        priority = list(range(len(arrays)))
    else:
        priority = list(reversed(range(len(arrays))))

    selected: list[FloatArray] = []
    claimed_ranges: list[tuple[float, float]] = []
    for index in priority:
        data = arrays[index][:, :2].astype(np.float64, copy=False)
        keep = np.ones(data.shape[0], dtype=bool)
        for lower, upper in claimed_ranges:
            keep &= (data[:, 0] < lower) | (data[:, 0] > upper)
        selected.append(data[keep])
        claimed_ranges.append((float(data[0, 0]), float(data[-1, 0])))

    merged = np.concatenate(selected)
    merged = merged[np.argsort(merged[:, 0], kind="stable")]
    if not np.all(np.diff(merged[:, 0]) > 0):
        raise ValueError("merged scan does not have strictly increasing angles")
    return merged
