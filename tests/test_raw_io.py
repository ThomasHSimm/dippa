"""Regression tests against raw scans shipped with the original DPPA tool."""

from pathlib import Path

import numpy as np
import pytest

from dippa import load_dat, load_udf, merge_scans, regrid_to_g

RAW_DATA = Path(__file__).parent / "fixtures" / "raw"


def test_udf_matches_original_udfcon_output():
    """The committed DAT was produced by the original DPPA ``udfcon.pl``."""
    scan = load_udf(RAW_DATA / "thom_co7.UDF")
    converted = load_dat(RAW_DATA / "thom_co7.dat")

    np.testing.assert_allclose(scan.data, converted, rtol=0.0, atol=1e-10)
    assert scan.data.shape == (2667, 2)
    assert scan.data[0, 0] == pytest.approx(45.005)
    assert scan.data[-1, 0] == pytest.approx(124.985)


def test_udf_metadata_and_loadb_gui_cobalt_wavelength():
    scan = load_udf(RAW_DATA / "thom_co7.UDF")

    assert scan.sample_id == "100gs Sample 2"
    assert scan.tube == "Co"
    assert scan.wavelength == pytest.approx(1.79091)


def test_merge_real_scan_overlap_prefers_finer_segment_then_regrids():
    raw = load_dat(RAW_DATA / "thom_co7.dat")
    coarse = raw[:201:2]
    fine = raw[100:301]

    merged = merge_scans([fine, coarse])

    assert np.all(np.diff(merged[:, 0]) > 0)
    np.testing.assert_array_equal(merged[merged[:, 0] < fine[0, 0]], coarse[:50])
    np.testing.assert_array_equal(
        merged[(merged[:, 0] >= fine[0, 0]) & (merged[:, 0] <= coarse[-1, 0])],
        fine[:101],
    )
    regridded = regrid_to_g(merged, wavelength=1.79091)
    assert np.all(np.diff(regridded[:, 0]) > 0)


def test_merge_overlap_policy_is_exposed():
    raw = load_dat(RAW_DATA / "thom_co7.dat")
    first = raw[:5].copy()
    last = raw[:5].copy()
    first[:, 1] = 1.0
    last[:, 1] = 2.0

    np.testing.assert_array_equal(merge_scans([first, last], overlap="first"), first)
    np.testing.assert_array_equal(merge_scans([first, last], overlap="last"), last)


def test_merge_rejects_non_increasing_segment():
    raw = load_dat(RAW_DATA / "thom_co7.dat")
    with pytest.raises(ValueError, match="strictly increasing"):
        merge_scans([raw[:5][::-1]])
