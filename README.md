# Dippa

An early-stage Python package for reproducible diffraction peak fitting
and (in progress) modified Williamson-Hall analysis of plastically
deformed cubic metals. It rebuilds selected methods from
[ThomasHSimm/DPPA](https://github.com/ThomasHSimm/DPPA) — the original
MATLAB DPPA tool — as independently testable numerical components,
cross-checked against the original tool's own saved reference data.

**Status: pre-alpha.** Implemented and tested so far: legacy `.mat` I/O,
the pseudo-Voigt profile functions (symmetric, asymmetric, Kα doublet),
quadratic background estimation, and the per-peak local-window fitter with
structured diagnostics. Not yet implemented: contrast factors and modified
Williamson-Hall — so **nothing here produces microstructural outputs yet**,
and nothing should be used for real analysis. See
[`AUDIT.md`](./AUDIT.md) for the technical audit of the original MATLAB
source, including what has been verified against real data and what
hasn't.

## Why this exists

The original DPPA/BIGdippa MATLAB app implements methods published in:

- T.H. Simm, P.J. Withers, J. Quinta da Fonseca, *An evaluation of
  diffraction peak profile analysis (DPPA) methods to study plastically
  deformed metals*, Materials & Design 111 (2016) 331–343.
  [doi:10.1016/j.matdes.2016.08.091](https://doi.org/10.1016/j.matdes.2016.08.091)
- T.H. Simm, P.J. Withers, J. Quinta da Fonseca, *Peak broadening anisotropy
  in deformed face-centred cubic and hexagonal close-packed alloys*, J. Appl.
  Cryst. 47(5) (2014) 1535–1551.
  [doi:10.1107/S1600576714015751](https://doi.org/10.1107/S1600576714015751)

The MATLAB code works but is a single-user GUI app: no automated tests, no
package boundary between fitting logic and GUI callbacks, hardcoded paths,
and mixed global/preference-based state. This project separates the
scientific core (testable, GUI-independent) from any future front end,
matching the original app's workflow (`BIGdippa` → peak fitting,
`dippaFC` → Williamson-Hall/Fourier analysis on the fitted peaks).

## Scope for v0.1

- Input: a diffraction pattern already expressed in the documented
  reciprocal-space coordinate (g = 2 sin θ / λ, in Å⁻¹ — see `AUDIT.md`
  §13), with user-supplied approximate peak positions.
- Pseudo-Voigt peak fitting matching the original tool (`pv_tv_aa` /
  `pk_alpha` lineage), including optional Kα1/Kα2 doublet handling for
  Cu/Co/Fe tubes, with per-peak fit diagnostics (optimiser status,
  bound-hit flags, local backgrounds) as first-class output.
- Cubic contrast factors and modified Williamson-Hall (`getWH.m`:
  variants A/B/C) — in progress.
- Reproduction of a worked example from the original repo's bundled
  `.mat` reference files within documented numerical tolerances.
- **Sample patterns only in v0.1** — no instrumental-standard handling
  yet. The original tool always subtracts an instrumental breadth before
  Williamson-Hall (`AUDIT.md` §4), so until that lands, any WH output
  from this package is instrument-uncorrected and will be labelled as
  such: a methodological preview, not analysis-ready numbers.
- No GUI, no Warren-Averbach Fourier analysis, no HCP/anisotropic
  contrast factors, no CMWP-equivalent functionality — see
  `docs/roadmap.qmd` for why each is deferred.

Future scope (not in v0.1, not promised for a date): Warren-Averbach,
HCP symmetry, instrumental correction, variance methods, an optional GUI.

## Install (not yet published)

```bash
pip install dippa           # not yet on PyPI — name is reserved/available
```

## Documentation

Full docs, including the [source audit](./AUDIT.md) and
[roadmap](./docs/roadmap.qmd), are published at
[thomashsimm.github.io/dippa](https://thomashsimm.github.io/dippa) (once
the first `docs/**` change lands on `main`).

## Examples

Runnable walkthroughs in `notebooks/` — `01_profiles_walkthrough.ipynb`
(peak-shape functions, including the real forward-model check against the
original tool's saved output) and `02_fitting_walkthrough.ipynb` (the
per-peak fitter recovering a real pattern from a rough start). Requires
the `examples` extra:

```bash
pip install -e ".[examples]"
jupyter notebook notebooks/
```

## Development

```bash
pip install -e ".[dev]"
pre-commit install    # ruff, detect-secrets, large-file checks
pytest
ruff check .
```

## Relationship to the original MATLAB tool

The original tool remains at
[ThomasHSimm/DPPA](https://github.com/ThomasHSimm/DPPA) (MIT licensed) and is
not being deprecated by this project. Dippa is a from-scratch Python
reimplementation of the same published methods, cross-checked against the
original tool's numerical output where possible — it is not a transpilation
and is not an official CMWP-family tool.

## Citing

If you use Dippa, please cite the original methods papers above until a
software-specific citation (Zenodo DOI) exists for this package.
