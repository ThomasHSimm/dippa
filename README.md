# Dippa

Diffraction peak profile analysis (DPPA) in Python — a rebuild of
[ThomasHSimm/DPPA](https://github.com/ThomasHSimm/DPPA), a MATLAB tool for
Warren–Averbach, Williamson–Hall, and modified Williamson–Hall analysis of
plastically deformed metals.

**Status: pre-alpha, scaffold only.** No fitting or analysis code has been
ported yet. Nothing here should be used for real analysis. See
[`AUDIT.md`](./AUDIT.md) for the technical audit of the original MATLAB
source that this port is based on.

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
scientific core (testable, GUI-independent) from an optional interactive
front end, matching the original app's workflow (`BIGdippa` → peak fitting,
`dippaFC` → Fourier/Williamson-Hall analysis on the fitted peaks).

## Scope for v0.1 (proposed, not yet built)

- Load a diffraction pattern and an instrumental-standard pattern.
- Fit a set of peaks with the pseudo-Voigt profile used in the original tool
  (`pv_tv_aa` / `pk_alpha` / `pk_voigt2`), including Kα1/Kα2 doublet handling
  for Cu/Co/Fe tubes.
- Reproduce (modified) Williamson-Hall analysis (`getWH.m`: variants A/B/C).
- Reproduce one worked example from the original repo's bundled `.mat`
  reference files within a documented numerical tolerance.
- No GUI, no Warren-Averbach Fourier analysis, no HCP/anisotropic contrast
  factors, no CMWP-equivalent functionality yet — see `AUDIT.md` for why
  these are deferred.

## Install (not yet published)

```bash
pip install dippa           # not yet on PyPI — name is reserved/available
pip install "dippa[gui]"    # future: adds PySide6 GUI
```

## Documentation

Full docs, including the [source audit](./AUDIT.md) and
[roadmap](./docs/roadmap.qmd), are published at
[thomashsimm.github.io/dippa](https://thomashsimm.github.io/dippa) (once
the first `docs/**` change lands on `main`).

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
