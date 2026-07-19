# AGENTS.md

Context for any agentic coding tool (Claude Code, Codex, or otherwise) picking
this up. Read by Codex CLI natively; Claude Code reads it as additional
context too, so this is the single source of truth rather than a
Codex-specific duplicate of anything else.

Read this before making changes — several decisions here came from getting
things wrong first, and re-litigating them wastes time. Full history is in
`AUDIT.md`; this file is the condensed, actionable version.

## What this project is

A Python port of `github.com/ThomasHSimm/DPPA` (MATLAB, `BIGdippa`/`dippaFC`)
— diffraction peak profile analysis (DPPA) for studying dislocation density
and crystal size in plastically deformed metals. The original author's own
papers and PhD thesis are the ground truth for both the physics and for
which methods are and aren't trustworthy — see `AUDIT.md` and the `docs/`
site for citations. This is a from-scratch reimplementation cross-checked
against the original tool's output, not a transpilation.

## Non-negotiable design principles (established, don't re-decide these)

1. **Nothing is "ported" until it reproduces a real number from the
   original tool.** Not a synthetic test, not "the math looks right" — an
   actual saved `.mat` fit against actual saved `.mat` data, with a
   pytest-enforced R² threshold. `tests/test_profiles.py` is the template
   for this pattern; follow it for every new module.
2. **Staged, inspectable fitting — never one global optimizer over
   background + size + strain + arrangement simultaneously.** This was a
   deliberate reaction to a documented failure mode in CMWP (the field's
   dominant closed tool): folding background into the same nonlinear fit as
   the physical parameters lets them trade off against each other
   invisibly. Fit in discrete, checkable steps, matching the original
   tool's INDI-style approach, not its log-ALL/whole-profile style.
3. **Report uncertainty and parameter sensitivity as first-class output,
   not an afterthought.** Don't return a bare number where the original
   literature shows that number is sensitive to background choice, fitting
   range, or starting values (see `AUDIT.md` on the Krivoglaz-Wilkens
   curvature problem and the M2/M4 disagreement in Variance-B). If a value
   is poorly constrained, say so in the output.
4. **v0.1 scope is deliberately narrow: cubic symmetry only, sample-only
   (no instrument-correction UI), pseudo-Voigt fitting, modified
   Williamson-Hall with the contrast factor.** Warren-Averbach, HCP,
   Variance-A/B, and anything GUI are explicitly out of scope for now — see
   `docs/roadmap.qmd` for why each one is deferred and what would need to
   be true before adding it. Don't expand scope without updating that page.

## Known traps (each one cost real debugging time — don't reintroduce them)

- **`genset.mat`-style settings files do not reliably describe how a saved
  fit was actually produced.** The original tool reads its `alpha2`
  (Kα-doublet) preference live from the GUI at evaluation time, not from
  the fit's own saved state. Always verify doublet-on vs doublet-off
  empirically against the real pattern rather than trusting an adjacent
  settings file. See `profiles.py` module docstring and `AUDIT.md` §9.
- **The `alpha2` flag in the original MATLAB source is inverted from what
  you'd guess**: `alpha2 == 1` means *no* doublet, `alpha2 == 0` means fit
  the doublet. This port avoids the inversion entirely by using
  `tube: str | None` instead (`None` = no doublet). Do not reintroduce a
  boolean flag with MATLAB's sense — it's the kind of thing that's easy to
  flip by accident and won't show up in a naive test.
- **The asymmetric peak parameter order cross-indexes FWHM/eta between
  sides** — `[x0, amplitude, fwhm_right, eta_left, fwhm_left, eta_right]`,
  not grouped by side. This is confirmed correct against real data, not a
  guess — don't "clean it up" without re-running the parity test.
- Only three X-ray tubes are supported (Cu, Co, Fe), hardcoded wavelength
  pairs, no arbitrary-wavelength path. This is a real limitation of the
  original tool, carried over deliberately (see `AUDIT.md` §3) rather than
  silently "fixed" — if you add a general-wavelength path, keep the
  hardcoded table as a named preset, don't remove it.

5. **The GUI is an optional install extra, never a core dependency.**
   `pip install dippa` must stay lightweight (numpy/scipy/pydantic only);
   `pip install "dippa[gui]"` pulls in PySide6. This is already reflected in
   `pyproject.toml`. Two rules that make the split real rather than
   aspirational, both still to be enforced once GUI code exists:
   - No module under `src/dippa/*.py` at the top level may `import PySide6`
     (or anything GUI-related), directly or transitively. GUI code lives
     entirely under `src/dippa/gui/`, and imports *from* core modules, never
     the other way round.
   - The `dippa-gui` console-script entry point must catch the `ImportError`
     from a missing PySide6 and print "run `pip install dippa[gui]`", not
     surface a raw traceback.
   - Add a CI job that installs `dippa` *without* `[gui]` and does
     `python -c "import dippa"` — this is the cheap, automatic way to catch
     a regression where some core module accidentally pulls in Qt.

## Notebooks

`notebooks/` holds runnable examples — `01_profiles_walkthrough.ipynb`
(profile functions) and `02_fitting_walkthrough.ipynb` (the fitter,
including reading the fit diagnostics). `03_fitting_demonstration.ipynb`
checks rough-start recovery across all nine nickel samples, and
`04_williamson_hall.ipynb` demonstrates integral-breadth instrumental
subtraction plus classical/mwhA/mwhB/mwhC fitting.
These serve two purposes at once: a fast manual sanity-check while developing
(rerun the notebook, eyeball the plots, no need to reach for a debugger for
"does this look right") and, once the API stabilises, the first thing a new
user is pointed at. Keep committed notebooks with cleared outputs (`jupyter
nbconvert --clear-output`) to keep diffs sane — run them locally to check,
don't commit the run.

## Current state (update this section as work lands)

- `src/dippa/io/legacy_io.py` — reads legacy `.mat` state files. Done, tested.
- `src/dippa/profiles.py` — pseudo-Voigt, asymmetric pseudo-Voigt, Kα
  doublet, multi-peak + background evaluator. Done, forward-model-verified
  (R²=0.994 evaluating the saved MATLAB fit against the saved data,
  `tests/fixtures/reference_{fit,data}.mat`). Coordinate convention: `x`
  is reciprocal-space g = 2sinθ/λ = 1/d in Å⁻¹ — established from the
  fixture itself, see `AUDIT.md` §13; the doublet displacement is only
  valid in that coordinate.
- `notebooks/01_profiles_walkthrough.ipynb` — runnable walkthrough of
  `profiles.py`, including the real parity check plotted visually. Add a
  new numbered notebook per module as it lands, don't keep bolting onto
  this one indefinitely.
- `src/dippa/background.py` — closed-form quadratic background fit. Done,
  tested against synthetic data and sanity-checked against the real
  reference fixture.
- `src/dippa/fitting.py` — the per-peak fitter. Done, with a
  single-fixture *recovery* result (deliberately not called "parity" — see
  `AUDIT.md` §11/§15): starting from a rough perturbation of the known
  answer (R²=-1.23), recovers R²=0.994, positions to ~0.001%, amplitudes
  to a few percent, integral breadth to ≤15%. Per-side FWHM/eta are NOT
  individually identified on this data (they trade off; one peak pins eta
  at the 1.3 bound — real, diagnosed, reported via the result object, see
  `AUDIT.md` §15). Returns `PatternFitResult` with per-peak diagnostics
  and `aabcg` local backgrounds — check `.warnings` before consuming
  widths downstream. Architecture is local-decontaminated-window fitting,
  not "freeze other peaks and fit against the whole pattern" — see
  `AUDIT.md` §11 before touching this module. Known simplification: fits
  one peak at a time even when windows overlap (the original fits
  window-sharing peaks jointly, with neighbour positions frozen but
  shapes loosely bounded — `AUDIT.md` §10) — see `TODO.md`.
- `src/dippa/contrast.py`, `breadth.py`, `williamson_hall.py` — cubic
  contrast factors, FW/IB breadths with nearest-position instrumental
  subtraction, and classical/mwhA/mwhB/mwhC fitting. Done, with typed
  explicit HKLs, exclusion policies for flagged/nonpositive breadths,
  Jacobian covariance and 95% confidence intervals. The contrast factor is
  parity-tested against the five stored C values in
  `SS316_logINDI_RES.mat`; see `AUDIT.md` §19. The nine-sample nickel
  demonstration reports that the half-percent sample retains only two
  trustworthy peaks and is not forced through a three-parameter fit.
- Repo is live at `github.com/ThomasHSimm/dippa` with CI
  (`.github/workflows/ci.yml`) and a Quarto docs publish workflow. The
  original baseline was built sandboxed and handed off as a zip; this
  section previously said "not pushed anywhere" — stale, corrected.

## House style

`ruff` (line length 100), `pytest`, pre-commit hooks already configured.
Docstrings should cite the specific original `.m` file(s) a function is
ported from, the way `profiles.py` does — this isn't decorative, it's what
makes future parity-checking possible without re-reading the MATLAB source
from scratch.
