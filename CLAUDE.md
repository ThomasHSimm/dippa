# CLAUDE.md

Context for any agentic coding tool (Claude Code, Codex) picking this up.
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

`notebooks/` holds runnable examples — currently `01_profiles_walkthrough.ipynb`.
These serve two purposes at once: a fast manual sanity-check while developing
(rerun the notebook, eyeball the plots, no need to reach for a debugger for
"does this look right") and, once the API stabilises, the first thing a new
user is pointed at. Keep committed notebooks with cleared outputs (`jupyter
nbconvert --clear-output`) to keep diffs sane — run them locally to check,
don't commit the run.

## Current state (update this section as work lands)

- `src/dippa/legacy_io.py` — reads legacy `.mat` state files. Done, tested.
- `src/dippa/profiles.py` — pseudo-Voigt, asymmetric pseudo-Voigt, Kα
  doublet, multi-peak + background evaluator. Done, parity-verified
  (R²=0.994 against `tests/fixtures/reference_{fit,data}.mat`).
- `notebooks/01_profiles_walkthrough.ipynb` — runnable walkthrough of
  `profiles.py`, including the real parity check plotted visually. Add a
  new numbered notebook per module as it lands, don't keep bolting onto
  this one indefinitely.
- **Not started**: the actual fitter (bounded nonlinear least squares from
  a real starting guess — everything above only evaluates *known*
  parameters, it doesn't find them yet). This is the next piece of work.
  Also not started: modified Williamson-Hall (`getWH.m` equivalent), the
  contrast factor (cubic case, `Ch00(1 - q*H²)`), and resolving whether
  `onepeak.m`'s fitting path is numerically identical to `pv_tv_aa.m`'s
  plotting path (open question, see `AUDIT.md` §7 and §9).
- Repo not yet pushed anywhere live — this baseline was built in a
  sandboxed environment and handed off as a zip. No CI has run for real yet.

## House style

`ruff` (line length 100), `pytest`, pre-commit hooks already configured.
Docstrings should cite the specific original `.m` file(s) a function is
ported from, the way `profiles.py` does — this isn't decorative, it's what
makes future parity-checking possible without re-reading the MATLAB source
from scratch.
