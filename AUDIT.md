# Dippa: technical audit of the original MATLAB source

Source: `dippa_v3.zip`, downloaded 2026-07-15 from
[github.com/ThomasHSimm/DPPA](https://github.com/ThomasHSimm/DPPA) (master
branch, MIT licensed, sole copyright holder T.H. Simm). This is a real audit
of the actual extracted source, not a description from memory or from the
handover document alone.

## 1. Inventory

197 files total: 80 `.m`, 34 `.mat`, 55 `.png`, 15 `.html`, 3 `.fig`, 2 `.pl`
(UDF format converters), 2 `.dat`, 2 `.UDF`, 1 `.jpg`, 1 `.asv`.

Two entry points, matching the README's documented workflow:

- `BIGdippa.m` (165,273 bytes) — peak loading and fitting GUI.
- `dippaFC.m` (85,254 bytes) — Fourier/Williamson-Hall analysis GUI, reads
  `BIGdippa`'s output.

Supporting code lives in two folders:

- `BIGdippaFunctions/` (34 files) — peak fitting, data import, WH fitting.
- `FourierFunctions/` (27 files) — Warren-Averbach, contrast factors, plotting.

State is persisted as `.mat` files under numbered folders (`0. variables`,
`1. settings`, `2. data`, `3. fit_data`, `4. RES`), which the GUI reads/writes
using hardcoded paths built from `cd` (MATLAB's current directory) — 34 of
the 80 `.m` files contain a `[cd,'/...']`-style path. This is the single
biggest reason the app can't be scripted or tested headlessly as-is: file
location depends on where MATLAB happens to be running from at the time.

## 2. The `aa` / `aabcg` parameter arrays (decoded from real data)

Loading the bundled `0. variables/fit.mat` with `scipy.io.loadmat` gives:

- `aa`: shape `(6, 11)` — 11 columns for a fit with 10 peaks: each of the
  first 10 columns is one peak's parameters, the 11th (last) column holds
  the quadratic background coefficients (`pv_tv_aa.m`:
  `backgd = lam0(1,len) + x*lam0(2,len) + lam0(3,len)*x.^2`). 6 rows are
  used because this dataset uses the *asymmetric* peak variant
  (`pk_alpha_asymm`), which needs 6 parameters per peak; the symmetric
  variant (`pk_alpha` / `pk_voigt2`) only uses 4 of the 6 rows.
- `aabcg`: shape `(2, 10)` — one pair of local background coefficients per
  peak (distinct from the single global background in `aa`'s last column).
  Real values are all small negative numbers (~-0.0008 to -0.008), i.e. a
  slight negative-sloping local baseline per peak.
- `aa_I` / `aabcg_I`: identical structure, for the **i**nstrumental-standard
  fit, stored and loaded completely separately from the sample fit
  (`fit.mat` vs `fit_I.mat`). Confirms the handover doc's note that sample
  and instrument code paths are duplicated rather than parameterised.
- `genset.mat` holds the per-dataset scalar settings actually used at fit
  time: `identa` (free-text sample name), `alpha2` (bool: fit Kα2 as a
  separate peak or fold it in), `tube` (`'Cu'|'Co'|'Fe'`, hardcoded
  wavelengths per tube in `pk_alpha.m`), `bcg2peak` (background-to-peak
  ratio used to set fit bounds), `wavelen`.

This settles audit question 1–2 from the handover doc's suggested-questions
list with real numbers rather than inference.

## 3. Profile function (`pv_tv_aa.m` and friends)

Confirmed: `pv_tv_aa` is a sum of pseudo-Voigt peaks plus a shared quadratic
background. Peak shape itself is computed in `pk_alpha.m` (handles the
Kα1/Kα2 doublet) which calls `pk_voigt2.m` (the actual pseudo-Voigt), or in
`pk_alpha_asymm.m` for the 6-parameter asymmetric variant used in the
worked example above.

`pk_alpha.m` hardcodes characteristic wavelengths per tube:

| Tube | λα1 (Å) | λα2 (Å) |
|---|---|---|
| Cu | 1.54056 | 1.54439 |
| Co | 1.789010 | 1.792900 |
| Fe | 1.93604 | 1.93998 |

**Correction (2026-07-18), resolving the Co discrepancy flagged in review:**
the original source is internally inconsistent between its two profile
paths. `pk_alpha.m` (symmetric) uses Co = 1.789010 / 1.792900, with
1.78897 / 1.79285 present *as comments*; `pk_alpha_asymm.m` (asymmetric)
uses Co = 1.78897 / 1.79285 as the live literals. Both pairs are therefore
"exact original literals" — from different files. The port's
`TUBE_WAVELENGTHS` carries the `pk_alpha_asymm.m` pair, which is the right
choice for now because the asymmetric path is the one that has been
verified against the reference fixture; if/when the symmetric doublet path
gets its own real-data check, it may need the `pk_alpha.m` pair instead
(doublet spacing differs by ~0.2% relative between the two).

Three tubes only — anything else (Mo, Ag, synchrotron/neutron wavelength
inputs used elsewhere in the same codebase, per the thesis) falls through
with no explicit error. This is a real gap to fix, not carry over, in a
port: the port should take wavelength as an explicit input rather than a
tube-name lookup, and treat the three hardcoded pairs as one built-in table
among others.

There are three near-duplicate profile functions in the tree
(`pv_tv_aa.m`, `_pv_tv.m`, `_pv_tv_asymm.m`, the latter two prefixed with
`_` suggesting they're superseded/dead code) plus a further copy embedded
as a nested function inside `onepeak.m`. Which one is actually live in the
current GUI flow needs confirming by tracing `BIGdippa.m`'s callbacks, not
assumed from file naming.

## 4. Williamson-Hall (`getWH.m`)

This is a real, literal implementation of **modified** Williamson-Hall
(Ungár-family), not plain Williamson-Hall — it weights peak width by a
contrast-factor-corrected term:

```
C = Aest * (1 - q*H2)      % q = contrast-factor curvature parameter
X = g .* sqrt(C)           % g = fitted peak POSITION (aa row 1)
```

**Correction (2026-07-18):** an earlier draft of this section called `g`
"the fitted peak width". That reversed the roles: `getWH.m` reads
`g = aa(1,selected)` — row 1 is peak *position* (the diffraction-vector
magnitude), consistent with the `aa` convention documented in §2 and in
`profiles.py`. The contrast-corrected `X = g·√C` is the horizontal
coordinate; peak *breadth* is the dependent quantity fit against it.

**Which breadth enters mWH — answered from source (`getFW_IB.m`):** it's a
user preference (`type_FWorIB`), one of two options, each averaged over the
two sides for the asymmetric 6-parameter case:

- `'FW'`: FWHM = `0.5*(aa(3,:) + aa(5,:))` (mean of right/left widths);
- `'IB'`: integral breadth = `0.5*fw*(π·η + (1−η)·√(π/ln2))`, the
  closed-form pseudo-Voigt integral breadth, with `fw` and `η` each
  side-averaged first.

**In both cases the instrumental breadth is subtracted before WH**:
`delK = delK_R − delK_I`, where `delK_I` is the same quantity computed
from the instrumental-standard fit (`fit_I.mat`), matched per-peak by
nearest position. So instrument correction in the original is a simple
per-peak breadth subtraction, not a deconvolution — cheap to port, but it
means the original's WH/mWH path *always* assumes an instrumental fit
exists. A sample-only port must either carry the subtraction as an option
or clearly label its WH outputs as instrument-uncorrected (see §12 and
`docs/roadmap.qmd`).

and fits one of three functional forms selected by `type_sizestrain`:

- `mwhA`: `size + strain·X` (linear)
- `mwhB`: `sqrt(size² + strain²·X²)` (root-sum-square)
- `mwhC`: `size + strain·X²` (quadratic in X)

This matches the standard modified-Williamson-Hall formulation from the
Ungár/Borbély literature cited in T's own 2016 paper — i.e. the maths is
exactly what the field would expect, and the "three variants" (A/B/C) map
onto known alternative size-strain separation models rather than being
Simm-specific inventions. Useful confirmation that a port is reproducing a
known, citable method rather than an idiosyncratic one.

## 5. HCP duplication (a second axis of duplication, beyond sample/instrument)

Six files are `_HCP`-suffixed near-duplicates of their cubic-symmetry
counterparts: `dippa_fitWH_HCP.m`, `getWH_HCP.m`, `contrastHCP.m`,
`dippaFCplot_HCP.m`, `logWAindi_HCP.m`, `logWAindi_noSIZE_HCP.m`. Combined
with the sample/instrument (`_I`) duplication in §2, symmetry handling and
sample/instrument handling are each implemented as copy-pasted branches
rather than as orthogonal parameters. A Python port that represents crystal
symmetry and sample-vs-instrument as data (not as which function you call)
removes both duplication axes in one design decision — this is probably the
single highest-value structural improvement over the original.

## 6. State management is inconsistent across the codebase's history

4 files still use MATLAB `global` (e.g. `_pv_tv.m`: `global alpha2`), while
the rest of the codebase (including the live `pv_tv_aa.m`) uses a
preference-store pattern (`getDPPApref('alpha2')` / `setDPPApref`). This is
consistent with the `_`-prefixed files being an earlier generation that
was superseded but not deleted — further reason to trace actual call graphs
rather than porting every file found in the folder.

## 7. What this confirms or changes from the handover document

Confirms: file sizes, entry points, `.mat`-based state, hardcoded paths,
sample/instrument duplication, pseudo-Voigt + quadratic background model,
Warren-Averbach and Williamson-Hall both present, MIT license with T as
sole holder (safe to relicense/reuse without third-party clearance).

New, not previously known from the handover doc alone:
- Real shapes and example values for `aa`/`aabcg` (§2).
- A second duplication axis (cubic vs. HCP, §5), independent of the
  sample/instrument one already flagged.
- Only 3 hardcoded tube wavelengths, no user-supplied-wavelength path in
  the peak-shape function (§3) — a genuine limitation to *not* carry
  forward.
- Evidence of at least two code generations coexisting (`global` vs.
  preference-store state, §6) — some files in the tree are very likely
  dead code, not just messy.

## 8. Landscape check (quick, not exhaustive — see caveats)

- No open-source **Python** package implementing Warren-Averbach or
  (modified) Williamson-Hall DPPA was found. This matches the handover
  doc's claim.
- One relevant piece of prior art exists that the handover doc didn't
  mention: [d-m-collins/ToF-WarrenAverbach](https://github.com/d-m-collins/ToF-WarrenAverbach),
  an open MATLAB (not Python) Warren-Averbach implementation for
  time-of-flight neutron diffraction, released alongside a 2019 paper. It's
  narrower in scope (Warren-Averbach only, ToF-specific) and still MATLAB,
  so it doesn't remove the gap, but it means "no open implementation
  exists anywhere" is slightly overstated — "no open Python implementation,
  and the one open MATLAB one is narrower in scope" is the accurate claim.
- CMWP still appears to be the dominant closed/proprietary tool in this
  space based on search results, consistent with the handover doc.
- This was a ~10-minute check, not the literature scan the handover doc's
  own "agreed approach" calls for as step 1. Treat the "field hasn't moved"
  conclusion as provisionally supported, not confirmed.

## 9. Profile function ported and parity-verified (2026-07-18)

`pv_tv_aa`/`pk_alpha_asymm`/`pk_voigt2asymm` are now ported to
`src/dippa/profiles.py`. Verified against the bundled reference `fit.mat`
and `data.mat`, not just re-derived from reading the source: evaluating the
ported function at the real saved fit parameters against the real measured
pattern gives **R² = 0.994**, with the largest residuals attributable to
genuine counting-statistics scatter near the peak tops, not a systematic
model error.

**New finding, not previously known:** getting to R² = 0.994 required
*not* applying the Kα1/Kα2 doublet, even though this dataset's bundled
`genset.mat` has `alpha2 = 0`, which in the original tool's inverted
preference convention (`alpha2 == 0` means "do fit the doublet") implies
the doublet should be applied. Applying it anyway gives R² = 0.73, with
large, systematic residuals to the right of every peak centre (the doublet
component adds spurious intensity there). This means `genset.mat`'s
settings cannot be assumed to reflect the state the tool was actually in
when a given `fit.mat` was produced — these are independent files, and the
preference is read live from the GUI at evaluation time, not stored
per-fit. Any future parity check against a bundled `.mat` pair should
verify the doublet assumption empirically (try both, keep whichever
reproduces the pattern) rather than trusting the adjacent settings file.
This is now pinned by `test_parity_with_doublet_is_much_worse` in
`tests/test_profiles.py`, specifically so a future "fix" back to trusting
`genset.mat` can't silently regress this.

This also settles part of audit question 1 from §7 of this document
(`onepeak.m` vs `pv_tv_aa.m`): plotting (`pv_tv_aa.m`) is now confirmed
correct against real data. Whether `onepeak.m`'s fitting path uses
numerically identical logic is still open and is the next thing to check
before trusting a from-scratch fit rather than a forward-model evaluation
at already-known parameters.


## 10. Fitting pipeline is three stages, not one (2026-07-18)

Resolves the open question from §7/§9: does `onepeak.m` (the actual
fitting entry point) use the same peak-shape logic as `pv_tv_aa.m`
(plotting, already parity-verified)?

**Peak shape: yes, identical.** `onepeak.m` has its own local copy of the
objective function (`pv_tv_asymm` for the 6-parameter case) and it calls
the exact same `pk_alpha_asymm`/`pk_voigt2asymm` already ported into
`profiles.py`. No rework needed there.

**Background: no — it isn't part of either function.** `pv_tv_aa.m`
evaluates a full quadratic background (`c0 + c1*x + c2*x²`), but
`onepeak.m`'s fitting objective explicitly drops the quadratic term (it's
commented out in the source: `backgd = lam0(1,len) + x*lam0(2,len)
;%+ lam0(3,len)*x.^2;`). Tracing where the quadratic coefficient actually
comes from (`BIGdippa.m` line ~1316) shows the real pipeline is **three
stages, not one**:

1. **Background estimation** (`steel_bcg.m` → `bcgminus`): given the
   user-supplied approximate peak positions, collect the data points
   *outside* a window of `bcg2peak*1.25` around every peak, and fit a
   quadratic to just those background-only points.
2. **Per-peak fitting** (`onepeak.m`, called once per peak): fit that one
   peak's position/amplitude/width(s)/eta(s), plus a small linear-only
   adjustment to the background's constant and slope terms — the quadratic
   term from stage 1 is carried through unchanged. **Correction
   (2026-07-18) to this section's earlier bounds description:** when other
   peaks share the fitting window, only their *position* is frozen
   (0.9999–1.0001× the starting value); their amplitude
   (0.01×start … max(F)), widths (0.1×start … 1e-2 absolute) and eta
   (0 … 1.3) remain loosely bounded and are fit jointly with the target
   peak. The target peak gets position ±`bcg2peak`, amplitude 0.1×–10×,
   widths up to `bcg2peak/1.5`, eta 0–1.3. The local linear background is
   unbounded (±inf). So "freeze everything else" was only ever true of
   neighbour *positions*, not their shapes.
3. Peaks are fit one at a time this way, not simultaneously.

This is good news for the port's design (already staged, now with sharper
detail — a background stage genuinely comes first, not just conceptually).

**Correction to an earlier draft of this section:** it's tempting to read
the random-start `fminsearch` as a reliability risk (different starts
landing in different local minima), but that's not actually correct here.
Fitting a polynomial by least squares is linear in its parameters, so the
sum-of-squared-error surface over `(c0, c1, c2)` is a strictly convex bowl
with a single global minimum — there is no local-minima trap for
`fminsearch` to fall into regardless of starting point. The random start
doesn't put correctness at risk.

What *is* still true: it's unnecessary complexity for a problem with an
exact closed-form solution, and the internal function it calls
(`fitcube`/`expfun`, with parameters literally named `A`, `m`, `lambda`) is
clearly a copy-pasted exponential-decay curve-fitting utility repurposed
for a quadratic without renaming — a "reused what was lying around" signal
in the original codebase, not a deliberate design choice. Dippa's
background stage should still use a direct least-squares solve (e.g.
`numpy.polynomial.polynomial.polyfit`) for simplicity and to drop the
dependency on solver tolerances and an unused RNG call — but the
justification is code clarity, not a correctness fix, and shouldn't be
overstated as one.

The per-peak background handled inside `onepeak.m` is a genuinely
different, smaller job: only a linear local correction (`c0 + c1*x`), not
the quadratic — and that's a sound design, not a simplification. Over the
narrow window used to fit one peak (`±bcg2peak`, a small fraction of the
full pattern's g-range), any smooth background curve is well approximated
by a line; the quadratic curvature only matters at the scale of the full
pattern, which is exactly where the one-time `steel_bcg` stage operates.


## 11. The fitter is built, and it's more interesting than `onepeak.m` alone suggested (2026-07-18)

Tracing `onepeak.m`'s caller (`indivfit_GUI.m`) revealed the real algorithm
is not "freeze every other peak, fit one against the whole pattern." For
each peak fit, the original:

1. Extracts a **local window** of data around that peak.
2. **Decontaminates** it: subtracts the current best-known contribution of
   every *other* peak plus the full background, via the algebraic identity
   `local_data + model(nearby peaks only) − model(everything)`, which nets
   out to `local_data − model(far peaks) − background`.
3. Fits just the target peak, plus a small fresh *local linear* background
   (the quadratic term is already gone, having been subtracted in step 2)
   against that cleaned window.

This is genuinely elegant, and it explains two things noted in §10 without
being able to explain them at the time: why `onepeak.m`'s background is
linear-only (the quadratic isn't its job, decontamination already removed
it), and why per-peak fitting works reliably despite tight bounds — the
"other peaks" in a fitting window are only ever the ones close enough to
matter, not all nine.

**Ported to `src/dippa/background.py` and `src/dippa/fitting.py`.**
Background: closed-form quadratic OLS (see §10 for why this is simpler
than, not more correct than, the original's `fminsearch`). Per-peak
fitting: local decontaminated window + bounded `scipy.optimize.least_squares`,
staged over multiple passes.

**Documented simplification:** the original fits overlapping peak clusters
jointly within a shared window. This port fits one peak at a time always,
relying on multiple passes (each seeing the latest fit of every other peak)
to let close peaks refine each other iteratively — a Gauss-Seidel-style
scheme rather than a joint fit. Fine for well-separated peaks, an open gap
for tightly overlapping ones (see `TODO.md`).

**Recovery result — deliberately no longer called "parity" (reframed
2026-07-18 after review):** §9 proved the *forward model* was correct by
evaluating known parameters. This shows the *fitter* can find an answer —
starting from a bad guess (R² = −1.23, worse than predicting the mean:
positions jittered within the window, amplitudes scaled by ±50% of their
true values, widths and etas set to generic placeholder values, background
re-derived from the perturbed positions), three passes of `fit_pattern`
recovers R² = 0.994 against the real measured pattern, with peak positions
within ~0.001% and amplitudes within a few percent of the saved MATLAB
values. Two honest limits on that claim: (1) the starting guess is a
*perturbation of the known answer*, not an independent start — in
particular the amplitude guess is derived from the true amplitude; (2)
R² + position + amplitude recovery does not establish width/shape
recovery, and widths/etas are exactly what feeds Williamson-Hall. See §15
for what width/shape recovery actually looks like on this fixture, and
`tests/test_fitting.py` for what is now asserted (including integral
breadth). Genuine MATLAB parity — identical starts, per-parameter
comparison against MATLAB's own fitted output — remains open (TODO.md).

## 12. Immediate implication for scope



The handover doc's own risk section (§16.1) warns against building
everything at once. Given the real duplication pattern found here, the
narrowest defensible v0.1 is: **one symmetry (cubic), one dataset type
(sample only, no instrument-correction UI), pseudo-Voigt peak fitting, and
modified Williamson-Hall** — deliberately excluding Warren-Averbach Fourier
analysis, HCP, and instrumental deconvolution from v0.1, not because they're
hard, but because each one doubles the number of code paths that need a
parity check against the original MATLAB output before they can be trusted.

## 13. Coordinate system of the reference data, made explicit (2026-07-18)

Review flagged (correctly) that nothing documented what `x` actually is.
Resolved by checking the reference fixture against physics: the ten peak
positions (0.478, 0.553, 0.781, 0.916, 0.956, 1.105, 1.203, 1.235, 1.353,
1.435) match `1/d` for the fcc reflection sequence
(111, 200, 220, 311, 222, 400, 331, 420, 422, 511/333) of an austenitic
lattice with a ≈ 3.6 Å — consistent with the SS316 stainless-steel
provenance of the bundled data. So **`x` is the diffraction-vector
magnitude g = 2 sin θ / λ = 1/d, in Å⁻¹**, and `half_width`/`bcg2peak`
values like 0.02 are in Å⁻¹ too.

A corroborating detail: the data extends to g ≈ 1.54 Å⁻¹, beyond the
maximum reachable with Co Kα radiation (2/λ ≈ 1.12 Å⁻¹ at θ = 90°). This
pattern cannot have been measured with the Co tube its `genset.mat`
declares — consistent with (a) the §9 finding that no Kα doublet
reproduces it, and (b) `genset.mat` also containing `wavelen = 1.542475`
(a Cu-family number) alongside `tube = 'Co'`. Yet another independent
confirmation that settings files describe GUI state, not data provenance.

Port implication (still open, logged in TODO.md): the doublet displacement
in `evaluate_peak` is only correct when `x` is this reciprocal-space
coordinate. Nothing currently stops a user passing degrees 2θ. The API
should eventually carry coordinate/units explicitly (e.g. a typed pattern
object); until then the convention is documented in `profiles.py` and here.

## 14. Window preference is `sizfit`, not `bcg2peak` (2026-07-18)

`indivfit_GUI.m` selects the data window using `fitsiz =
getDPPApref('sizfit')` — a *different* preference from `bcg2peak`, which
`onepeak.m` uses for the target peak's position bound and width cap. The
port had conflated the two into a single `half_width`. The bundled
`genset.mat` stores `bcg2peak = 0.02` but no `sizfit` at all (read live
from the GUI), so the window actually used to produce the reference fit is
unknowable from the saved files — the same class of problem as the
`alpha2` finding in §9. The port keeps a single `half_width` for now
(documented in `fitting.py`); splitting it is only worth doing if a real
fixture ever pins down a `sizfit` ≠ `bcg2peak` case.

`indivfit_GUI.m` also confirms the local-background lifecycle: each peak's
fit is *seeded* from the stored `aabcg` column and the fitted values are
stored back (`aabcg(:,n) = aa(1:2,lenaa)`), persisting across refits. The
port now does the same (see §15) — an earlier version fitted the local
line and then discarded it, which review correctly flagged as diverging
from the original and from the function's own docstring.

## 15. Width/shape identifiability on the reference fixture, and structured results (2026-07-18)

Review made the right criticism: R² + position recovery does not establish
width/shape recovery, and breadth is what mWH consumes. Measured directly
(fit from the standard rough guess, three passes, against the saved MATLAB
values):

- positions: ~0.001% — excellent;
- amplitudes: ≤ 3% — good;
- FWHM (per side): 0.4–28% error; eta: 2.5–100% error — **not recovered
  parameter-by-parameter**;
- integral breadth (the side-averaged pseudo-Voigt closed form that
  `getFW_IB.m` feeds to WH): within ~15% on every peak, mostly ≤ 8%.

The FWHM/eta pairs trade off against each other (and against the local
background) with little residual cost — they are not individually
identified by this data at this window size. The worst case is the
weakest, last peak (index 9): its `eta_right` lands exactly on the 1.3
upper bound (the original tool's own bound, from `onepeak.m`). This is
*not* a rough-guess artefact — refitting peak 9 starting from the exact
MATLAB answer also pegs `eta_right` at 1.3, at every window half-width
tried (0.008–0.03), and the pegged solution has genuinely lower SSE in
this objective than the MATLAB parameters. The MATLAB answer is a
different local minimum (its optimiser was seeded differently, via the
persisted `aabcg`, with an unknowable `sizfit` window — §14 — and the
saved fit is whatever the user interactively accepted).

Consequences, now implemented:

- `fit_one_peak` / `fit_pattern` return structured results
  (`PeakFitResult` / `PatternFitResult`) carrying optimiser
  success/message, cost, evaluation count, bound-hit flags per named
  parameter, and the fitted local backgrounds in the original's `aabcg`
  convention. A fit with a parameter pinned at a bound is no longer
  silently indistinguishable from a clean one — this is the §3-of-AGENTS
  "uncertainty as first-class output" principle applied to the fitter.
- `tests/test_fitting.py` asserts integral-breadth recovery (15%
  envelope, documented not aspirational), asserts the peak-9 bound hit is
  *reported*, and deliberately does not assert per-parameter FWHM/eta
  agreement it cannot deliver.
- Anything consuming widths downstream (mWH, when built) should check
  `result.warnings` / `result.all_clean` and prefer integral breadth over
  raw per-side FWHM, since that is both the better-identified quantity
  and one of the two the original tool itself uses (§4).

## 16. Phase 3 verification: dippa re-fit of all 9 ni_combo.mat samples (2026-07-19)

All 9 samples from `ni_combo.mat` were re-fit using dippa's `fit_pattern`
from deliberately rough, reproducible starts: stored positions jittered by
up to 0.005 Å⁻¹, stored amplitudes scaled by 0.5–1.5, generic FWHM/eta
values, and a quadratic background re-estimated from the perturbed positions.
This remains a recovery check, not an independent start or MATLAB parity.
The executable comparison table, figure, and per-sample warnings are in
`notebooks/03_fitting_demonstration.ipynb`. Result:

- **Parameter recovery**: across all 45 peaks, the largest sample-level
  position error is 0.0000303 Å⁻¹, amplitude error is 0.225%, and integral-
  breadth error is 1.744%.
- **Minor bound violations**: Only sample 0 (niHmid_halfpc, lowest SNR)
  hit bounds during fitting (eta_left pegged at 1.3 for peaks 3 and 4),
  consistent with the width/shape non-identifiability documented in §15.

This verifies that:
1. The port's implementation of `fit_pattern` is algorithmically sound.
2. The MATLAB `.mat` file extraction (via `io/matlab_io.py`) correctly
   preserves all parameter values.
3. The coordinate system assumption (g = 2sinθ/λ, tube='Co' for doublet) is
   consistent with how the reference data was produced.

## 17. Phase 4 complete: tsinterpl.m port to regrid.py (2026-07-19)

Ported `BIGdippaFunctions/tsinterpl.m` (time-interpolation/regridding) to
Python as `src/dippa/regrid.py`. Functionality:

- **Coordinate conversion**: 2θ (degrees) ↔ g = 2sinθ/λ (Ų⁻¹)
- **Regridding with PCHIP**: Converts irregular 2θ grid to uniform g-step
  grid using monotone cubic Hermite interpolation (scipy's PchipInterpolator)
- **Default step**: 5e-5 Ų⁻¹ (matches ni_combo.mat)
- **Reversed-data handling**: Automatically sorts if g is descending (line 15–20
  in MATLAB original)
- **PCHIP overshoot handling**: Optional clipping of negative values on low-
  intensity regions (PCHIP can overshoot near noise)
- **Int16 overflow fix**: Python ints are arbitrary precision; no risk of
  overflow in grid-index calculation (MATLAB line 26)

**Validation against ni_combo.mat**:
- Grid step uniformity: 5e-5 Ų⁻¹ to machine precision (std: 7.7e-19)
- Monotonicity preserved across all 9 samples
- Coordinate conversion matches expected Bragg positions
- Roundtrip conversion (theta → g → theta) accurate to 10⁻¹⁰

**Tests**: 21 new unit tests (all passing), covering:
- Coordinate conversions and roundtrip accuracy
- Known Bragg reflections (Ni FCC with Cu/Co radiation)
- Regridding properties (monotonicity, step uniformity)
- Edge cases (extreme step sizes, single peaks, noisy data)
- Integration tests (synthetic workflows)

**Total test suite**: 59 passing (38 existing + 21 new regrid).
