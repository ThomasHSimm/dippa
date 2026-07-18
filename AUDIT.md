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
X = g .* sqrt(C)           % g = fitted peak width (aa row 1)
```

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


## 10. Immediate implication for scope

The handover doc's own risk section (§16.1) warns against building
everything at once. Given the real duplication pattern found here, the
narrowest defensible v0.1 is: **one symmetry (cubic), one dataset type
(sample only, no instrument-correction UI), pseudo-Voigt peak fitting, and
modified Williamson-Hall** — deliberately excluding Warren-Averbach Fourier
analysis, HCP, and instrumental deconvolution from v0.1, not because they're
hard, but because each one doubles the number of code paths that need a
parity check against the original MATLAB output before they can be trusted.
