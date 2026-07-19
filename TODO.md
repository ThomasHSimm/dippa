# TODO

Running list of what's next, what's deferred-but-real, what's exploratory,
and what's an open question worth resolving before it quietly becomes a
bug. Not a roadmap (see `docs/roadmap.qmd` for v0.1 scope) — this is the
working list underneath it. Update as things land or get dropped; a stale
TODO is worse than none.

## Now / next up

### Completed (Phases 1–3)

- [x] **matlab_io module** (`src/dippa/io/matlab_io.py`): `MatlabSample` dataclass
      and extraction functions for scipy.io.loadmat data. Tested on minimal
      fixture (ni_combo_minimal.mat: 37.8 KB, 1 sample, decimated).
- [x] **The fitter.** Built as `background.py` + `fitting.py`.
      Single-fixture recovery result (not MATLAB parity — see `AUDIT.md`
      §11/§15): R²=0.994 from a rough start (R²=-1.23), positions ~0.001%,
      amplitudes ≤3%, integral breadth ≤15%. Per-side FWHM/eta not
      individually identified on this data.
- [x] **Structured fit results** (`PeakFitResult`/`PatternFitResult`):
      optimiser success, cost, bound-hit flags, local backgrounds in the
      original's `aabcg` convention (persisted across passes and seeded,
      matching `indivfit_GUI.m`). Added after review flagged that the
      bare-array API violated the uncertainty-first principle and that
      local backgrounds were fitted then discarded.
- [x] **Phase 3 validation**: All 9 ni_combo.mat samples re-fit with dippa
      from rough reproducible starts. Position, amplitude, and integral-
      breadth recovery are documented in `AUDIT.md` §16 and the fitting
      demonstration notebook.
### Completed (Phase 2 notebook documentation)

- [x] **Demonstration notebook fixes** (`notebooks/01_fitting_demonstration.ipynb`):
      - [x] Fix coordinate axis label: 2θ (degrees) → g (Å⁻¹) in axis labels and text
      - [x] Add dataset documentation cell explaining coordinate system, material, radiation
      - [ ] Remove claim about peak 4 artifact/overlap — confirmed distinct via spacing analysis
      - [ ] Add hkl indexing (111/200/220/311/222) to peak labels (partial)
      - [x] Compute all narrative numbers as f-strings from actual data (done in plots)
      - [ ] Fix eta interpretation (>1.3 is fitting artifact, not physical)
      - [x] Relabel samples with (strain) notation in sample names
      - [x] Add residual/√counts normalization alongside absolute residuals
      - [x] Drop any strain-position-shift claims unsupported by low-SNR low-strain data

### Completed (Phase 4 regridding/interpolation)

- [x] **Phase 4: Interpolation and regridding** (`tsinterpl.m` port, `regrid.py`):
      - [x] Port `tsinterpl.m` time-interpolation logic to Python (`regrid.py`)
      - [x] Handle 2θ ↔ g coordinate conversion (theta_to_g, g_to_theta)
      - [x] Implement PCHIP (Piecewise Cubic Hermite Interpolation) regridding
      - [x] Fix int16 overflow edge cases (automatic via Python arbitrary-precision ints)
      - [x] Add comprehensive test suite (21 tests, all passing)
      - [x] Validate against ni_combo.mat (grid uniformity to machine precision)

### Future (Phases 5+)
- [x] Modified Williamson-Hall (`getWH.m` equivalent) — mWH-1/2/3 per
      Equation 1 / Table 2 of the 2016 Mater. Des. paper. Prerequisites
      now settled from source (`AUDIT.md` §4): abscissa `X = g·√C` with
      g = peak *position*; breadth is side-averaged FWHM or integral
      breadth (`getFW_IB.m`), and the original always subtracts the
      instrumental breadth. Implemented with explicit exclusion policies,
      covariance/CIs, and integral breadth as the default.
- [x] Contrast factor, cubic case only: `Ch00 * (1 - q*H²)`; exact stored
      C-vector parity against `SS316_logINDI_RES.mat`.
- [ ] **Second-material end-to-end candidate** — run the staged breadth +
      mWH workflow against `data/SSnew_interpBCG_fit.mat` once its peak and
      instrument-fit provenance is pinned down.
- [ ] **Paper-value comparison (2016 Mater. Des. tables) — author task.**
      Compare computed nickel strain-series trends and intervals manually;
      do not encode literature q values as test truth.
- [ ] **Warren–Averbach remains out of scope.** The SS316 WA-track result is
      used only as a stored contrast-factor parity oracle, not as authority
      to port the Fourier workflow.
- [x] Decided (not "silently picked"): bound-setting logic in
      `fitting.py` is simplified from `onepeak.m`'s, not a literal
      transcription — documented in the module docstring and `AUDIT.md`
      §10/§11. No joint fitting of window-sharing peaks (the original
      freezes neighbour positions but fits their shapes jointly).

## Open technical questions (not blockers, but need resolving before they bite)

- [ ] **Real multi-step overlap fixture for `merge_scans`** — the committed
      `thom_co7` instrumental-standard scan has one constant step. The merge
      policy is tested with overlapping slices of that real scan, but needs
      validation against a genuine fine-step peak-window plus coarse survey
      scan once one is available.
- [ ] **Genuine MATLAB parity for the fitter** — the current test is
      single-fixture *recovery* from a perturbed-truth start. Real parity
      needs MATLAB and Python runs from identical starting parameters,
      per-peak width/shape/amplitude/background comparison with tolerances
      sized by downstream effect, and multiple starting points. Requires
      running the original tool, so parked until MATLAB access.
- [ ] **Peak 9 eta identifiability** — the weakest peak pins `eta_right`
      at the 1.3 bound even when refit from the exact MATLAB answer, at
      every window width tried; the pegged solution beats the MATLAB
      values on SSE in this objective (`AUDIT.md` §15). Now *reported* via
      diagnostics rather than silent. Options if it needs fixing rather
      than flagging: joint refit with a robust loss, penalising η > 1,
      or accepting integral breadth (stable to ≤12%) as the deliverable.
- [ ] **Coordinate/units are convention, not code** — `x` is g in Å⁻¹
      (`AUDIT.md` §13) but nothing stops a caller passing 2θ degrees,
      which would silently break the doublet displacement and every
      half-width default. A typed pattern object (coordinate kind, units,
      wavelength, intensities) is the eventual fix; before the public API
      settles, not after.
- [ ] **`sizfit` vs `bcg2peak`** — the original uses two preferences
      (data-window half-width vs bounds/width-cap); the port conflates
      them into one `half_width` (`AUDIT.md` §14). Only worth splitting if
      a fixture ever pins down a case where they differ.
- [ ] The **symmetric** (4-parameter, non-asymmetric) peak case has never
      been parity-tested against real data — only the 6-parameter
      asymmetric case has (R²=0.994, see `AUDIT.md` §9). Don't assume it
      also works; check it the same way once there's a real symmetric-peak
      reference case.
- [ ] Confirm `_pv_tv.m` / `_pv_tv_asymm.m` (underscore-prefixed, using
      `global` rather than `getDPPApref`) are genuinely dead code, not a
      second live path — likely dead per the state-management finding in
      `AUDIT.md` §6, but "likely" isn't "confirmed."
- [ ] `genset.mat`-style settings files cannot be trusted to reflect how a
      saved fit was actually produced (confirmed once, see `AUDIT.md` §9)
      — every new bundled reference fixture needs its doublet/background
      assumptions checked empirically, not read off the settings file.
- [ ] **Overlapping peak clusters are fit one-at-a-time with iterative
      refinement (multiple passes), not jointly**, unlike the original
      tool. Works fine for well-separated peaks (validated on the real
      reference pattern, whose closest peaks are ~0.04 apart). Untested on
      genuinely close/overlapping peaks — find or construct a test case
      before trusting this on a busier pattern than the reference one.

## Deferred but real (Stage 2 candidates — evidenced, not yet scheduled)

- [ ] **SNIP background algorithm** as an alternative to `steel_bcg`'s
      approach — doesn't require peak positions to be known first (unlike
      the current windowed-quadratic method), well-established across
      XRD/XRF/gamma spectroscopy, simple to implement. Would let
      background estimation run *before* peak-finding rather than after.
- [ ] **Rachinger-style Kα2 stripping** (1948, refined 1975, still being
      revisited as recently as 2026) as an alternative to fitting the
      doublet directly during peak fitting. Different philosophy
      (preprocessing vs model fitting), not a strict improvement — and
      *not* assumption-free: classical Rachinger assumes the Kα1 and Kα2
      components share the same line-profile shape, a fixed intensity
      ratio (≈2:1) and a known angle-dependent separation, i.e. much the
      same physical assumptions as the current linked-doublet fit, applied
      to raw data instead of inside the model. (An earlier note here
      claimed Rachinger avoids the same-shape assumption — wrong,
      corrected 2026-07-18.) Known distortions on the Kα2 side. Needs its
      own comparison, not a swap-in.
- [ ] **Variance-B method** — real accuracy (comparable to the best
      Fourier methods) but historically unreliable (2/15 peaks converged
      in the thesis). Likely explanation: needs synchrotron/neutron-grade
      peak-to-background ratio, not fixable by a better optimizer alone
      (M4 is a 4th-order, tail-weighted moment). If built: check
      peak-to-background on load and warn/refuse rather than silently
      returning a bad number on lab-quality data.
- [ ] **τ-plot visualization** — plot any DPPA output (broadening, lattice
      strain, planar-fault fraction) against measurement angle τ rather
      than just axial/transverse. Self-authored, already validated
      (τ-plot paper), low-risk since it's a visualization layer on
      existing outputs, not new fitting math. Credit explicitly as T's own
      method.

## Exploratory / uncertain (may or may not happen — do not scope hours against these yet)

- [ ] Crystal-plasticity-derived, slip-system-resolved contrast factor.
      Real literature precedent, real known limitation (CP-predicted
      contrast factor only partly explains real broadening anisotropy —
      T's own 2018 Crystals review). Check whether thesis-era
      per-slip-system contrast factor code still exists before assuming
      this is cheaper than starting from scratch.
- [ ] Read the 2024 DDD + crystal-plasticity + X-ray-ray-tracing paper
      (closest existing prior art to the OpenDiS idea) before committing
      any effort to synthetic-diffraction validation — may change scope or
      just become the citation for "here's the closest existing work."
- [ ] OpenDiS scoping spike (run the Frank-Read example, confirm it works,
      note what a dislocation-network export looks like) — cheap, can
      happen anytime, doesn't touch Dippa's codebase. The genuinely hard
      part (dislocation network → synthetic diffraction pattern, Bertin
      spectral framework) has no existing code to build on and would be
      implemented from the paper.
- [ ] Kikuchi-band-width DPPA (EBSD analogue of line-broadening analysis)
      — explicitly *not* a build target (see prior discussion: confounds
      like detector MTF and structure-factor-dependent band-minimum shifts
      make this a fifteen-year-stalled problem, not a quick extension).
      Keeping here only as a note that it's real background knowledge for
      conversations, not something to schedule.

## Docs & literature

- [ ] Once v0.1 is real: replace the `docs/audit.qmd`/`roadmap.qmd`-led
      front page with a working example leading, per the earlier framing
      decision — a project that leads with "here's what's not done yet"
      undersells itself once there's something that works.
- [ ] Read the 2014 J. Appl. Cryst. paper directly (only seen secondhand
      via the 2018 Crystals review so far) — fills the last gap in the
      core methods-paper cluster.
- [ ] Pull the 2010 Materials Science Forum paper (Simm, Withers, Ungár,
      Fonseca) — short conference paper, but a real citable artifact of
      the Ungár collaboration, worth having in hand for the eventual
      literature-review writeup.
- [ ] Spike `quartodoc` for API reference generation on the docs site —
      candidate, fit/effort unverified, worth a short check before
      committing to it over hand-written reference pages.
- [ ] MATLAB→Python migration notes page, for anyone coming from the
      original tool.

## Repo structure (later — noted, not scheduled)

- [ ] Move the `docs/` Quarto site out of the `dippa` repo into its own
      repo, hosted at dmata.co.uk (re-purchased) rather than GitHub Pages.
      Likely bigger than just Dippa's docs — dmata.co.uk as a hub covering
      more than one project is plausible, so don't scope this as a Dippa
      task specifically until the shape of that is clearer. Revisit once
      there's real docs content worth moving, not before.

## Housekeeping / packaging (mostly decided, not all implemented)

- [ ] `dippa-gui` console-script entry point — must catch the `ImportError`
      from a missing PySide6 and point at `pip install dippa[gui]`, not
      surface a raw traceback. Not written yet.
- [ ] CI job that installs `dippa` *without* `[gui]` and does
      `python -c "import dippa"` — cheap automatic guard against a core
      module accidentally importing Qt. Not written yet.
- [x] GUI packaging split decided and reflected in `pyproject.toml`
      (`[gui]` optional extra) — mechanism done, enforcement (above two
      items) not.
- [x] PyPI name `dippa` — confirmed available, low risk of confusion
      (the existing `dppa` package is an unrelated protein-analysis tool).
- [x] Domain — dmata.co.uk re-purchased (see "Repo structure" above);
      GitHub Pages (`thomashsimm.github.io/dippa`) stays the front door
      for Dippa's own docs until there's a working v0.1, and the earlier
      "not spending anything on a domain" note here is superseded.
