# TODO

Running list of what's next, what's deferred-but-real, what's exploratory,
and what's an open question worth resolving before it quietly becomes a
bug. Not a roadmap (see `docs/roadmap.qmd` for v0.1 scope) — this is the
working list underneath it. Update as things land or get dropped; a stale
TODO is worse than none.

## Now / next up

- [x] **The fitter.** Built as `background.py` + `fitting.py`. Parity
      result: recovers R²=0.994 on the real reference pattern from a
      starting guess with R²=-1.23. See `AUDIT.md` §11.
- [ ] Modified Williamson-Hall (`getWH.m` equivalent) — mWH-1/2/3 per
      Equation 1 / Table 2 of the 2016 Mater. Des. paper.
- [ ] Contrast factor, cubic case only: `Ch00 * (1 - q*H²)`.
- [x] Decided (not "silently picked"): bound-setting logic in
      `fitting.py` is simplified from `onepeak.m`'s, not a literal
      transcription — documented in the module docstring and `AUDIT.md`
      §11. No "other peaks in this window" freezing logic, because
      decontamination already removes their contribution before the local
      fit runs.

## Open technical questions (not blockers, but need resolving before they bite)

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
      doublet directly during peak fitting. Different philosophy, not a
      strict improvement — Rachinger doesn't assume Kα1/Kα2 share a peak
      shape (current approach does), but is a raw-data preprocessing step
      with its own known distortions. Needs its own comparison, not a
      swap-in.
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
- [ ] Domain name — deliberately not spending anything on this yet;
      GitHub Pages (`thomashsimm.github.io/dippa`) is the front door until
      there's a working v0.1 worth pointing people at.
