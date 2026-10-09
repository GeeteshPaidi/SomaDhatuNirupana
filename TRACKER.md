# TRACKER.md — ASR rework progress

Branch: `asr-rework`. Plan: [PLAN.md](PLAN.md). Rules: [CLAUDE.md](CLAUDE.md).
New rework code lives in `src/asr/`; rework outputs go to `outputs/asr/` (gitignored).

---

## Phase 0 — Orientation

- [ ] **0.1** Map the repo (summary below).
  - Status / Notes: in progress
- [ ] **0.2** Confirm the environment runs on one month / a few hundred spectra and reproduces the old behaviour; record old ratios as baseline.
  - Status / Notes:
- [ ] **0.3** Check FITS headers (angles, exposure, gain, channels, response/efficiency info); check for local XSM data.
  - Status / Notes:

### 0.1 Repo map

_(filled in by task 0.1)_

## Phase 1 — Forward model and fitting

- [ ] **1.1** Solar-spectrum model (reference flare spectrum, configurable temperature).
  - Status / Notes:
- [ ] **1.2** Rewrite primary fluorescence per PLAN §4; unit-test against hand-computed values for one element.
  - Status / Notes:
- [ ] **1.3** Secondary fluorescence: concentration-coupled, or drop with justification.
  - Status / Notes:
- [ ] **1.4** Detector response: binned Gaussian over ±5σ, normalised; unit test for line-area conservation.
  - Status / Notes:
- [ ] **1.5** Reparameterise for identifiability (Si ≡ 1, ratios + amplitude, O fixed/excluded, 1.0–7.5 keV, continuum term).
  - Status / Notes:
- [ ] **1.6** Poisson/Cash fit on raw counts with background model; covariance, uncertainties, GoF, quality flag.
  - Status / Notes:
- [ ] **1.7** GATE — synthetic recovery test (bias/scatter per ratio, condition numbers, cross-temperature systematic). Do not proceed until Al/Si and Mg/Si recover acceptably.
  - Status / Notes:
- [ ] **1.8** Smoke test on one real month; check sanity ranges; diagnostic plots (high/median/low-quality fits).
  - Status / Notes:

## Phase 2 — Full run (ask Geetesh before starting; ~8 h with 14 workers last time)

- [ ] **2.1** Apply new selection criterion; record counts at each stage (raw → dayside → QC → fitted).
  - Status / Notes:
- [ ] **2.2** Full fit; save per-spectrum results with uncertainties and flags.
  - Status / Notes:
- [ ] **2.3** Selection-threshold sensitivity (subset OK) and solar-temperature sensitivity.
  - Status / Notes:

## Phase 3 — Maps

- [ ] **3.1** Rebuild element/Si maps; observed vs smoothed/infilled coverage; fitted vs all footprints; equal-area stats.
  - Status / Notes:
- [ ] **3.2** Quality layers per map: observation count, mean fit uncertainty, infill mask.
  - Status / Notes:
- [ ] **3.3** State effective resolution honestly (footprint scale, not grid spacing).
  - Status / Notes:

## Phase 4 — Validation

- [ ] **4.1** Apollo/Luna (9 sites): reference table with sources; pre-fixed radius sensitivity (1°, 2°, 5°); bootstrap 95% CIs; leave-one-site-out; propagated per-site uncertainties.
  - Status / Notes:
- [ ] **4.2** Lunar Prospector GRS comparison (common 5° grid; Fe, Ti, Al, Ca, Mg vs X/Si; Spearman + Pearson; block bootstrap; scatter plots).
  - Status / Notes:
- [ ] **4.3** Decide which elements the paper claims, based only on 4.1 + 4.2. Log decision.
  - Status / Notes:

## Phase 5 — Paper rewrite (Elsevier `elsarticle`, ASR)

- [ ] **5.1** Check current ASR Guide for Authors; record requirements here.
  - Status / Notes:
- [ ] **5.2** New title (2–3 options in Open questions).
  - Status / Notes:
- [ ] **5.3** Remove all super-resolution content.
  - Status / Notes:
- [ ] **5.4** Rewrite methods to match new code exactly.
  - Status / Notes:
- [ ] **5.5** Results: pipeline stats, multiple fit examples, synthetic recovery table, sensitivity tests, maps + quality layers, validation tables.
  - Status / Notes:
- [ ] **5.6** Limitations section (no XSM, relative units, footprint resolution, small n).
  - Status / Notes:
- [ ] **5.7** TGRS minor points: 24 archives + non-contiguity; rename "Flux Integration"; Fig. 1 caption/legend; 1024 vs 2048 channels; gain/calibration handling.
  - Status / Notes:
- [ ] **5.8** Highlights, cover letter draft, data availability statement (new Zenodo release → Open questions, no publish without approval).
  - Status / Notes:
- [ ] **5.9** Consistency check: script dumps headline numbers to JSON; cross-check against text.
  - Status / Notes:

## Phase 6 — Final review

- [ ] **6.1** Hostile-reviewer pass; write `RESPONSE_MAP.md`.
  - Status / Notes:
- [ ] **6.2** Compile PDF; check figures, references, template compliance.
  - Status / Notes:
- [ ] **6.3** Hand off to Geetesh: summary, final PDF, open questions.
  - Status / Notes:

---

## Decisions log

| Date | Task | Decision | Rationale |
|---|---|---|---|
| 2026-10-10 | 0 | New rework code goes in `src/asr/`, outputs in gitignored `outputs/asr/`; old `src/fitting/` left intact for baseline reproduction. | Keeps the old pipeline runnable for the 0.2 baseline and for reviewers' comparison. |

## Open questions for Geetesh

_(none yet)_
