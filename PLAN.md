# PLAN.md — Rework of the CLASS XRF lunar mapping paper for Advances in Space Research

## 0. First thing to do (before any other work)

1. Create `TRACKER.md` at the repo root, built from the task list in section 5 (checkbox format, grouped by phase, with a "Status / Notes" line per task and a "Decisions log" + "Open questions for Geetesh" section at the bottom).
2. Add this block to `CLAUDE.md` (create the file if it doesn't exist; keep any existing content):

   ```
   ## Working rules for the ASR rework
   - At the start of every session, read PLAN.md and TRACKER.md before doing anything else.
   - After finishing (or abandoning) any task, update TRACKER.md immediately: tick it, add a one-line note
     with the result (numbers, file paths, commit hash), and log any decision in the Decisions log.
   - Never mark a task done if its verification step hasn't passed.
   - Put anything that needs Geetesh's input under "Open questions" and keep working on unblocked tasks.
   ```
3. Create a branch `asr-rework` and work there.

## 1. Situation

This repo (SomaDhatuNirupana) holds the pipeline and paper "Global Fundamental-Parameter Mapping of Lunar Elemental Abundances from Chandrayaan-2 CLASS XRF Data with Cross-Modal Super-Resolution" (Kanjilal, Maji, Chandra; Geetesh Chandra is the corresponding author; all three listed as equal contribution).

It was **rejected by IEEE TGRS** (TGRS-2026-07721) and cannot go back there. The new target is **Advances in Space Research (Elsevier/COSPAR), free standard (non-open-access) route**. The aim is a defensible, honest paper submitted in about 3–4 weeks. This is not a full physics rebuild.

**Everything needed is local.** The uncommitted CLASS Level-1 FITS files from PRADAN (24 monthly archives, Dec 2020–May 2026, ~1.5M exposures) are in the working tree. **Never commit or delete them.** Add them to `.gitignore` if they aren't already there. Find the paper source (LaTeX), figure scripts, and the fitting code yourself. Summarise where everything lives in TRACKER.md as task 0.

## 2. What the reviewers said (summary — the rework must answer these)

**Physics / forward model**
- Eq. 7: primary intensity is just μᵢ(E_ℓ). There is no excitation integral over the solar spectrum, no photoionisation cross-section, and no matrix attenuation or geometry. The elemental constant K (jump × fluorescence yield × rate) only appears in the secondary term, so primary intensity ignores ω entirely. *(We found this one ourselves.)*
- Secondary fluorescence is ad hoc (product of μ's × β). The geometry is absorbed into the global scale s, which can't absorb energy-dependent effects.
- No solar-flux normalisation over 24 months. Ratios to Si do NOT cancel solar-spectrum hardness.

**Fitting**
- Scale degeneracy: s multiplies everything, and so do all Cᵢ. Only s·Cᵢ is identifiable.
- The Gaussian detector response is truncated at ±1σ, which discards ~32% of each line.
- Denoising (SWT + Savitzky–Golay + clipping) is done before an unweighted least-squares fit. It should be a Poisson likelihood on raw counts.
- No goodness-of-fit criterion, no uncertainties, and only one example fit shown.
- O is fitted at ~13.6 relative units because it absorbs the low-energy continuum.

**Data / mapping**
- The top-50%-flux selection may bias composition.
- Centre-point binning ignores the footprint. Observed vs smoothed/infilled coverage must be clearly separated, and grid spacing ≠ resolution.
- 1,056,619 footprints are used for coverage but only 855,190 are fitted. These need to be distinguished.

**Validation**
- Only Al/Si validates (r = 0.74, n = 9). Mg, Fe, and Ti are ~0. The 5° pooling radius is huge. There are no CIs and no leave-one-out.
- Needs comparison against an established orbital product.

**Super-resolution** — it is optical texture injection (R² ≈ 0.19), α = 1 is unjustified, and the held-out check shows no gain. **Decision: remove SR from the paper entirely.** Leave the code in the repo.

**Sanity check we found:** the current fitted ratios are physically implausible. Mg/Si is ~1–2.2, but lunar values are ~0.15–0.3 by weight. Ca/Si is ~0.01 vs ~0.35–0.55. Fe/Si is ~0.03 vs ~0.15–0.7. Only Al/Si (~0.3–0.8) is in range. This is consistent with missing excitation physics. Use these ranges as a smell test after every model change.

The full review PDF can be dropped into `reviews/` if more detail is needed.

## 3. Scope decisions (fixed — don't relitigate)

- **Do:** fix the forward model to a defensible simplified FP form, fix fitting and statistics, add uncertainties and quality layers, strengthen validation, cut SR, and rewrite with honest claims.
- **Don't:** build a full XSM time-resolved solar correction (see task 1.1 for the exception), add new ML models, or chase a bigger venue.
- **Honesty rule:** report whatever the fixed pipeline produces. Do not tune thresholds, radii, or parameters to recover the old r = 0.738 or to make more elements "validate". If results get worse, that goes in the paper and in TRACKER.md. Every analysis choice is fixed *before* looking at validation numbers, and logged.
- Zero cost: no paid services, no APC route.

## 4. Technical spec for the model fixes

- **Primary fluorescence (Sherman-type, simplified):**
  Iᵢ,ℓ ∝ Cᵢ · Kᵢ,ℓ · ∫_{E>edgeᵢ} S(E) · τᵢ(E) / [μ_matrix(E)/sinψ₁ + μ_matrix(E_ℓ)/sinψ₂] dE

  - τᵢ is the photoionisation cross-section, from xraylib (`CS_Photo` / shell-partial).
  - μ_matrix = Σⱼ Cⱼ μⱼ.
  - ψ₁ and ψ₂ are the incidence and emergence angles. Use the per-exposure solar incidence angle from the FITS header if available; take emergence as nadir.
- **Solar spectrum S(E):** check whether any Chandrayaan-2 XSM data exists locally.
  - If none exists (expected), use a fixed reference flare spectrum: an isothermal thin-thermal / bremsstrahlung shape at a nominal temperature, documented in the paper.
  - Run a sensitivity test with two or three temperatures, e.g. low-B, C, and M class. Report how much the ratios move. That number is the honest systematic for "no XSM correction".
- **Secondary fluorescence:** keep it, but make it concentration-coupled (depends on C_j of the exciter as well as Cᵢ). Alternatively, drop it and quantify its size on synthetic spectra. Pick whichever is defensible and log the decision.
- **Detector response:** a Gaussian integrated per channel bin (erf difference) over ±5σ, normalised. FWHM free, or as σ(E) = √(a + bE). Apply an energy-dependent detector efficiency if any CLASS response or efficiency info is in the FITS files or docs. Otherwise, note it as a limitation.
- **Identifiability:** fix C_Si ≡ 1, fit ratios rᵢ = Cᵢ/C_Si plus one amplitude. Drop O as a free parameter: either fix it stoichiometrically or exclude it. Fit range ~1.0–7.5 keV so the <1 keV continuum isn't absorbed by a fake element. Add a simple scattered-continuum/background term.
- **Statistics:** fit raw counts (no denoising) with a Cash / Poisson likelihood. The model is source + exposure-scaled nightside background, rather than subtracting the background. Keep denoising for plots only.
- **Per-fit outputs:**
  - ratio uncertainties from the inverse Hessian / Jacobian covariance
  - parameter correlation matrix and condition number
  - a goodness-of-fit statistic
  - a quality flag
- **Selection:** replace the top-50%-flux union with a pre-declared minimum-counts + goodness-of-fit criterion. Show a sensitivity test with 2–3 thresholds.

## 5. Task list (seed for TRACKER.md)

### Phase 0 — Orientation
- [ ] 0.1 Map the repo: data location, pipeline entry points, fitting module, mapping module, SR module, paper source, figure scripts. Write the summary into TRACKER.md.
- [ ] 0.2 Confirm the environment runs (one month, a few hundred spectra) and reproduces the old behaviour. Record the old ratios as a baseline.
- [ ] 0.3 Check FITS headers for incidence/emergence angles, exposure, gain, channel count, and any response or efficiency info. Check for any local XSM data.

### Phase 1 — Forward model and fitting
- [ ] 1.1 Implement the solar-spectrum model (reference flare spectrum, configurable temperature).
- [ ] 1.2 Rewrite primary fluorescence per section 4. Unit-test against hand-computed values for one element.
- [ ] 1.3 Fix secondary fluorescence (coupled) or drop it with justification.
- [ ] 1.4 Detector response: binned Gaussian over ±5σ, normalised. Add a unit test that the line area is conserved.
- [ ] 1.5 Reparameterise for identifiability (Si = 1, ratios + amplitude, O fixed or excluded, 1.0–7.5 keV range, continuum term).
- [ ] 1.6 Implement the Poisson/Cash fit on raw counts with the background model, plus covariance, uncertainties, GoF, and quality flag.
- [ ] 1.7 **Gate — synthetic recovery test:** generate synthetic CLASS-like spectra (realistic counts, Poisson noise, background, a range of lunar compositions and incidence angles). Fit them and report the bias and scatter of each ratio, and the condition numbers. Also generate with one solar temperature and fit with another to get the systematic. **Do not proceed to full runs until recovery is acceptable for at least Al/Si and Mg/Si.** Log the results.
- [ ] 1.8 Smoke test on one real month. Check the ratios against the lunar sanity ranges in section 2. Make diagnostic plots of high-, median-, and low-quality fits.

### Phase 2 — Full run (ask Geetesh before starting; it took ~8 h with 14 workers last time)
- [ ] 2.1 Apply the new selection criterion. Record counts at each stage (raw → dayside → passes QC → fitted).
- [ ] 2.2 Run the full fit. Save per-spectrum results with uncertainties and flags.
- [ ] 2.3 Run the selection-threshold sensitivity (subset is fine) and the solar-temperature sensitivity.

### Phase 3 — Maps
- [ ] 3.1 Rebuild the element/Si maps. Report **observed** coverage (directly populated cells) separately from smoothed/infilled coverage, and separately for fitted vs all footprints. Note the latitude-dependent pixel area, or use an equal-area grid for the statistics.
- [ ] 3.2 Produce quality layers per map: observation count, mean fit uncertainty, infill mask.
- [ ] 3.3 State the effective resolution honestly (footprint scale, not grid spacing).

### Phase 4 — Validation
- [ ] 4.1 Apollo/Luna (9 sites):
  - Add a table of reference compositions with sources.
  - Use radius sensitivity (e.g. 1°, 2°, 5°), fixed in advance.
  - Report bootstrap 95% CIs on r, leave-one-site-out, and propagated per-site uncertainties.
- [ ] 4.2 **Lunar Prospector GRS comparison** (public PDS products, e.g. Prettyman et al. 2006 elemental maps). If there's no network access, add an open question asking Geetesh to download them.
  - Regrid both datasets to a common coarse grid (e.g. 5°).
  - Compare Fe, Ti, Al, Ca, and Mg (where available) with the corresponding X/Si ratio.
  - Use Spearman + Pearson with spatially blocked uncertainty (block bootstrap).
  - Make scatter plots.
- [ ] 4.3 Decide which elements the paper claims, based only on the 4.1 + 4.2 evidence. Log the decision.

### Phase 5 — Paper rewrite (Elsevier `elsarticle` template for ASR)
- [ ] 5.1 Check the current ASR Guide for Authors (abstract limit, highlights, keywords, reference style, data statement, any length guidance). Record the requirements in TRACKER.md.
- [ ] 5.2 New title, without "fundamental-parameter" unless section 4 physics justifies it, and without "abundances". Something like *"Whole-spectrum fitting of Chandrayaan-2 CLASS XRF data for global relative element/Si mapping of the Moon"*. Offer 2–3 options in Open questions.
- [ ] 5.3 Remove all super-resolution content (sections, figures, tables, supplementary figures, abstract and conclusion claims).
- [ ] 5.4 Rewrite the methods to match the new code exactly. Equations must match the implementation.
- [ ] 5.5 Write the results:
  - pipeline statistics
  - multiple fit examples
  - the synthetic recovery table
  - sensitivity tests
  - maps with quality layers
  - validation tables
- [ ] 5.6 Write the limitations section honestly: no time-resolved XSM correction (with the quantified systematic), relative not absolute units, footprint-scale resolution, small n for landing sites.
- [ ] 5.7 Fix the TGRS minor points that still apply:
  - list the 24 monthly archives and explain the non-contiguity
  - rename the "Flux Integration" heading
  - fix the Fig. 1 caption/legend mismatch
  - explain the 1024 vs 2048 channels
  - document gain/calibration handling
- [ ] 5.8 Write the highlights (Elsevier format), a cover letter draft, and update the data availability statement. A new Zenodo release (e.g. v2.0.0) should go under Open questions — **don't publish or push it without Geetesh's approval.**
- [ ] 5.9 Run a consistency check: every number in the text matches the generated outputs. Write a script that dumps all headline numbers to a JSON file, and cross-check it against the text.

### Phase 6 — Final review
- [ ] 6.1 Do an internal "hostile reviewer" pass. Re-read the TGRS reviews point by point and write `RESPONSE_MAP.md`, mapping each concern → what changed / where it's addressed / or why it's out of scope.
- [ ] 6.2 Compile the PDF. Check the figures, references, and template compliance.
- [ ] 6.3 Hand off to Geetesh: summary of changes, final PDF, and the open questions list.

## 6. Ground rules

- Small commits on `asr-rework` with clear messages. No force-pushes, no pushing to `main`, no changes to the Zenodo record without approval.
- Never commit the FITS data or large intermediate outputs. Keep outputs in a gitignored `outputs/` directory.
- Ask before anything that takes over ~1 hour of compute.
- Test changes on a small subset before scaling up.
- When unsure about physics, pick the simpler defensible option, state it as an assumption in the paper, and log it.
