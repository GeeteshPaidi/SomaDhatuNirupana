# TRACKER.md — ASR rework progress

Branch: `asr-rework`. Plan: [PLAN.md](PLAN.md). Rules: [CLAUDE.md](CLAUDE.md).
New rework code lives in `src/asr/`; rework outputs go to `outputs/asr/` (gitignored).

---

## Phase 0 — Orientation

- [x] **0.1** Map the repo (summary below).
  - Status / Notes: done — see "0.1 Repo map". Phase 0 committed in 712dcaa (setup) + 08bd8a4 (results).
- [x] **0.2** Confirm the environment runs on one month / a few hundred spectra and reproduces the old behaviour; record old ratios as baseline.
  - Status / Notes: done (`src/asr/baseline_check.py` → `outputs/asr/baseline/baseline_summary.json`, `baseline_ratio_percentiles.csv`, `baseline_validation_expanded_5deg.csv`). Python 3.12.10, xraylib 4.2.1.
    - 2025-10: 27,237 files, 16,137 dayside, 1,406 nightside bg. Re-preprocessing 300 random dayside files: max relative diff vs stored catalogue = **0.0**. Re-fitting them: all 300 succeed; for the 146 also in the stored result, max |ΔC| = **0.0** for all 7 elements → bit-exact reproduction.
    - Old global baseline (855,190 successful fits), median [5–95%] X/Si: **Mg 1.41 [0.008–3.08]**, **Al 0.49 [0.004–2.97]**, **Ca 0.0093 [0.0005–0.20]**, **Fe 0.025 [0.002–0.27]**, **Ti 0.0054 [0.0001–0.15]**, **O 136 [5.7–731]**. 0.85% of fits have C_Si < 1e-3; 1.3–3.1% of Al/Mg/Ca/Ti fits pinned at the 0 bound. Confirms the PLAN §2 smell test (Mg too high, Ca/Fe/Ti ~10–50× too low, O absorbing continuum).
    - Old validation re-run on the 855,190-row expanded result (5°, median): Al r = **0.738**, Ca 0.263, Fe −0.002, Mg 0.006, Ti −0.117 → identical to `outputs/figures/validation_stats.csv`, so the paper's number did come from the expanded result (despite `validate.py`'s default path pointing at the 34,748-row `catalogue_full_result`).
- [x] **0.3** Check FITS headers (angles, exposure, gain, channels, response/efficiency info); check for local XSM data.
  - Status / Notes: done (`src/asr/survey_headers.py`, 200 files/month × 24 = 4,800 sampled; `outputs/asr/header_survey_summary.txt`). Total L1 files = **1,515,239** (matches paper). Findings:
    - GAIN = 13.5 eV/ch and DETCHANS = 2048 in **all** sampled files → the 1024-ch / 27.7 eV branch in the old code is dead; no per-file gain drift recorded, no channel offset keyword.
    - Incidence: `SOLARANG` (= `PHASEANG` in 100% of files, i.e. phase angle with nadir boresight). Emergence: `EMISNANG` ≈ 0 (max 2.8e-6°) → nadir confirmed. Use ψ₁ = 90° − SOLARANG, ψ₂ = 90°.
    - EXPOSURE = 8 s in 99.4%; a few partial exposures (1.8–7.8 s) → must use per-file exposure.
    - `SCD_USED`: 16 SCDs in ~74% of samples, but 26% use 1–15 SCDs (2022-03 median 12). **Effective area scales with n_SCD** → background model must be scaled by n_SCD (old pipeline subtracted a single 16-SCD-mixed nightside mean). Ratios unaffected, amplitude/QC thresholds are.
    - `TEMP` −54 to −22 °C (median −41) → gain/resolution may vary; no per-file calibration provided. Fit FWHM per spectrum and note as limitation.
    - `SW_VERSN` 2.0 (2020-12 … 2021-02) vs 2.1 (rest) → note in paper.
    - `RESPFILE`/`ANCRFILE` empty, `HDUCLAS2=TOTAL` (background not subtracted), `POISSERR=T`. No RMF/ARF or efficiency info anywhere in the archives (zips contain only .fits + .xml).
    - Some headers have SOLARANG = 0 and SAT_ALT = 0 (likely missing geometry) → QC must reject SAT_ALT ≤ 0.
    - **No XSM data locally** (searched repo + archive listing) → fixed reference solar spectrum per PLAN §4.

### 0.1 Repo map

**Data (never commit / delete)**
- `data/archive/ch2_cla_l1_YYYY_MM.zip` — 24 PRADAN monthly archives (14 GB): 2020-12, 2021-01/02, 2022-03, 2024-03, 2024-11 … 2026-05 (contiguous from 2024-11). Each zip = `cla/data/calibrated/YYYY/MM/DD/*.fits` + PDS4 `.xml` labels only (no XSM, no RMF/ARF).
- `data/raw/YYYY/MM/DD/ch2_cla_l1_<start>_<end>.fits` (+ `.xml`) — extracted L1 spectra, one 8 s exposure per file, 2048-ch OGIP PHA table (`CHANNEL`, `COUNTS`), all 16 SCDs summed.
- `data/processed/months/cat_YYYY_MM.parquet` — old per-month preprocessed catalogues (bg-subtracted, **denoised** count *rates*, float32) + `flux_<El>` window sums; `coverage_hit.npy`.
- `data/processed/catalogue_expanded_result.parquet` — old fit result used by the paper (855,190 rows = union top-50% of 1,056,619 dayside). Other `catalogue_full*` files are earlier, smaller runs (34,748 / 97,244 / 142,083 spectra).
- `data/reference/validation.xlsx` — Apollo/Luna site wt% (rows 34–39, Ti/Al/Fe/Mg/Ca/Si; sources not recorded in the sheet → needs a citation table in 4.1).
- `assets/moon_real.jpg` — optical guide image for SR only.

**Old pipeline (`src/fitting/`, kept intact as baseline)**
- Entry point: `run_pipeline.py` (month loop: `prepare_data_parallel.prepare_month` → union top-50% flux thresholds → `run_fit.fit_single_row` per spectrum → `catalogue_expanded_result.parquet`). Older entry points: `prepare_data.py` (serial), `run_full_fit.py`, `batch_run.py`, `build_fitset.py`, `enrich_result.py`.
- Preprocessing: `prepare_data_parallel.py` — counts/EXPOSURE → minus mean nightside rate (SOLARANG ≥ 90°, first 3000 files of the month) → clip ≥0 → SWT sym8 L2 + Savitzky–Golay(7,3) → clip ≥0. Energy scale hard-coded `0.0135 keV/ch × channel` (no offset).
- Forward model: `element_model.py` (`ElementModel`: xraylib line energies/RadRate, Kα1+Kβ1 only, `CS_Total`, `JumpFactor`, `FluorYield`) + `element_handler.py` (`ElementHandler.calculate_folded_intensity`): line intensity = `C_i·[μ_i(E_ℓ) + Σ μ_j(E_j)·μ_i(E_j)·K_i·β]`, Gaussian evaluated only within ±1σ window → confirms reviewer points (primary = μ only, no ω in primary, ±1σ truncation). Note `std_dev` param is σ but bounds/initial (0.1 keV) are not tied to the detector.
- Fitter: `xrf_fitter.py` (`XRFFitter`): 7 free C (Fe, Al, Mg, Si, Ca, Ti, O; bounds 0–200) + scale s (1e-8–1) + σ (0.01–0.5); unweighted TRF `least_squares` on denoised rates, full 0–27.6 keV range. `run_fit.py` wraps it per spectrum.
- 1024-channel branch: `energy_factor = 0.0277` if `num_channels != 2048` (never triggered — all L1 files are 2048 ch, GAIN = 13.5 eV/ch).

**Mapping (`src/mapping/`)**
- `build_maps_ipynb.py` — the paper's maps: X/Si ratio, log1p+IQR tail cut, **centre-point binning** at 1080×540, gaussian_filter smoothing + 3×3 infill.
- `build_maps.py` + `map_making_v1.py` (`GaussianArray`, footprint-polygon rasterisation, adaptive Gaussian) — described in an older version of Sec. 5; not what made the final maps.

**SR (to be cut from paper, code stays)**: `src/superres.py`, `src/sr_site_validate.py`, `src/figures/fig3_sr_sidebyside.py`, `src/figures/appendix_sr_maps.py`, `outputs/figures/sr_*`.

**Validation**: `src/validate.py` — 9 sites, centre-point within 5° (cos-lat scaled), median X/Si vs reference X/Si (wt%), Pearson r only. Output `outputs/figures/validation_{table,stats}.csv` (Al r = 0.738, n = 9; Fe −0.002, Mg 0.006, Ca 0.263, Ti −0.117).

**Figures (`src/figures/`)**: `fig1_preprocessing.py` (reads `outputs/figures/preprocessing_example.npz`), `fig2_fitted_spectrum.py`, `fig3_sr_sidebyside.py`, `fig_coverage_map.py` (1,056,619 dayside footprints, 99.0% 1° coverage), `appendix_sr_maps.py`; Fig. 4 is written by `src/validate.py`. Tracked outputs in `outputs/figures/final/`, `outputs/figures/supplementary/`.

**Paper**: `paper/main.tex` (IEEEtran, 1241 lines; SR = Sec. "Cross-Modal Guided Super-Resolution…" + results/ablation subsections), `paper/supplementary.tex`, compiled `main.pdf`/`supplementary.pdf`; untracked `paper/main_manuscript_source.zip`, `outputs/cover_letter.{tex,pdf}` (TGRS submission artefacts, left untouched). Reviews not in repo (`reviews/` is gitignored if added).

**New rework code**: `src/asr/` (header survey, baseline check, and the new model/fitter as they are written).

## Phase 1 — Forward model and fitting

- [x] **1.1** Solar-spectrum model (reference flare spectrum, configurable temperature).
  - Status / Notes: `src/asr/xrf_model.py::solar_spectrum` — isothermal bremsstrahlung, S(E) ∝ exp(−E/kT)/E, Gaunt = 1, continuum only (no flare lines). Reference kT = **1.07 keV** (single-T fit to a C3 flare measured by Ch-1 XSM, Narendranath et al. 2013 as quoted in arXiv:1305.6218). Sensitivity temperatures 0.6 keV (~7 MK, B-class) and 1.5 keV (~17 MK, M-class). Line-intensity ratios vs kT 0.6→1.5 keV (45°, mare comp.): Mg/Si 0.46→0.24, Al/Si 0.43→0.30, Ca/Si 0.028→0.29, Ti/Si 0.0008→0.022, Fe/Si 0.0004→0.092 → Ca/Ti/Fe are dominated by the assumed solar shape.
- [x] **1.2** Rewrite primary fluorescence per PLAN §4; unit-test against hand-computed values for one element.
  - Status / Notes: `ForwardModel.line_intensities` / vectorised `primary_fast` (Sherman integral with CS_Photo_Partial(K), ω_K, RadRate, matrix μ from oxide closure, ψ₁ = 90° − SOLARANG, ψ₂ = 90°). Kα + Kβ for Mg, Al, Si, Ca, Ti, Fe (Mg Kβ absent in xraylib). Tests (`tests/test_asr_model.py`, 13 pass): Si Kα vs independent scipy.quad Sherman integral within 1%; Q·τ_K vs xraylib CS_FluorLine within 5%; closure gives O 40–46 wt%; intensity monotonic in incidence. Grid 1.30–15 keV, ΔE = 5 eV → 0.3% integration error vs 1 eV grid; 0.08 ms/eval.
- [x] **1.3** Secondary fluorescence: concentration-coupled, or drop with justification.
  - Status / Notes: implemented the coupled thick-target formula (de Boer 1990) in `line_intensities(secondary=True)`; test confirms coupling (Ca secondary grows with Fe/Si). Size: sec/prim ≤ **1.3%** for all Kα lines over kT 0.6–1.5 keV and mare/highland compositions (Mg 0.5–1.2%, Al 0.8–1.2%, Ti ≤ 1.3%, Ca ≤ 0.6%, Si ≤ 0.2%). **Dropped from the fit** (30× slower); synthetic spectra are generated *with* it, so the recovery test measures the omission bias directly (see 1.7).
- [x] **1.4** Detector response: binned Gaussian over ±5σ, normalised; unit test for line-area conservation.
  - Status / Notes: `Detector` — erf-integrated per channel, clipped at ±5σ, renormalised; σ(E)² = σ₀² + F·w·E (F = 0.115, w = 3.68 eV). Tests: area = 1 to 1e-9 and centroid exact for 4 line energies × 2 σ₀. Efficiency: unity (no CLASS QE/filter data in the release or found in literature search) → limitation + open question. Calibration found on 2024-12 dayside stack (20,881 spectra): lines sit ~2 ch below `0.0135·ch`; linear fit to Mg/Si/Ca centroids gives **E ≈ 0.01363·ch + 0.017 keV**, FWHM ≈ 112 eV (Mg) – 157 eV (Si). Per-month calibration to be fitted in 1.8. Also: nightside background contains a strong instrumental **Al Kα** line → background must be a scaled template, not a fixed subtraction.
- [x] **1.5** Reparameterise for identifiability (Si ≡ 1, ratios + amplitude, O fixed/excluded, 1.0–7.5 keV, continuum term).
  - Status / Notes: `src/asr/fitter.py` — θ = (A, Mg/Si, Al/Si, Ca/Si, Ti/Si, Fe/Si, B, β): C_Si from oxide closure (fixed minor Na₂O 0.5, K₂O 0.1, MnO 0.15, Cr₂O₃ 0.2 wt%), O stoichiometric (not free). Continuum = physical scattered-solar term (Rayleigh+Compton / matrix μ) with free amplitude B. Fit window **1.15–7.5 keV** (avoids Na Kα 1.04 keV and the steep <1.1 keV excess seen in the stack).
- [x] **1.6** Poisson/Cash fit on raw counts with background model; covariance, uncertainties, GoF, quality flag.
  - Status / Notes: raw counts, no denoising; model = source + β × nightside-template × exposure × n_SCD/16. Cash statistic minimised by TRF on signed deviance residuals. Uncertainties = inverse Poisson Fisher matrix; correlation matrix + condition number stored. GoF z = (C − E[C])/√Var[C] with exact per-channel Poisson expectations (Kaastra 2017). Quality bits: 1 no-conv, 2 counts < 300, 4 |z| > 3, 8 ratio at bound, 16 cond > 1e4, 32 bad geometry, 64 covariance failure. ~60 ms/fit.
- [x] **1.7** GATE — synthetic recovery test (bias/scatter per ratio, condition numbers, cross-temperature systematic). Do not proceed until Al/Si and Mg/Si recover acceptably.
  - Status / Notes: **GATE PASSED for Mg/Si and Al/Si** (`src/asr/synth_recovery.py`; outputs `outputs/asr/synth/summary_gen*_fit*.csv`). Generator: forward model *with* secondary fluorescence, real 2024-12 nightside background template, 8 s / 16 SCD, ψ₁ = 10–80°, random lunar ratios, scatter level matched to the 2024-12 stack, Poisson noise. Fit: primary-only model.
    - Same kT (1.07/1.07), 200 spectra per level — median bias / robust scatter / pull width, by Si Kα counts 100 · 300 · 1k · 3k · 10k:
      - Mg/Si: bias +2.5%, +1.1%, −0.1%, 0.0%, +0.6% · scatter 35%, 18%, 7.7%, 4.2%, 2.3% · pulls 0.92–1.03
      - Al/Si: bias +4.9%, +0.4%, −0.4%, +0.9%, +0.8% · scatter 36%, 13%, 5.7%, 3.3%, 1.7% · pulls 0.85–0.97
      - Ca/Si: bias −5.4%, −1.1%, −0.6%, +0.9%, +0.3% · scatter 31%, 20%, 10%, 5.4%, 2.9%
      - Ti/Si: scatter 138%, 58%, 27%, 14%, 7.6%; Fe/Si: scatter 149%, 97%, 32%, 19%, 9.1% (unconstrained below ~1k Si counts)
      - Converged 100%; median cond(corr) 115 → 25; median GoF z ≈ −0.3; quality = 0 for 49% / 85% / 96% / 99% / 99.5%. 60 ms/fit.
    - Solar-shape systematic (fit kT 1.07 keV; 100 spectra each at 1k and 10k counts), median ratio bias:
      - true kT 0.6 keV: Mg **+40%**, Al **+24%**, Ca **−86%**, Ti −99%, Fe −99%. At 10k counts GoF z = 5.9 → wrong kT is detectable in bright spectra.
      - true kT 1.5 keV: Mg +3%, Al 0%, Ca **+79%**, Ti +140%, Fe **+520%**; GoF does not flag it.
    - Conclusion: Mg/Si and Al/Si carry a −3…+40% solar-shape systematic. Ca/Si is uncertain by about a factor of 2, and Ti/Si and Fe/Si by an order of magnitude, without XSM. This is the honest systematic for the paper.
- [ ] **1.8** Smoke test on one real month; check sanity ranges; diagnostic plots (high/median/low-quality fits).
  - Status / Notes: **BLOCKED on a scope decision (Open questions 3–4).** Pipeline runs (`src/asr/run_month.py`, `src/asr/diag_plots.py`; 2024-12, 3,000 random dayside spectra, ~60 s with 12 workers). Findings:
    - **Calibration (2024-12 stack, 75,707 spectra):** gain 13.714 eV/ch, offset +2.6 eV, σ₀ 73.5 eV → FWHM 183 eV @ 1.49 keV, 209 eV @ 5.9 keV (within the pre-launch "<250 eV @ 1.48 keV" spec). The gain is 1.6% above the header GAIN = 13.5.
    - **Solar-driven instrumental Al Kα.** Nightside spectra show a strong Al Kα line, plus continuum below 1.2 keV and a 6.7 keV bump, *only* when the spacecraft is sunlit. Median Al-line rate: 31.8 c/s (n = 8,555) when SOLARANG < shadow boundary (≈109° at 100 km), vs 0.16 c/s (n = 2,254) in the lunar umbra. It correlates with solar proxies (r = 0.61 with 0.6–0.9 keV counts, 0.58 with 6.7 keV) and not with particle background (r = 0.11 with 8–20 keV). It changes on less than 30-min scales, so orbit-interpolation is not viable. Its size is comparable to the lunar Si line. The old pipeline subtracted a month-mean of it, so old Al/Si is contaminated by an amount that varies with solar activity.
    - **New background model:** particle template P from umbra nightside (scaled per spectrum by a single 8–20 keV summed bin) plus instrumental template T from sunlit nightside (sunlit − P) with free amplitude γ. Synthetic test with γ_true = 0.5 (`summary_gen1.07_fit1.07_inst.csv`): γ is identifiable. Mg/Si and Al/Si biases stay ≤ 1%; Al/Si scatter at 300 / 1k / 3k / 10k Si counts = 31% / 14% / 5.2% / 3.5% (vs 13% / 5.7% / 3.3% / 1.7% without γ).
    - **Real-data ratios fail the sanity ranges at kT = 1.07 keV:** median Mg/Si 0.69, Al/Si 1.07, Ca/Si 0.011 (lunar 0.15–0.3, 0.3–0.8, 0.35–0.55). γ is mostly fitted at 0. This is the same pattern as the synthetic "true sun softer than assumed" case.
    - **kT scan on the same 3,000 spectra** (median ratio; ΣC):
      - kT 0.4: Mg 0.33, Al 0.71, Ca 0.47; ΣC 1,558,407
      - kT 0.6: Mg 0.50, Al 0.86, Ca 0.084; ΣC 1,528,197
      - kT 1.07: Mg 0.69, Al 1.07, Ca 0.011; ΣC 1,519,619
      - kT 1.5: Mg 0.81, Al 1.18, Ca 0.004; ΣC 1,522,719
    - Likelihood barely discriminates kT (and mildly prefers 1.07). Only kT ≈ 0.4 keV gives lunar-like ratios, and picking it for that reason would be tuning (honesty rule). Bright spectra (>3k counts) have median GoF z ≈ 4.7–5 at every kT, i.e. a residual misfit that kT does not fix (candidates: solar line emission missing from S(E), the instrumental-template shape on the dayside, efficiency).
    - Old pipeline note: Mg/Si high and Ca/Fe/Ti tiny in the TGRS results is the same solar-shape signature.

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
| 2026-10-10 | 0.3 | Geometry: ψ₁ = 90° − SOLARANG (incidence elevation), ψ₂ = 90° (nadir, EMISNANG ≈ 0). Background scaled per spectrum by exposure × n_SCD/16. QC rejects SAT_ALT ≤ 0 and SOLARANG ≥ 90°. | Header survey; n_SCD varies in 26% of spectra. |
| 2026-10-10 | 1.1 | Reference solar spectrum: isothermal bremsstrahlung kT = 1.07 keV (Ch-1 XSM C3 flare); sensitivity at 0.6 and 1.5 keV. Fixed before seeing any validation numbers. | Only citable flare temperature found; no XSM data locally. |
| 2026-10-10 | 1.3 | Secondary fluorescence dropped from the fit; kept in the synthetic generator to quantify omission bias. | ≤1.3% of primary for every Kα line; 30× cost. |
| 2026-10-10 | 1.4 | Detector efficiency = 1 (no data); per-month gain/offset/σ₀ from stacked dayside spectra, fixed in per-spectrum fits. | Single 8 s spectra cannot constrain calibration; stack shows ~2-ch offset from nominal. |
| 2026-10-10 | 1.5 | Fit window 1.15–7.5 keV; Na and K lines not modelled (matrix only). O stoichiometric via oxide closure. | Avoids Na Kα and the <1.1 keV excess; Na/K lines are minor. |
| 2026-10-10 | 0 | New rework code goes in `src/asr/`, outputs in gitignored `outputs/asr/`; old `src/fitting/` left intact for baseline reproduction. | Keeps the old pipeline runnable for the 0.2 baseline and for reviewers' comparison. |

## Open questions for Geetesh

1. **CLASS detector efficiency / filter.** No QE curve, filter thickness, or depletion depth ships with the L1 release, and a literature search found none (C1XS-era papers mention an Al filter; the nightside background shows a strong Al Kα line). Efficiency is set to 1, which biases Fe/Si and Ca/Si by an unknown multiplicative factor (Mg/Al/Si much less). If you have the CLASS calibration paper or an ARF (e.g. Pillai et al. 2021 / Narendranath et al. 2024), point me to it.
2. **Per-spectrum solar temperature (optional).** *(superseded by Q3)* Synthetic tests show GoF z ≈ 6 when bright spectra are fitted with the wrong kT. Fitting kT freely for bright spectra (Si counts ≳ 3k) could reduce the solar-shape systematic. This is not a full XSM correction, but it goes beyond the PLAN's "fixed reference spectrum". Do you want me to test it as a sensitivity analysis only?
3. **Solar spectrum source (blocking 1.8 → Phase 2).** Real ratios depend on the assumed kT by a factor of 1.5–2.5 for Mg/Si and Al/Si, and by about 100× for Ca/Si. The data cannot pin kT per spectrum. Options:
   - (a) **GOES XRS** 1–8 Å / 0.5–4 Å ratio → isothermal T(t) and EM(t) per exposure. Public NOAA data, zero cost, standard method (White et al. 2005). This gives a time-resolved solar shape without XSM. *Recommended.*
   - (b) Chandrayaan-2 **XSM** L1/L2 spectra from PRADAN (needs your PRADAN login). This is the proper answer, but the PLAN puts a full XSM correction out of scope.
   - (c) Keep the fixed reference spectrum (kT = 1.07 keV) and report Mg/Si and Al/Si as relative maps with a −3…+40% (or larger) solar systematic, with absolute ratios explicitly not interpretable.
4. **Instrumental Al (affects Al/Si, the only element that validated before).** Should the γ-template approach above be the method, or should Al/Si be restricted to spectra/time periods where the sunlit-nightside Al level is low?
