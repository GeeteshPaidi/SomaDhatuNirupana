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
    - **Addendum (final configuration: GOES-T CHIANTI spectrum, two-component background, Ti/Fe fixed)**, 80 spectra/level at Si counts 300 / 1k / 3k / 10k:
      - T_gen = T_fit = 5 MK: Mg bias −1 / −4 / −2 / −1%, scatter 19 / 10 / 7 / 4%; Al bias −3 / 0 / 0 / +0.5%, scatter 29 / 15 / 7 / 4%; pulls 0.85–1.19. Ca bias +30% even at the right T, scatter 58–136% → **Ca/Si not recoverable per spectrum.**
      - GOES T error +30% (true 6.5, fit 5): Mg −13…−16%, Al +9…+16%, Ca ×9. Error −20% (true 4, fit 5): Mg +32…+36%, Al +1…+6%, Ca → 0.
      - → Claimed products limited to **Mg/Si and Al/Si**. Systematic from the GOES-T uncertainty: Mg/Si −15/+35%, Al/Si 0/+15% for ±20–30% T error.
    - (Original continuum-model) Conclusion: Mg/Si and Al/Si carry a −3…+40% solar-shape systematic. Ca/Si is uncertain by about a factor of 2, and Ti/Si and Fe/Si by an order of magnitude, without XSM. This is the honest systematic for the paper.
- [x] **1.8** Smoke test on one real month; check sanity ranges; diagnostic plots (high/median/low-quality fits).
  - Status / Notes: done after the Q3/Q4 decisions. 2024-12, 3,000 random dayside spectra (`outputs/asr/months/2024_12/`, `diag_fits.png`).
    - **Calibration** (stack of 75,707): gain 13.714 eV/ch, offset +2.6 eV, σ₀ 73.5 eV (FWHM 183 eV @ 1.49 keV, 209 eV @ 5.9 keV).
    - **Solar-driven instrumental Al Kα found.** It is present in nightside spectra only when the spacecraft is sunlit (31.8 vs 0.16 c/s), tracks solar proxies and not particles, and varies on <30 min. Handled by a two-component background: particle template P (umbra nightside, β pinned by a single 8–20 keV summed bin) + instrumental template T (sunlit − P, free γ). It cannot be fixed from GOES: Spearman 0.39, ×2.2 scatter, and the scale depends on SOLARANG (`outputs/asr/inst_vs_goes_2024_12.csv`).
    - **Solar spectrum:** CHIANTI isothermal lines + continuum (`sunkit_spex`, coronal abundances) at the per-exposure GOES XRS temperature (`src/asr/solar.py`; GOES-16 for 2020–22, GOES-18 for 2024–26; 1-min science data, White et al. 2005 via `sunkit_instruments`). Yearly median T: 2020 2.8 MK, 2021 8.0, 2022 6.0, 2024 5.5, 2025 5.1, 2026 5.3. Solar lines carry 55–80% of the exciting flux above the Mg/Si edges at 3–15 MK, so the continuum-only model was badly wrong.
    - Continuum-only kT scan (superseded): Mg/Si 0.33–0.81, Ca/Si 0.47–0.004 for kT 0.4–1.5 keV; ΣC nearly flat.
    - **GOES + CHIANTI result** (Ti/Si, Fe/Si fixed): quality = 0 for 47%. Median Mg/Si **0.30** [5–95%: 0.12–0.54], bright spectra 0.235 → inside the lunar range. Al/Si **1.12** [0.39–1.61], median rel. unc. 67% → high; instrumental-Al degeneracy plus unknown efficiency. Ca/Si unconstrained in most spectra (rel. unc. ~300%).
    - Final-config synthetic gate (`summary_gen*_fit5.0_final.csv`, T = 5 MK, γ = 0.5, 80/level): see 1.7 addendum.
    - Fit speed 96 ms/spectrum single-core (tol 1e-6; identical results to 1e-8).
## Phase 2 — Full run (ask Geetesh before starting; ~8 h with 14 workers last time)

- [x] **2.1** Apply new selection criterion; record counts at each stage (raw → dayside → QC → fitted).
  - Status / Notes: `outputs/asr/maps/summary_2deg.json`: raw L1 **1,515,239** → dayside with valid geometry **1,053,989** (all fitted, no flux selection) → GOES T reliable (≥ B1) **1,017,226** → converged **1,037,841** → quality-0 mapping set **429,217** (40.7% of dayside; was 440,234 before the boundary-flag fix). Per-month quality-0 fraction ranges 2% (2021-02) to 66% (2025-02); see `run_all_log.txt`.
- [x] **2.2** Full fit; save per-spectrum results with uncertainties and flags.
  - Status / Notes: done, v2 (with filter): all 24 months OK (`outputs/asr/months/<YYYY_MM>/fits.parquet`, calib.json with per-month gain/offset/σ₀ and template sources). ~5.5 h wall, 12 workers. v1 (no filter) log archived.
- [x] **2.3** Selection-threshold sensitivity (subset is fine) and solar-temperature sensitivity.
  - Status / Notes: count thresholds done (`summary_2deg_min*.json`, `validation/*_min*.csv`). Mapping set 429,217 / 122,552 / 10,672 for ≥300 / 1000 / 3000 counts. GRS r: Mg 0.60 / 0.59 / 0.62, Al 0.35 / 0.38 / 0.40 (n pixels 1790 / 1533 / 134). Median Mg/Si 0.231 / 0.250 / 0.247, Al/Si 0.816 / 0.785 / 0.855. B3 GOES cut: 428,915 spectra, results unchanged. Solar-T ×0.8 / ×1.2 refits (2024-12, 2025-06, 2025-12; 10k each, both-quality-0 pairs): median Mg/Si change −24% / +26%, Al/Si −8% / +7%.
## Phase 3 — Maps

- [x] **3.1** Rebuild element/Si maps; observed vs smoothed/infilled coverage; fitted vs all footprints; equal-area stats.
  - Status / Notes: `src/asr/maps.py` → `outputs/asr/maps/maps_2deg.npz`, `map_{Mg,Al}_Si_2deg.png`. Observed coverage (no smoothing or infill anywhere): all dayside footprints 100% of 2° cells / 100% of area; quality-0 set 99.99% cells / 100% area; ≥5 footprints in cells covering 99.0% of area. Footprint rasterised from 4 corners (5×5 samples).
- [x] **3.2** Quality layers per map: observation count, mean fit uncertainty, infill mask.
  - Status / Notes: count, SE of weighted mean, weighted scatter per cell (npz + figure panels). No infill is used, so the infill mask is trivially empty (stated in paper).
- [x] **3.3** State effective resolution honestly (footprint scale, not grid spacing).
  - Status / Notes: footprint corners span median 25.3 km (across track) × 37.2 km (along track; FOV + 8 s motion); 5–95%: 20.6–31.0 × 32.0–42.7 km. Collimator FWHM 12.5 km (Pillai et al. 2021). Maps at 2° (61 km). In paper §Mapping.
## Phase 4 — Validation

- [x] **4.1** Apollo/Luna (9 sites): reference table with sources; pre-fixed radius sensitivity (1°, 2°, 5°); bootstrap 95% CIs; leave-one-site-out; propagated per-site uncertainties.
  - Status / Notes: `src/asr/validate_asr.py` → `outputs/asr/validation/sites_{table,stats}.csv`. 2° (primary): Mg/Si r = 0.37 [−0.22, 0.88], LOO 0.16–0.50; Al/Si r = 0.43 [−0.08, 0.96], Spearman 0.68 (p = 0.04), LOO 0.30–0.71; RMSE Mg 0.061, Al 0.183; bias Mg +0.046, Al +0.100. 1°: Mg −0.68, Al 0.23. 5°: Mg 0.00, Al −0.09. Per-site SEs are tiny (0.005–0.014) compared with the site-to-site spread, so the scatter is systematic (mixing, footprint ≫ landing-site scale, regional heterogeneity), not statistical. **Neither passes** the pre-registered CI criterion. v1 (pre-fix) numbers archived (`*_v1_prereg.csv`): 2° Mg 0.38, Al 0.44.
- [x] **4.2** Lunar Prospector GRS comparison (common 5° grid; Fe, Ti, Al, Ca, Mg vs X/Si; Spearman + Pearson; block bootstrap; scatter plots).
  - Status / Notes: PDS LP-L-GRS-5-ELEM-ABUNDANCE-V1.0 5° equal-area (1,790 pixels, all with ≥20 CLASS footprints). Block bootstrap 20° blocks. **Mg/Si r = 0.60 [0.54, 0.65], ρ = 0.61. Al/Si r = 0.35 [0.26, 0.44], ρ = 0.37** (v1 before the boundary fix: 0.21 [0.13, 0.28]). Fe/Ti/Ca not compared: not products (fixed or unconstrained). `validation_scatter.png`.
- [x] **4.3** Decide which elements the paper claims, based only on 4.1 + 4.2. Log decision.
  - Status / Notes: per the pre-registered rule, **no element is "claimed" as validated**: both pass GRS, both fail the 9-site CI test. Paper framing: Mg/Si and Al/Si are presented as relative maps whose large-scale pattern agrees with GRS (Mg strongly, Al moderately). The landing-site test is reported as inconclusive (n = 9, power-limited), without changing radius or thresholds. Ca/Ti/Fe explicitly not retrievable without XSM.
## Phase 5 — Paper rewrite (Elsevier `elsarticle`, ASR)

- [x] **5.1** Check current ASR Guide for Authors; record requirements here.
  - Status / Notes: Elsevier/ScienceDirect guide pages return 403 to automated fetch. Known/confirmed: Elsevier `elsarticle` (installed locally), free subscription route (APC only for open access), ASR 2025 IF 3.2 / CiteScore 5.6 (COSPAR news), ~51% acceptance (SciRev), ~13 printed pages per paper (2022 issue: 825 pp / 64 papers). Assumed standard Elsevier items: highlights 3–5 × ≤ 85 characters (ours: 5, max 71), keywords, data availability statement, competing-interest declaration, author-year references (elsarticle-harv). **To verify in Editorial Manager at submission:** abstract word limit (ours ≈ 230 words), graphical abstract (optional?), review model.
- [x] **5.2** New title (2–3 options in Open questions).
  - Status / Notes: working title in `paper_asr/main.tex`; options in Open question 6.
- [x] **5.3** Remove all super-resolution content.
  - Status / Notes: new manuscript `paper_asr/` contains no SR content (sections, figures, tables, abstract, conclusions). Old `paper/` left untouched as the TGRS record.
- [x] **5.4** Rewrite methods to match new code exactly.
  - Status / Notes: `sections/{data,model,fitting,synthetic,mapping}.tex`. Equations correspond to `xrf_model.py` (closure, primary, scatter, response), `fitter.py` (model counts, Cash, Fisher, GoF, flags) and `maps.py`.
- [x] **5.5** Results: pipeline stats, multiple fit examples, synthetic recovery table, sensitivity tests, maps + quality layers, validation tables.
  - Status / Notes: written (`results.tex`, `validation.tex`); all numbers via `numbers.tex`. Final synthetic gate (with filter, `summary_gen*_fit5.0_final.csv`): Mg/Al |bias| ≤ 3%, pulls 0.78–1.21; T error +30% → Mg −22…−14%, Al +7…+11%; −20% → Mg +31…+35%, Al 0…+4%; Ca/Si bias up to +51% with 52–156% scatter even at the right T.
- [x] **5.6** Limitations section (no XSM, relative units, footprint resolution, small n).
  - Status / Notes: `discussion.tex` §Limitations (6 items, incl. SCD QE, flare exclusion by GoF noted in results).
- [x] **5.7** TGRS minor points: 24 archives + non-contiguity; rename "Flux Integration"; Fig. 1 caption/legend; 1024 vs 2048 channels; gain/calibration handling.
  - Status / Notes: archives and non-contiguity in data.tex §CLASS; "Flux Integration" and the old Fig. 1 no longer exist (no denoising stage); 1024 vs 2048 explained (all files 2048, 1024 branch unused); gain fitted per month (Table months, model.tex §Detector).
- [x] **5.8** Highlights, cover letter draft, data availability statement (new Zenodo release → Open questions, no publish without approval).
  - Status / Notes: highlights.tex, `paper_asr/cover_letter.tex`, data_availability.tex now cites the Zenodo concept DOI 10.5281/zenodo.21993207. v2.0.0 released 2026-10-11 (merged to master, GitHub release → Zenodo) on Geetesh's request; `products/` holds the released mapping set, maps and validation tables.
- [x] **5.9** Consistency check: script dumps headline numbers to JSON; cross-check against text.
  - Status / Notes: `src/asr/make_numbers.py` → `outputs/asr/headline_numbers.json` + `paper_asr/numbers.tex`; text uses only these macros, so values cannot drift. Verified: no "??" in numbers.tex or the compiled PDF (commit b34a894).
## Phase 6 — Final review

- [x] **6.1** Hostile-reviewer pass; write `RESPONSE_MAP.md`.
  - Status / Notes: `RESPONSE_MAP.md` maps every TGRS concern (P1–3, F1–5, D1–3, V1–2, S1) plus the self-found issues.
- [ ] **6.2** Compile PDF; check figures, references, template compliance.
  - Status / Notes: `paper_asr/main.pdf` compiles cleanly (27 pp elsarticle preprint, no undefined refs). Remaining: verify the 4 bib entries marked VERIFY (de Boer 1990, Narendranath 2014, 2020, Mithun 2020); confirm the abstract limit in Editorial Manager.
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
| 2026-10-10 | 1.8 / Q3 | **Geetesh approved option (a): per-exposure solar temperature from GOES XRS** (1–8 Å / 0.5–4 Å ratio → isothermal T, EM). Fixed kT = 1.07 keV retained only as a sensitivity case. | Real ratios depend strongly on kT; GOES is free and public; not an XSM correction. |
| 2026-10-10 | 1.8 / Q4 | **Instrumental Al: two-component background (particle umbra template P + solar-driven instrumental template T with free γ)** adopted as the method. | Recommended option; synthetic recovery shows γ identifiable with ≤1% bias. |
| 2026-10-10 | 1.8 | Ti/Si = 0.04 and Fe/Si = 0.37 fixed (global lunar-soil means, matrix only); conditioning computed on the core block (A, Mg/Si, Al/Si, B, β, γ). | At GOES T 3–8 MK the solar flux above the Ti/Fe edges is too small; free Ti/Fe made 92% of fits ill-conditioned. |
| 2026-10-10 | 1.8 | Ca/Si is fitted but **not claimed** per spectrum. | Final-config synthetic: +30% bias at the right T, ×9 / →0 under ±20–30% T error. |
| 2026-10-10 | 2 | **Pre-declared (before any validation):** mapping set = quality == 0 (Si-window counts ≥ 300, abs(GoF z) ≤ 3, Mg/Si & Al/Si not at bound, core cond ≤ 1e4, valid geometry, GOES T available). Threshold sensitivity: min counts 300 / 1000 / 3000. Solar-T sensitivity: GOES T × 0.8 and × 1.2 on a 3-month subset. | PLAN §4 selection; fixed before 4.1/4.2. |
| 2026-10-10 | 3–4 | **Pre-declared:** maps on a 2° equal-angle grid, plus 5° for the GRS comparison. Per-cell statistic = inverse-variance-weighted mean of footprint ratios, footprints rasterised by their 4 corners. Observed coverage = cells hit by ≥1 quality-0 footprint; no infill for statistics. Validation radii 1° / 2° / 5° (primary = 2°), bootstrap 95% CI (10,000), leave-one-site-out, block bootstrap with 20° blocks for GRS. Element claimed only if Pearson r > 0 with 95% CI excluding 0 in BOTH Apollo/Luna (primary radius) and GRS. | Honesty rule: fix analysis choices before seeing numbers. |
| 2026-10-10 | 2 | **Pre-declared:** GOES temperature reliable only if XRS-B (1–8 Å) ≥ 1e-7 W m⁻² (B1). Spectra below → quality bit 256 (excluded from maps). Sensitivity: B3 (3e-7) cut. | GOES-16 T below B1 is pinned at the 2.81 MK table floor (A-class) or diverges (A3–B1, 95th pct 22 MK, max 2.6e5 MK); see `goes_tem.parquet` by class. |
| 2026-10-10 | 4.1 | Reference compositions in `validation.xlsx` identified as average soil compositions (Lunar Sourcebook, Heiken et al. 1991, Table 7.15) — to be confirmed by Geetesh. | Values match the Sourcebook table; sheet has no source column. |
| 2026-10-10 | 2 | Months with < 100 umbra (or sunlit) nightside exposures borrow the 2024-12 particle (or instrumental) template shape; recorded per month in calib.json (`*_template_source`). | 2022-03 has 0 umbra exposures (empty template broke every fit); β, γ are fitted, so only shapes are borrowed. |
| 2026-10-10 | 2 | GOES-16 used to fill GOES-18 gaps in 2024–2025 (24,694 minutes). Cross-check over 1–15 Jun 2024 (20,145 minutes ≥ B1): median T16/T18 = 1.005 (16–84%: 0.991–1.030), XRS-B ratio 1.002. | GOES-18 good-minute coverage only 80% in 2024-03 and 85% in 2024-04. |
| 2026-10-10 | 1.4 / 2 | **Detector efficiency = published CLASS optical-blocking filter transmission**: 2 layers × (0.2 µm polyimide + 0.2 µm Al) = 0.4 µm Al + 0.4 µm polyimide (Narendranath et al. 2014; Pillai et al. 2021). Applied at photon energy (lines, scatter fine grid). Transmission: Mg Kα 0.858, Al Kα 0.908, Si Kα 0.691, Ca Kα 0.949. SCD QE is still not modelled (no published curve). Full run v1 (no filter) stopped at 2024-11 and restarted; v1 log kept as `run_all_log_v1_nofilter.txt`. | Published instrument parameter, found via literature search after Q1. Not tuned: without it Al/Si is inflated ×1.31 and Mg/Si ×1.24. Smoke test after the change: Mg/Si median 0.23, Al/Si 0.83. |
| 2026-10-10 | 1.8 | Instrumental Al vs SOLARANG (`paper_asr/figures/fig_background`): ~0.5 c/s at 90–93°, 30–100 c/s at 95–109°, drops in umbra → needs the Sun behind the spacecraft; expected small on the dayside, consistent with dayside γ ≈ 0. γ stays free. | Data. |
| 2026-10-10 | 3/4 | **Post-validation bug fix (disclosed):** the at-bound threshold for Mg/Si, Al/Si raised from 1e-6 to 1e-3 (fitter + maps.py). 254 Al/Si fits pinned at ~1.5e-5 had collapsed Fisher SDs (<1e-3) and carried 93% of all Al/Si inverse-variance weight. Pre-fix results archived as `*_v1_prereg.*` and reported in the paper alongside the fixed ones. | Defect in uncertainty validity at a boundary (Fisher invalid), visible in the map/GRS scatter as a spurious Al/Si ≈ 0 population; not a parameter tuned toward agreement. Mg unaffected (top 0.1% fits carry 2.7% of the weight). |
| 2026-10-10 | 4.3 | Claim decision per pre-registered rule: neither element meets "both tests"; paper reports GRS agreement and inconclusive site test without changing choices. | Honesty rule. |
| 2026-10-10 | 4 (post hoc) | **Score diagnostics** (`src/asr/score_diagnostics.py` → `validation/score_diagnostics.json`). GRS itself vs the 9 sites: Mg −0.11, Al −0.34 → the site test is non-diagnostic. CLASS split-half reliability (alternate months, Spearman–Brown): Mg 0.82, Al 0.42; GRS reliability from its error matrix: 0.73 / 0.78; ceilings 0.77 / 0.57 (observed 0.60 / 0.35). Old pipeline vs GRS on the same pixels: Mg 0.48, Al 0.19. | Requested by Geetesh: why are the scores low. |
| 2026-10-10 | 4 (post hoc) | **Weighted median** (`alt_median_map.py`): GRS Mg 0.58, Al 0.31; Al reliability 0.34 → no improvement. **γ ≈ 0 subset** (γ < 0.05; 295,814 = 69%): GRS Al 0.72 [0.66, 0.77], Mg 0.56 [0.49, 0.61]; Al reliability 0.96; sites 2° Al 0.80 [0.55, 0.97] (1°: 0.75, 5°: 0.43), Mg 0.15. Diagnosis: monthly farside Al/Si vs fraction γ>0, Spearman −0.88. All reported in the paper as post hoc, next to the pre-registered results. | Instrumental-Al degeneracy, not fluorescence physics, limits Al/Si. Next step (pre-specifiable): fix γ = 0 on the dayside. |
| 2026-10-10 | 5.2 | Title (b): *Global Mg/Si and Al/Si maps of the Moon from Chandrayaan-2 CLASS with time-resolved solar excitation from GOES XRS*. Abstract trimmed to 246 words. | Geetesh's choice. |
| 2026-10-10 | v3 (PRE-DECLARED before running) | **v3 model: instrumental amplitude γ fixed to 0 for all dayside spectra** (particle template + β only). Everything else is identical to v2: GOES/CHIANTI, filter, calibration, quality rules, 2° grid, radii, statistics, claim rule. Justification is physical: the nightside geometry (fig_background) shows the instrumental Al needs the Sun behind the spacecraft (SOLARANG > ~93°), i.e. it is absent on the dayside. **v3 becomes the primary product whatever its validation outcome**; v2 (free γ) is reported as the alternative. No further model changes after seeing v3 numbers. Outputs: `outputs/asr/months_v3/`. | Removes the post-hoc status of the γ≈0 Al/Si result. Approved by Geetesh 2026-10-10. |
| 2026-10-11 | v3 T-sens | v3 solar-T ×0.8 / ×1.2 (3 months × 10k): Mg/Si −25% / +27%, Al/Si −9% / +8%. Final manuscript compiles: 30 pp, no undefined refs, no '??'. | |
| 2026-10-11 | v3 result | v3 run complete (`outputs/asr/months_v3`, `run_all_log_v3.txt`). Mapping set 578,530. GRS: **Al/Si 0.72 [0.66, 0.76], Mg/Si 0.61 [0.56, 0.66]**. Sites 2°: **Al 0.81 [0.48, 0.97]** (1°: 0.71, 5°: 0.47), Mg 0.38 (fails). **Claim rule: Al/Si CLAIMED (both tests); Mg/Si not (site test).** Reliability Al 0.97, Mg 0.76; ceilings 0.87 / 0.75. Thresholds: GRS Al 0.72 / 0.72 / 0.69, Mg 0.61 / 0.58 / 0.58 (300 / 1k / 3k); B3 unchanged. Synthetic v3: abs(bias) ≤ 4%, pulls 0.84–1.25. Paper switched to v3 primary; v2 reported as the alternative model with full history. | Per the pre-declaration. |
| 2026-10-11 | 5.8 | Repo cleaned for release: products/ added (v3 mapping set 40 MB, 2° maps, validation tables, headline numbers); unused v2 weighted-median figure removed; TGRS submission artefacts gitignored; .zenodo.json / CITATION.cff → v2.0.0; asr-rework merged into master; release v2.0.0. | Geetesh asked to publish on Zenodo now; text revisions will follow as a new version. |
| 2026-10-10 | 0 | New rework code goes in `src/asr/`, outputs in gitignored `outputs/asr/`; old `src/fitting/` left intact for baseline reproduction. | Keeps the old pipeline runnable for the 0.2 baseline and for reviewers' comparison. |

## Open questions for Geetesh

1. **[PARTLY ANSWERED by literature: filter = 2×(0.2 µm polyimide + 0.2 µm Al), now in the model; SCD QE still unknown]** **CLASS detector efficiency / filter.** No QE curve, filter thickness, or depletion depth ships with the L1 release, and a literature search found none (C1XS-era papers mention an Al filter; the nightside background shows a strong Al Kα line). Efficiency is set to 1, which biases Fe/Si and Ca/Si by an unknown multiplicative factor (Mg/Al/Si much less). If you have the CLASS calibration paper or an ARF (e.g. Pillai et al. 2021 / Narendranath et al. 2024), point me to it.
2. **Per-spectrum solar temperature (optional).** *(superseded by Q3)* Synthetic tests show GoF z ≈ 6 when bright spectra are fitted with the wrong kT. Fitting kT freely for bright spectra (Si counts ≳ 3k) could reduce the solar-shape systematic. This is not a full XSM correction, but it goes beyond the PLAN's "fixed reference spectrum". Do you want me to test it as a sensitivity analysis only?
3. **[ANSWERED 2026-10-10: option (a) GOES]** **Solar spectrum source.** Real ratios depend on the assumed kT by a factor of 1.5–2.5 for Mg/Si and Al/Si, and by about 100× for Ca/Si. The data cannot pin kT per spectrum. Options:
   - (a) **GOES XRS** 1–8 Å / 0.5–4 Å ratio → isothermal T(t) and EM(t) per exposure. Public NOAA data, zero cost, standard method (White et al. 2005). This gives a time-resolved solar shape without XSM. *Recommended.*
   - (b) Chandrayaan-2 **XSM** L1/L2 spectra from PRADAN (needs your PRADAN login). This is the proper answer, but the PLAN puts a full XSM correction out of scope.
   - (c) Keep the fixed reference spectrum (kT = 1.07 keV) and report Mg/Si and Al/Si as relative maps with a −3…+40% (or larger) solar systematic, with absolute ratios explicitly not interpretable.
4. **[ANSWERED 2026-10-10: γ-template method]** **Instrumental Al (affects Al/Si, the only element that validated before).** Should the γ-template approach above be the method, or should Al/Si be restricted to spectra/time periods where the sunlit-nightside Al level is low?
5. **[ANSWERED: official Apollo/Luna returned-sample data; verified = Lunar Sourcebook Table 7.15 best-average soils (Basu & Riegsecker), Apollo rows identical to Basu & Riegsecker 1998, LPI abstract 6042; both cited]** **Reference table source.** `data/reference/validation.xlsx` has no citation. Its values match the Lunar Sourcebook (Heiken, Vaniman & French 1991), Table 7.15, average soil compositions. Please confirm, or tell me the real source.
6. **[ANSWERED: option (b)]** **Title (pick one).**
   - (a) *Whole-spectrum Poisson fitting of Chandrayaan-2 CLASS X-ray fluorescence spectra with GOES-constrained solar excitation: global Mg/Si and Al/Si maps of the Moon* (current).
   - (b) *Global Mg/Si and Al/Si maps of the Moon from Chandrayaan-2 CLASS with time-resolved solar excitation from GOES XRS.*
   - (c) *What can CLASS measure without a solar monitor? Poisson whole-spectrum fitting, instrumental aluminium and global Mg/Si and Al/Si maps.*
7. **[ANSWERED 2026-10-11: publish v2.0.0 on Zenodo now; paper text to be revised by Geetesh, then a new version]** **Zenodo v2.0.0 release + GitHub push of `asr-rework`.** Not done; needs your approval. The data availability statement currently says "archived on Zenodo upon acceptance".
8. **[ANSWERED: yes → `fig_maps_wmedian_posthoc`, labelled post-hoc]** **Al/Si map speckle.** A few isolated cells are pulled low by near-boundary fits with small uncertainties (inverse-variance mean, pre-declared). A weighted-median map would remove them, but it would be a post-validation change. Currently kept as is and disclosed only implicitly. Do you want a weighted-median version added as a clearly labelled sensitivity figure?

