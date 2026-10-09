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
| 2026-10-10 | 0.3 | Geometry: ψ₁ = 90° − SOLARANG (incidence elevation), ψ₂ = 90° (nadir, EMISNANG ≈ 0). Background scaled per spectrum by exposure × n_SCD/16. QC rejects SAT_ALT ≤ 0 and SOLARANG ≥ 90°. | Header survey; n_SCD varies in 26% of spectra. |
| 2026-10-10 | 0 | New rework code goes in `src/asr/`, outputs in gitignored `outputs/asr/`; old `src/fitting/` left intact for baseline reproduction. | Keeps the old pipeline runnable for the 0.2 baseline and for reviewers' comparison. |

## Open questions for Geetesh

_(none yet)_
