# RESPONSE_MAP.md — TGRS-2026-07721 concerns → ASR manuscript

Each reviewer concern (as summarised in PLAN.md §2) is mapped to what changed, where it is addressed in `paper_asr/`, or why it is out of scope. Section names refer to `paper_asr/sections/*.tex`.

## Physics / forward model

| # | Concern | Response | Where |
|---|---|---|---|
| P1 | Eq. 7: primary intensity was just μᵢ(E_ℓ); no excitation integral, no photo cross-section, no matrix attenuation or geometry; ω only in the secondary term. | Replaced by the Sherman primary integral: CS_Photo_Partial(K) × ω_K × radiative rate, matrix μ from oxide closure, incidence ψ₁ = 90° − SOLARANG and nadir emergence. Unit-tested against an independent quadrature (1%) and against xraylib CS_FluorLine (5%). | model.tex §Primary fluorescence, Eq. (primary); tests/test_asr_model.py |
| P2 | Secondary fluorescence ad hoc (μ products × β); geometry absorbed in the global scale. | Implemented the concentration-coupled thick-target formula (de Boer 1990). It is ≤ 1.3% of primary for every Kα line, so it is omitted from the fit but kept in the synthetic generator; the omission bias is part of the measured recovery bias. β and the global geometric scale are removed. | model.tex §Secondary; synthetic.tex |
| P3 | No solar-flux normalisation over 24 months; ratios do not cancel solar hardness. | Agreed, and quantified. Per-exposure CHIANTI isothermal spectrum at the GOES XRS temperature (lines + continuum). Synthetic and real-data sensitivity show Mg/Si and Al/Si move by tens of % with T, and Ca/Ti/Fe by orders of magnitude, so those three are no longer reported. | data.tex §GOES; model.tex §Solar; discussion.tex |

## Fitting

| # | Concern | Response | Where |
|---|---|---|---|
| F1 | Scale degeneracy (s × all Cᵢ). | Reparameterised: ratios rᵢ = Cᵢ/C_Si, C_Si fixed by oxide closure, one amplitude A. O is stoichiometric, not free. | model.tex §Closure; fitting.tex |
| F2 | Gaussian truncated at ±1σ (~32% of each line lost). | Gaussian integrated per channel (erf), ±5σ, renormalised; area conserved to 1e-9 (unit test); σ(E) with the Fano term. | model.tex §Detector, Eq. (response) |
| F3 | Denoising before unweighted LSQ; should be Poisson on raw counts. | Raw counts, Cash statistic; background modelled, not subtracted. No denoising anywhere. | fitting.tex §Likelihood |
| F4 | No GoF, no uncertainties, one example fit. | Fisher covariances, correlation matrix and condition number, exact Poisson C-stat GoF z (Kaastra 2017), quality flags; three example fits (bright / median / faint); synthetic pull widths 0.85–1.19 confirm the uncertainties. | fitting.tex §Uncertainties; results.tex Fig. fits; synthetic.tex |
| F5 | O fitted at ~13.6 because it absorbs the low-energy continuum. | O is not fitted. The fit window starts at 1.15 keV, and the continuum is a physical scattered-solar term plus measured background templates. | fitting.tex; model.tex §Scattered |

## Data / mapping

| # | Concern | Response | Where |
|---|---|---|---|
| D1 | Top-50%-flux selection may bias composition. | Removed. All 1,053,989 dayside spectra fitted. The pre-declared selection uses counts ≥ 300, abs(GoF z) ≤ 3, no parameter at a bound, conditioning, and GOES ≥ B1. Threshold sensitivity 300 / 1000 / 3000 reported. | fitting.tex; validation.tex §GRS |
| D2 | Centre-point binning ignores footprint; observed vs smoothed/infilled coverage; grid ≠ resolution. | Footprint quadrilaterals rasterised (5×5 samples). No smoothing or infill anywhere. Coverage reported separately for all footprints and for the mapping set, in cell and area terms. Resolution stated as the footprint scale (25 × 37 km; collimator FWHM 12.5 km); maps at 2°. | mapping.tex; results.tex §Maps |
| D3 | 1,056,619 footprints for coverage vs 855,190 fitted. | Stage table: 1,515,239 files → 1,053,989 dayside (all fitted) → 1,017,226 GOES-reliable → 429,217 mapped. Coverage is quoted for both the dayside and mapping sets. | results.tex Table stages |

## Validation

| # | Concern | Response | Where |
|---|---|---|---|
| V1 | Only Al/Si validated (r = 0.74, n = 9); 5° radius huge; no CIs; no LOO. | Radii 1/2/5° pre-declared (2° primary), bootstrap 95% CIs, leave-one-out, per-site SEs. Result reported honestly: neither element passes the 9-site test (Mg r = 0.37, Al r = 0.43, CIs include 0). The old r = 0.74 was not reproduced and was not chased. | validation.tex §Landing sites, Table sites |
| V2 | Needs comparison with an established orbital product. | Lunar Prospector GRS (Prettyman et al. 2006), 1,790 equal-area 5° pixels, 20° block bootstrap: Mg/Si r = 0.60 [0.54, 0.65], Al/Si r = 0.35 [0.26, 0.44]; stable under threshold sensitivity. | validation.tex §GRS, Fig. validation |

## Super-resolution

| # | Concern | Response | Where |
|---|---|---|---|
| S1 | SR is optical texture injection; α = 1 unjustified; no held-out gain. | Removed from the paper entirely (code left in the repo, unused). | — |

## Issues found during the rework (not raised by reviewers, disclosed in the paper)

- **Solar-driven instrumental Al Kα** present only when the spacecraft is sunlit. The old mean-nightside subtraction removed a time-variable amount of Al. Now a separate background template with a free amplitude (model.tex §Background, Fig. background).
- **Optical-blocking filter transmission** (0.4 µm Al + 0.4 µm polyimide, published) applied. Without it, Al/Si and Mg/Si were inflated by 31% and 24% (model.tex §Detector).
- **In-flight gain** is 13.74 eV/ch (not the header's 13.5), fitted per month (Table months).
- **Boundary-flag bug** found after the first validation (254 pinned Al/Si fits carried 93% of the weight). Fixed; pre-fix numbers reported (validation.tex §Disclosed correction).

## Out of scope (stated as limitations)

- Simultaneous XSM solar spectra (the model accepts them; future work).
- SCD quantum efficiency / dead layer (unpublished).
- Absolute abundances; Ca/Si, Ti/Si, Fe/Si.
- Particle-size and roughness effects; phase-angle effects beyond plane-parallel geometry.
