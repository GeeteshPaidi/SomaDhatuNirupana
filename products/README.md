# Data products (revised analysis, pipeline v5)

Products of the revised analysis (ancillary response applied at photon energy, detector-temperature-binned calibration, 8–17 keV particle proxy band, instrumental amplitude fixed to 0 on the dayside). All ratios are elemental weight ratios (oxide-closure model, Si reference). `products/v5/` has the same layout as `outputs/asr/`, where the pipeline writes them.

The products of the first submission (`*_v3*`) are in release v2.0.3 (Zenodo concept DOI 10.5281/zenodo.21993207, which resolves to the latest version).

| Path | Content |
|---|---|
| `v5/maps/mapping_set_v5.parquet` | 638,117 quality-0 CLASS footprints of the primary mapping set (24 monthly archives). Columns: `source_file` (L1 FITS name), `month`, `clat`/`clon` (footprint centre, deg), `solarang` (deg), `goes_T_mk` (GOES isothermal temperature, MK), `goes_xrsb` (W m⁻²), `n_counts_fit` (counts in 1.15–7.5 keV), `r_Mg`, `r_Al` (Mg/Si, Al/Si) and `sd_r_Mg`, `sd_r_Al` (1σ Fisher uncertainties). |
| `v5/maps/mapping_set_v5_eB.parquet` | Quality-0 footprints of the second epoch (February 2022 – December 2023). The first epoch (November 2024 – May 2026) is the subset of the primary set with those months; `src/asr/reproduce_release.py` derives it. |
| `v5/maps/maps_2deg_v5*.npz` | 2° equal-angle maps (90 × 180, row 0 = −90° lat, column 0 = −180° lon): `{Mg,Al}_mean` (weighted mean with footprint sharing), `_se` (formal standard error, √(Σf²w)/Σfw), `_wstd`, `_n` (effective footprint count Σf), `hit_all` / `hit_used`. Tags: `_v5` primary; `_v5_eA`, `_v5_eB` epochs; `_v5_b3` GOES ≥ B3; `_v5_min1000`, `_v5_min3000` count thresholds. No smoothing or infill. |
| `v5/maps/summary_2deg_v5*.json` | Pipeline stage counts and coverage statistics for each map set. |
| `v5/validation/` | Landing-site (`sites_*`) and Lunar Prospector GRS (`grs_*`) comparisons, claim-rule outcomes (`claims_*`), split-half reliabilities (`score_diagnostics_v5.json`), weighted-median check (`alt_median_posthoc_v5.json`); `*_v3fix*` are the first-stage γ ≡ 0 results used as the reference for the free-γ model. |
| `v5/epochs/epoch_compare.json` | Epoch comparison (2°, 6° and 10° blocks). |
| `v5/xsm/` | XSM comparison: XSM 1T/2T fits per minute (`xsm_fits_*`), CLASS refits with GOES, XSM 1T and XSM 2T spectra (`class_refits_*`), monthly summaries. |
| `v5/months_v5/<month>/` | Per-month calibration (`calib.json`: gain, offset, line width per detector-temperature bin), particle-background and instrumental templates (`background.npy`, `instrumental.npy`) and the stacked spectra used for calibration (`stack.npz`). |
| `figures_frozen/fig_background.png` | Figure 1, frozen from the first submission; its right-panel input was not retained. |
| `headline_numbers.json` | Every number quoted in the manuscript (source of `paper_asr/numbers.tex`). |

Not included: the per-spectrum fit tables of all dayside spectra (`outputs/asr/months_v5/<month>/fits.parquet`, about 1.1 GB, including flagged fits), which are available from the authors on request.

To regenerate Table 5 and Figures 1, 4 and 7 of the paper from a clean checkout:

    python src/asr/reproduce_release.py
