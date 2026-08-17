# Lunar Elemental Abundance Mapping from Chandrayaan-2 CLASS XRF Data

This repository implements a full, physics-based pipeline for turning raw Chandrayaan-2 CLASS
(Chandrayaan-2 Large Area Soft X-ray Spectrometer) X-ray fluorescence spectra into global maps of
lunar elemental abundance. It performs fundamental-parameter (FP) forward modelling and
trust-region-reflective (TRF) least-squares spectral fitting to recover per-footprint element
concentrations (Fe, Al, Mg, Si, Ca, Ti, O), projects those fits onto a spatial grid with
adaptive Gaussian weighting, and fuses the resulting low-resolution abundance maps with a
co-registered optical reflectance image via a cross-modal, FFT-based super-resolution step to
sharpen spatial detail. Fitted and mapped abundances are validated against Apollo/Luna landing-site
ground-truth geochemistry.

This repository accompanies a paper submitted to **IEEE Transactions on Geoscience and Remote
Sensing (TGRS)**:

> *Global Fundamental-Parameter Mapping of Lunar Elemental Abundances from Chandrayaan-2 CLASS
> XRF Data with Cross-Modal Super-Resolution*

The paper source and compiled PDF are in [`paper/`](paper/); a supplementary PDF with additional
per-element maps is at [`paper/supplementary.pdf`](paper/supplementary.pdf).

## Data source

Raw CLASS Level-1 XRF spectra are publicly available via the **PRADAN** (Policy-based Data
Retrieval, Analytics, Dissemination, and Notification) portal, hosted by the Indian Space Science
Data Centre (ISSDC), ISRO: **https://pradan.issdc.gov.in/ch2/**

This repository does not redistribute the raw FITS archive (see [Data layout](#data-layout) below)
— fetch it from PRADAN directly.

## Repository layout

```
├── paper/                     IEEE TGRS paper source, compiled PDF, supplementary material
├── data/
│   ├── raw/                    Raw CLASS L1 FITS spectra (not tracked — fetch from PRADAN)
│   ├── archive/                Original PRADAN zip downloads (not tracked)
│   ├── processed/              Preprocessed/fitted parquet catalogues (not tracked, regenerable)
│   └── reference/
│       └── validation.xlsx      Apollo/Luna landing-site ground-truth geochemistry
├── assets/
│   └── moon_real.jpg           Optical reflectance guide image used by the super-resolution step
├── src/
│   ├── fitting/                 FP forward model + TRF spectral fitting
│   ├── mapping/                 Spatial abundance mapping (adaptive Gaussian weighting)
│   ├── postprocess/             Utilities that merge/combine fit outputs
│   ├── figures/                 Scripts that render the paper's figures
│   ├── superres.py              Cross-modal guided super-resolution (Pearson correlation,
│   │                            ridge regression, FFT high-pass fusion, alpha ablation)
│   ├── validate.py              Ground-truth validation against Apollo/Luna sites
│   └── sr_site_validate.py      Held-out check of the super-resolved product itself
├── outputs/figures/
│   ├── final/                    The paper's main figures
│   ├── supplementary/            The paper's supplementary figures
│   └── *.csv                     Validation and ablation result tables
└── requirements.txt
```

## Installation

Requires Python 3.10+.

```bash
pip install -r requirements.txt
```

## Running the pipeline

All commands below assume you are in the relevant `src/` subdirectory (each script's default
paths are relative to its own location) and that `data/raw/` has been populated from PRADAN.

### 1. Preprocessing + fitting (`src/fitting/`)
Fits the fundamental-parameter forward model to every raw spectrum via TRF least-squares,
recovering 7 elemental concentrations per footprint.

```bash
cd src/fitting
python run_pipeline.py --raw_dir ../../data/raw --out ../../data/processed/catalogue_result.parquet
```

Use `--prepare_only` to stop after preprocessing/coverage reporting, or `--skip_prepare` to reuse
already-preprocessed month catalogues. See `element_model.py` / `element_handler.py` for the FP
forward model itself, and `xrf_fitter.py` for the TRF fit.

### 2. Spatial mapping (`src/mapping/`)
Projects the per-footprint fits onto a lat/lon grid with adaptive Gaussian footprint weighting.

```bash
cd src/mapping
python build_maps.py --fit_result ../../data/processed/catalogue_result.parquet \
                      --out_dir ../../outputs/figures/maps --ratio
```

### 3. Cross-modal super-resolution (`src/`)
Fuses a mapped abundance array with the optical guide image (`assets/moon_real.jpg`) to sharpen
spatial detail, and computes the Pearson correlation / ridge-regression / ablation numbers
reported in the paper.

```bash
cd src
python superres.py --element Al --maps_dir ../outputs/figures/maps --optical ../assets/moon_real.jpg
```

### 4. Validation (`src/`)
Compares fitted/mapped abundances against Apollo/Luna landing-site ground truth
(`data/reference/validation.xlsx`), and checks whether the super-resolved product itself moves
estimates toward or away from ground truth (leakage-free, since no ground-truth data enters the
SR fusion step).

```bash
cd src
python validate.py --fit_result ../data/processed/catalogue_result.parquet
python sr_site_validate.py --sr_result ../outputs/figures/sr_result_Al.npz
```

### 5. Figures (`src/figures/`)
Each script renders one paper figure from the outputs of the stages above, saving to
`outputs/figures/final/` or `outputs/figures/supplementary/`.

```bash
cd src/figures
python fig1_preprocessing.py
python fig2_fitted_spectrum.py
python fig3_sr_sidebyside.py
python appendix_sr_maps.py
python fig_coverage_map.py
```

## Citation

If you use this code or the accompanying paper, please cite as described in
[`CITATION.cff`](CITATION.cff).
