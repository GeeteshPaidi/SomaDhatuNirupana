# Global Mg/Si and Al/Si maps of the Moon from Chandrayaan-2 CLASS

Processing chain, data products and manuscript source for:

> *Global Mg/Si and Al/Si maps of the Moon from Chandrayaan-2 CLASS with time-resolved solar excitation from GOES XRS*
> Smarak Kanjilal, Subarno Maji, Paidi Geetesh Chandra.

The chain fits every Chandrayaan-2 CLASS Level-1 X-ray fluorescence spectrum on its raw counts with a Poisson (Cash) likelihood. The forward model is a simplified fundamental-parameter model:

- Sherman primary fluorescence with matrix attenuation and incidence geometry;
- a scattered solar spectrum;
- an erf-binned Gaussian detector response;
- the published CLASS filter transmission;
- a particle-background template.

The incident solar spectrum is an isothermal CHIANTI spectrum at the GOES XRS temperature of each exposure. Only public data are used.

## Repository layout

```
src/asr/
  xrf_model.py            Forward model (composition closure, primary/secondary fluorescence,
                          CLASS filter transmission, detector response)
  fitter.py               Cash-statistic fitter: Fisher uncertainties, GoF, quality flags
  solar.py                GOES XRS -> T, EM per minute; CHIANTI spectrum table (sunkit-spex)
  run_month.py            One month: read FITS, background templates, energy calibration, fits
  run_all.py              All months (one subprocess per month)
  maps.py                 2-deg maps, quality layers, pipeline stage counts
  reference.py            Returned-soil reference compositions and landing-site coordinates
  validate_asr.py         Apollo/Luna sites (1/2/5 deg) and Lunar Prospector GRS (5 deg) validation
  score_diagnostics.py    Split-half reliabilities and correlation ceilings
  alt_median_map.py       Weighted-median sensitivity map
  synth_recovery.py       Synthetic recovery tests
  survey_headers.py       FITS header survey
  make_numbers.py         Every number in the paper -> numbers.tex / headline_numbers.json
  paper_figures.py, diag_plots.py   Figures
tests/                    Unit tests (pytest)
products/                 Mapping-set footprint ratios, 2-deg maps, validation tables
                          (see products/README.md)
paper_asr/                Manuscript source (elsarticle); numbers.tex is auto-generated
data/reference/           Returned-soil averages (Lunar Sourcebook Table 7.15)
```

## Input data (not redistributed)

| Data | Source | Local path |
|---|---|---|
| CLASS Level-1 spectra (24 monthly archives) | ISRO PRADAN, https://pradan.issdc.gov.in/ch2/ | `data/raw/YYYY/MM/DD/*.fits` |
| GOES XRS 1-min science data (GOES-16 2020–2025, GOES-18 2024–2026) | NOAA NCEI | `data/goes/sci_xrsf-l2-avg1m_g1?_yYYYY_v2-2-1.nc` |
| Lunar Prospector GRS elemental abundances | PDS Geosciences, LP-L-GRS-5-ELEM-ABUNDANCE-V1.0 | `data/grs/lpgrs_high1_elem_abundance_5deg.{tab,lbl}` |

## Installation

Python 3.12:

```bash
pip install -r requirements.txt
```

## Reproducing the results

```bash
python src/asr/solar.py --table --goes          # CHIANTI table + GOES temperatures
python src/asr/run_all.py --gamma_zero --out_root outputs/asr/months_v3   # all months (~5.5 h, 12 CPU workers)
python src/asr/maps.py --months_dir outputs/asr/months_v3 --tag _v3
python src/asr/validate_asr.py --mapping_set outputs/asr/maps/mapping_set_v3.parquet --tag _v3
python src/asr/synth_recovery.py --solar --kT_gen 5 --kT_fit 5 --n 80 \
       --levels 300,1000,3000,10000 --bkg outputs/asr/months/2024_12/background.npy \
       --inst outputs/asr/months/2024_12/instrumental.npy --tag _final
python src/asr/make_numbers.py                  # numbers.tex + tables
python src/asr/paper_figures.py && python src/asr/paper_figures.py results
cd paper_asr && pdflatex main && bibtex main && pdflatex main && pdflatex main
python -m pytest tests
```

Everything runs on CPU (NumPy/SciPy, multiprocessing).

## Citation

See [`CITATION.cff`](CITATION.cff). Zenodo: [10.5281/zenodo.21993207](https://doi.org/10.5281/zenodo.21993207) (concept DOI, resolves to the latest version).

## License

MIT — see [`LICENSE`](LICENSE).
