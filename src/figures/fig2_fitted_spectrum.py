"""
Figure 2: measured spectrum vs. best-fit FP forward model, with residuals.
Picks a representative well-fit spectrum from the top-flux subset actually
used for mapping/SR (flux_total >= 75th percentile -- the same population
the paper's maps are built from, per src/mapping/build_maps_ipynb.py's
flux selection), then within that subset takes the one closest to median
TRF cost (typical convergence, not cherry-picked for a flattering shape).
Selecting from the whole dataset (including the many low-count, low-SNR
exposures) tends to land on a near-featureless spectrum, since cost and
total flux are strongly correlated (rho=0.88) -- median cost overall is
mostly picking a median-FLUX (i.e. low-flux) exposure, not a median-FIT-
QUALITY one at representative signal levels.
"""
import glob

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from xrf_fitter import XRFFitter

# Use the expanded fit result (cost + source_file); recover the spectrum
# from the corresponding per-month catalogue by matching source_file.
result = pd.read_parquet("../../data/processed/catalogue_expanded_result.parquet")
result = result[result["success"] == True]

flux_thresh = result["flux_total"].quantile(0.75)
high_flux = result[result["flux_total"] >= flux_thresh]

median_cost = high_flux["cost"].median()
idx = (high_flux["cost"] - median_cost).abs().idxmin()
row = high_flux.loc[idx]
src = row["source_file"]

# source_file looks like ch2_cla_l1_YYYYMMDDT..._...fits -> month YYYY/MM
stamp = src.split("_")[3]           # YYYYMMDDT......
year, month = stamp[:4], stamp[4:6]
cat_path = f"../../data/processed/months/cat_{year}_{month}.parquet"
mcat = pd.read_parquet(cat_path, columns=["counts", "source_file"])
counts = mcat.loc[mcat["source_file"] == src, "counts"].iloc[0]

fitter = XRFFitter(counts=counts, std_dev=0.1)
params = [row[f"conc_{el}"] for el in fitter.elements] + [row["scale"], row["std_dev"]]
model = fitter.calculate_model(params)
residual = counts - model

plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white",
                     "savefig.facecolor": "white", "font.size": 9})
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(3.5, 3.4), dpi=600, sharex=True,
                                gridspec_kw={"height_ratios": [3, 1]})

ax1.plot(fitter.energies, counts, "k.", markersize=1.5, label="Measured spectrum")
ax1.plot(fitter.energies, model, "r-", linewidth=1.0, label="FP forward-model fit")
ax1.set_ylabel("Count rate (counts s$^{-1}$)", fontsize=8.5)
ax1.tick_params(labelsize=8)
ax1.legend(fontsize=7, frameon=True, edgecolor="0.7", handlelength=1.8)
for side in ("top", "right"):
    ax1.spines[side].set_visible(False)

ax2.plot(fitter.energies, residual, color="0.4", linewidth=0.6)
ax2.axhline(0, color="black", linewidth=0.5)
ax2.set_xlabel("Energy (keV)", fontsize=8.5)
ax2.set_ylabel("Residual", fontsize=8.5)
ax2.tick_params(labelsize=8)
for side in ("top", "right"):
    ax2.spines[side].set_visible(False)

# CLASS operates in the 0.8-15 keV band (Sec 2.1); none of the 7 modeled
# elements have lines above ~7 keV (Fe Kbeta), so higher energies are
# background/noise outside the model's physical scope, not a fit failure.
# Zooming to 0-8 keV (rather than 0-16) keeps the visible line structure
# from being compressed by a long, empty high-energy tail.
ax1.set_xlim(0, 8)

plt.tight_layout()
import os
os.makedirs("../../outputs/figures/final", exist_ok=True)
plt.savefig("../../outputs/figures/final/fig2_fitted_spectrum.pdf", bbox_inches="tight")
plt.savefig("../../outputs/figures/final/fig2_fitted_spectrum.png", dpi=600, bbox_inches="tight")

print("Converged concentrations for caption:")
for el in fitter.elements:
    print(f"  {el}: {row[f'conc_{el}']:.5f}")
print(f"  scale: {row['scale']:.5e}, std_dev: {row['std_dev']:.4f}, cost: {row['cost']:.3f}")
print("Saved fig2_fitted_spectrum.png/.pdf")
