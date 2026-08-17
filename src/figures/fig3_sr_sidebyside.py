"""
Figure 3: element/Si ratio map before/after cross-modal guided super-
resolution, plus an explicit residual (difference) panel implementing
Sec 6.5's Delta(x,y) = A_SR(x,y) - A_hat(x,y) directly as a figure
(previously only defined as an equation with no accompanying visual).

NOTE on units: fitted concentrations are relative scale factors, not wt%
(see README "Known gaps") - the colorbar is labeled accordingly rather than
mislabeled as wt%, per the explicit instruction to flag rather than
silently misrepresent this.
"""
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from scipy.ndimage import maximum_filter

ELEMENT = sys.argv[1] if len(sys.argv) > 1 else "Al"
NPZ_DIR = sys.argv[2] if len(sys.argv) > 2 else "../../outputs/figures/sr_expanded_lat60"

d = np.load(f"{NPZ_DIR}/sr_result_{ELEMENT}.npz")
low_res = d["abundance_interp"]
sr = d["A_sr"]
alpha = float(d["alpha"])
residual = sr - np.nan_to_num(low_res, nan=np.nanmean(low_res))


def find_top_residual_features(low_res, sr, n=3, min_separation_px=40):
    """Locate the n most prominent local maxima of |A_sr - low_res|: the
    actual pixels where guided SR injected the most new spatial detail,
    rather than picking illustrative coordinates by eye."""
    res_abs = np.abs(np.nan_to_num(sr, nan=0.0) - np.nan_to_num(low_res, nan=0.0))
    local_max = (res_abs == maximum_filter(res_abs, size=min_separation_px)) & (res_abs > 0)
    ys, xs = np.where(local_max)
    vals = res_abs[ys, xs]
    order = np.argsort(vals)[::-1]
    picked = []
    for i in order:
        y, x = ys[i], xs[i]
        if all((y - py) ** 2 + (x - px) ** 2 > min_separation_px ** 2 for py, px in picked):
            picked.append((y, x))
        if len(picked) == n:
            break
    return picked


h, w = sr.shape
feature_px = find_top_residual_features(low_res, sr, n=3)
features = []
for (y, x) in feature_px:
    lon = (x / w) * 360 - 180
    lat = (y / h) * 180 - 90
    features.append((lon, lat))

vmin = np.nanpercentile(low_res, 2)
vmax = np.nanpercentile(low_res, 98)
norm = Normalize(vmin=vmin, vmax=vmax)
cmap = plt.get_cmap("viridis")

res_lim = np.nanpercentile(np.abs(residual), 98)
res_norm = TwoSlopeNorm(vmin=-res_lim, vcenter=0.0, vmax=res_lim)
res_cmap = plt.get_cmap("RdBu_r")

plt.rcParams.update({"font.size": 7.5})
# Sized to the IEEEtran double-column text width (7.16in) since this figure
# must span both columns (figure* in LaTeX) to stay legible -- three
# side-by-side maps do not fit in a single 3.5in column at any readable
# font size.
fig, axes = plt.subplots(1, 3, figsize=(7.16, 2.35), dpi=600)
fig.subplots_adjust(wspace=0.32)
extent = [-180, 180, -90, 90]

for j, (ax, arr, title, (this_norm, this_cmap)) in enumerate(zip(
    axes,
    [low_res, sr, residual],
    [f"(a) {ELEMENT}/Si (pre-SR)",
     f"(b) {ELEMENT}/Si (super-resolved)",
     f"(c) Residual $\\Delta = A_{{SR}} - \\hat{{A}}$"],
    [(norm, cmap), (norm, cmap), (res_norm, res_cmap)],
)):
    im = ax.imshow(arr, origin="lower", extent=extent, cmap=this_cmap, norm=this_norm, aspect="equal")
    ax.set_xlabel("Longitude (deg)", fontsize=7.5)
    if j == 0:
        ax.set_ylabel("Latitude (deg)", fontsize=7.5)
    ax.set_xticks(range(-180, 181, 90))
    ax.set_yticks(range(-90, 91, 45))
    ax.tick_params(labelsize=7)
    ax.set_title(title, fontsize=8)

    # Mark the pixels with the largest |A_sr - low_res| residual: these are
    # where guided SR actually injected the most new spatial detail relative
    # to the input map, found data-driven (see find_top_residual_features),
    # not chosen by eye.
    for i, (lon, lat) in enumerate(features):
        lbl = chr(ord("A") + i)
        ax.plot(lon, lat, marker="o", markersize=5, markerfacecolor="none",
                markeredgecolor="black" if this_cmap is res_cmap else "red", markeredgewidth=1.0)
        ax.annotate(lbl, (lon, lat), color="black" if this_cmap is res_cmap else "red",
                    fontsize=7, fontweight="bold", xytext=(3, 3), textcoords="offset points")

cbar1 = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes[:2], orientation="horizontal",
                      fraction=0.06, pad=0.32, aspect=45)
cbar1.set_label(f"{ELEMENT}/Si ratio (relative units, see caveat in text)", fontsize=7.5)
cbar1.ax.tick_params(labelsize=7)

cbar2 = fig.colorbar(plt.cm.ScalarMappable(norm=res_norm, cmap=res_cmap), ax=axes[2], orientation="horizontal",
                      fraction=0.06, pad=0.32, aspect=20)
cbar2.set_label("Residual (SR $-$ input)", fontsize=7.5)
cbar2.ax.tick_params(labelsize=7)

import os
os.makedirs("../../outputs/figures/final", exist_ok=True)
plt.savefig(f"../../outputs/figures/final/fig3_sr_sidebyside.pdf", bbox_inches="tight")
plt.savefig(f"../../outputs/figures/final/fig3_sr_sidebyside.png", dpi=600, bbox_inches="tight")
print(f"Saved fig3_sr_sidebyside.pdf/.png to outputs/figures/final/ for {ELEMENT} (alpha={alpha})")
