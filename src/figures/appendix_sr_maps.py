"""
Supplementary-material figures: same before/after/residual 3-panel format
as main-paper Fig. 3, for the non-headline elements (Mg, Ca, Fe, Ti) whose
correlation/regression/validation numbers are reported in main-text
Tables II/III/V but whose SR maps were moved out of the 10-page main
manuscript (see supplementary.tex).
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from scipy.ndimage import maximum_filter


def find_top_residual_features(low_res, sr, n=3, min_separation_px=40):
    residual = np.abs(np.nan_to_num(sr, nan=0.0) - np.nan_to_num(low_res, nan=0.0))
    local_max = (residual == maximum_filter(residual, size=min_separation_px)) & (residual > 0)
    ys, xs = np.where(local_max)
    vals = residual[ys, xs]
    order = np.argsort(vals)[::-1]
    picked = []
    for i in order:
        y, x = ys[i], xs[i]
        if all((y - py) ** 2 + (x - px) ** 2 > min_separation_px ** 2 for py, px in picked):
            picked.append((y, x))
        if len(picked) == n:
            break
    return picked


SR_DIR = "../../outputs/figures/sr_ipynb_lat60"


def make_figure(element):
    d = np.load(f"{SR_DIR}/sr_result_{element}.npz")
    low_res = d["abundance_interp"]
    sr = d["A_sr"]
    alpha = float(d["alpha"])
    residual = sr - np.nan_to_num(low_res, nan=np.nanmean(low_res))

    h, w = sr.shape
    feature_px = find_top_residual_features(low_res, sr, n=3)
    features = [((x / w) * 360 - 180, (y / h) * 180 - 90) for (y, x) in feature_px]

    vmin = np.nanpercentile(low_res, 2)
    vmax = np.nanpercentile(low_res, 98)
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap("viridis")

    res_lim = np.nanpercentile(np.abs(residual), 98)
    res_norm = TwoSlopeNorm(vmin=-res_lim, vcenter=0.0, vmax=res_lim)
    res_cmap = plt.get_cmap("RdBu_r")

    # Sized to the IEEEtran double-column text width (7.16in): this figure
    # spans both columns (figure* in LaTeX), same treatment as Fig. 3.
    plt.rcParams.update({"font.size": 7.5})
    fig, axes = plt.subplots(1, 3, figsize=(7.16, 2.35), dpi=600)
    fig.subplots_adjust(wspace=0.32)
    extent = [-180, 180, -90, 90]
    titles = [f"(a) {element}/Si (pre-SR)",
              f"(b) {element}/Si (super-resolved)",
              f"(c) Residual $\\Delta = A_{{SR}} - \\hat{{A}}$"]
    norms_cmaps = [(norm, cmap), (norm, cmap), (res_norm, res_cmap)]

    for j, (ax, arr, title, (this_norm, this_cmap)) in enumerate(
            zip(axes, [low_res, sr, residual], titles, norms_cmaps)):
        ax.imshow(arr, origin="lower", extent=extent, cmap=this_cmap, norm=this_norm, aspect="equal")
        ax.set_xlabel("Longitude (deg)", fontsize=7.5)
        if j == 0:
            ax.set_ylabel("Latitude (deg)", fontsize=7.5)
        ax.set_xticks(range(-180, 181, 90))
        ax.set_yticks(range(-90, 91, 45))
        ax.tick_params(labelsize=7)
        ax.set_title(title, fontsize=8)
        mcolor = "black" if this_cmap is res_cmap else "red"
        for i, (lon, lat) in enumerate(features):
            lbl = chr(ord("A") + i)
            ax.plot(lon, lat, marker="o", markersize=5, markerfacecolor="none",
                    markeredgecolor=mcolor, markeredgewidth=1.0)
            ax.annotate(lbl, (lon, lat), color=mcolor, fontsize=7, fontweight="bold",
                        xytext=(3, 3), textcoords="offset points")

    cbar1 = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes[:2], orientation="horizontal",
                          fraction=0.06, pad=0.32, aspect=45)
    cbar1.set_label(f"{element}/Si ratio (relative units)", fontsize=7.5)
    cbar1.ax.tick_params(labelsize=7)

    cbar2 = fig.colorbar(plt.cm.ScalarMappable(norm=res_norm, cmap=res_cmap), ax=axes[2], orientation="horizontal",
                          fraction=0.06, pad=0.32, aspect=20)
    cbar2.set_label("Residual (SR $-$ input)", fontsize=7.5)
    cbar2.ax.tick_params(labelsize=7)

    import os
    os.makedirs("../../outputs/figures/supplementary", exist_ok=True)
    plt.savefig(f"../../outputs/figures/supplementary/supp_sr_{element}.pdf", bbox_inches="tight")
    plt.savefig(f"../../outputs/figures/supplementary/supp_sr_{element}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved supp_sr_{element}.pdf/.png to outputs/figures/supplementary/")


if __name__ == "__main__":
    for el in ["Mg", "Ca", "Fe", "Ti"]:
        make_figure(el)
