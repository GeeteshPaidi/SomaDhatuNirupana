"""
POST-HOC SENSITIVITY (requested after validation; clearly labelled as such):
inverse-variance-weighted MEDIAN instead of weighted mean, per 2-deg cell
(footprint centres) and per GRS 5-deg pixel. Reports the GRS correlation and the
split-half reliability for both statistics, and writes the alternative map figure.
The pre-declared product remains the weighted mean (maps.py).
"""
import json
import os
import sys

import matplotlib
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(__file__))
import validate_asr as va  # noqa: E402

ROOT = va.ROOT
ELS = ["Mg", "Al"]
M = {"Mg": 24.305 / 40.304, "Al": 2 * 26.982 / 101.961, "Si": 28.085 / 60.084}


def wmedian(x, w):
    ok = np.isfinite(x) & np.isfinite(w) & (w > 0)
    if ok.sum() == 0:
        return np.nan
    x, w = x[ok], w[ok]
    o = np.argsort(x)
    c = np.cumsum(w[o])
    return x[o][np.searchsorted(c, 0.5 * c[-1])]


def group_stat(keys, x, w, stat):
    df = pd.DataFrame({"k": keys, "x": x, "w": w}).dropna()
    if stat == "median":
        return df.groupby("k").apply(lambda g: wmedian(g.x.values, g.w.values), include_groups=False)
    return df.groupby("k").apply(lambda g: np.sum(g.w * g.x) / np.sum(g.w), include_groups=False)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--mapping_set", default=os.path.join(ROOT, "outputs", "asr", "maps", "mapping_set.parquet"))
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    ms = pd.read_parquet(args.mapping_set)
    g = va.read_grs(os.path.join(ROOT, "data", "grs", "lpgrs_high1_elem_abundance_5deg.tab"))
    lat, lon = ms.clat.values, ms.clon.values
    pix = np.full(len(ms), -1)
    for i, p in g.iterrows():
        inside = (lat >= p.MIN_LAT) & (lat < p.MAX_LAT) & (lon >= p.MIN_LON) & (lon < p.MAX_LON)
        pix[inside] = i
    months = sorted(ms.month.unique())
    half = ms.month.map({m: i % 2 for i, m in enumerate(months)}).values
    out = {}
    for stat in ["mean", "median"]:
        for e in ELS:
            x = ms[f"r_{e}"].values
            w = 1 / np.maximum(ms[f"sd_r_{e}"].values, 1e-6) ** 2
            v = group_stat(pix, x, w, stat).drop(-1, errors="ignore")
            grs = g.loc[v.index, "W_MGO" if e == "Mg" else "W_AL2O3"] * M[e] / (g.loc[v.index, "W_SIO2"] * M["Si"])
            out[f"{stat}_{e}_grs_r"] = float(pearsonr(grs.values, v.values)[0])
            h = [group_stat(pix[half == k], x[half == k], w[half == k], stat) for k in (0, 1)]
            j = pd.concat(h, axis=1).dropna()
            j = j.drop(-1, errors="ignore")
            rh = pearsonr(j.iloc[:, 0], j.iloc[:, 1])[0]
            out[f"{stat}_{e}_split_half_r"] = float(rh)
            out[f"{stat}_{e}_reliability"] = float(2 * rh / (1 + rh))
    # bootstrap CI (20-deg blocks) for the median version
    rng = np.random.default_rng(20261010)
    for e in ELS:
        x = ms[f"r_{e}"].values
        w = 1 / np.maximum(ms[f"sd_r_{e}"].values, 1e-6) ** 2
        v = group_stat(pix, x, w, "median").drop(-1, errors="ignore")
        gg = g.loc[v.index]
        grs = gg["W_MGO" if e == "Mg" else "W_AL2O3"] * M[e] / (gg["W_SIO2"] * M["Si"])
        clat = (gg.MIN_LAT + gg.MAX_LAT) / 2
        clon = (gg.MIN_LON + gg.MAX_LON) / 2
        blk = (np.floor((clat + 90) / 20) * 100 + np.floor((clon + 180) / 20)).values
        s = va.corr_with_ci(grs.values, v.values, rng, n_boot=2000, groups=blk)
        out[f"median_{e}_grs_ci95"] = s["pearson_ci95"]
    # 2-deg weighted-median map (footprint centres)
    res = 2.0
    ri = np.clip(((lat + 90) / res).astype(int), 0, 89)
    ci = np.clip(((lon + 180) / res).astype(int), 0, 179)
    cell = ri * 180 + ci
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.4))
    for ax, (e, lo, hi) in zip(axes, [("Mg", 0.1, 0.4), ("Al", 0.2, 1.3)]):
        x = ms[f"r_{e}"].values
        w = 1 / np.maximum(ms[f"sd_r_{e}"].values, 1e-6) ** 2
        v = group_stat(cell, x, w, "median")
        a = np.full(90 * 180, np.nan)
        a[v.index.values.astype(int)] = v.values
        im = ax.imshow(a.reshape(90, 180), origin="lower", extent=[-180, 180, -90, 90], cmap="viridis",
                       vmin=lo, vmax=hi, aspect="auto")
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
        cb.set_label(f"{e}/Si (weighted median)", fontsize=7)
        ax.set_title(f"{e}/Si — POST-HOC weighted-median alternative", fontsize=7)
        ax.tick_params(labelsize=6)
    fig.tight_layout()
    fd = os.path.join(ROOT, "paper_asr", "figures")
    fig.savefig(os.path.join(fd, f"fig_maps_wmedian_posthoc{args.tag}.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(fd, f"fig_maps_wmedian_posthoc{args.tag}.png"), dpi=200, bbox_inches="tight")
    with open(os.path.join(va.OUT, f"alt_median_posthoc{args.tag}.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
