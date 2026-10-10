"""
Phase 4 validation (TRACKER 4.1, 4.2), with analysis choices fixed in the
decisions log before any numbers were computed.

4.1 Apollo/Luna (9 sites, average soil compositions, data/reference/validation.xlsx):
    radii 1, 2, 5 deg (primary 2 deg, great-circle on footprint centres),
    per-site inverse-variance-weighted mean ratio + propagated standard error
    (and site scatter), Pearson/Spearman r, bootstrap 95% CI over sites
    (10,000 resamples), leave-one-site-out r range.
4.2 Lunar Prospector GRS (Prettyman et al. 2006; PDS LP-L-GRS-5-ELEM-ABUNDANCE-V1.0),
    native 5-deg equal-area pixels: CLASS weighted mean of footprint centres in
    each pixel (>= MIN_FP footprints) vs GRS Mg/Si, Al/Si weight ratios from
    oxides. Pearson + Spearman with spatial block bootstrap (20x20 deg blocks).
Claim rule: element claimed only if Pearson r > 0 with 95% CI excluding 0 in
BOTH the primary-radius site test and the GRS test.

Outputs: outputs/asr/validation/*.csv, *.json, figures.
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import validate as oldval  # noqa: E402  (reference loader + site coordinates)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.join(ROOT, "outputs", "asr", "validation")
ELS = ["Mg", "Al"]
RADII = [1.0, 2.0, 5.0]
PRIMARY_RADIUS = 2.0
N_BOOT = 10000
MIN_FP = 20
BLOCK_DEG = 20.0
OX = {"Mg": ("W_MGO", 24.305 / 40.304), "Al": ("W_AL2O3", 2 * 26.982 / 101.961),
      "Si": ("W_SIO2", 28.085 / 60.084), "Ca": ("W_CAO", 40.078 / 56.077)}


def gc_dist_deg(lat1, lon1, lat2, lon2):
    a1, o1, a2, o2 = map(np.radians, (lat1, lon1, lat2, lon2))
    return np.degrees(2 * np.arcsin(np.sqrt(np.sin((a2 - a1) / 2) ** 2 +
                                            np.cos(a1) * np.cos(a2) * np.sin((o2 - o1) / 2) ** 2)))


def wmean(r, sd):
    ok = np.isfinite(r) & np.isfinite(sd) & (sd > 0)
    if ok.sum() == 0:
        return np.nan, np.nan, 0
    w = 1 / sd[ok] ** 2
    m = np.sum(w * r[ok]) / w.sum()
    return m, 1 / np.sqrt(w.sum()), int(ok.sum())


def corr_with_ci(x, y, rng, n_boot=N_BOOT, groups=None):
    """Pearson/Spearman with bootstrap CI; resample points, or blocks if groups given."""
    x, y = np.asarray(x), np.asarray(y)
    pr, pp = pearsonr(x, y)
    sr, sp = spearmanr(x, y)
    bp, bs = [], []
    if groups is None:
        idx_sets = [rng.integers(0, len(x), len(x)) for _ in range(n_boot)]
    else:
        ug = np.unique(groups)
        members = {g: np.where(groups == g)[0] for g in ug}
        idx_sets = [np.concatenate([members[g] for g in rng.choice(ug, len(ug))]) for _ in range(n_boot)]
    for idx in idx_sets:
        xs, ys = x[idx], y[idx]
        if np.std(xs) == 0 or np.std(ys) == 0:
            continue
        bp.append(np.corrcoef(xs, ys)[0, 1])
        bs.append(spearmanr(xs, ys)[0])
    return {"n": int(len(x)), "pearson_r": float(pr), "pearson_p": float(pp),
            "pearson_ci95": [float(np.percentile(bp, 2.5)), float(np.percentile(bp, 97.5))],
            "spearman_rho": float(sr), "spearman_p": float(sp),
            "spearman_ci95": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]}


def sites(ms, rng):
    ref = oldval.load_reference(os.path.join(ROOT, "data", "reference", "validation.xlsx"))
    rows = []
    for site, info in ref.items():
        d = gc_dist_deg(ms["clat"].values, ms["clon"].values, info["lat"], info["lon"])
        for rad in RADII:
            sel = ms[d <= rad]
            for e in ELS:
                m, se, n = wmean(sel[f"r_{e}"].values, sel[f"sd_r_{e}"].values)
                rows.append({"site": site, "lat": info["lat"], "lon": info["lon"], "radius_deg": rad,
                             "element": e, "n_footprints": n, "class_ratio": m, "class_se": se,
                             "class_median": float(np.median(sel[f"r_{e}"])) if n else np.nan,
                             "ref_ratio": info["ratio_to_Si"][e]})
    t = pd.DataFrame(rows)
    stats = []
    for (rad, e), g in t.groupby(["radius_deg", "element"]):
        g = g.dropna(subset=["class_ratio"])
        if len(g) < 3:
            continue
        s = corr_with_ci(g["ref_ratio"].values, g["class_ratio"].values, rng)
        loo = [pearsonr(g["ref_ratio"].values[np.arange(len(g)) != i],
                        g["class_ratio"].values[np.arange(len(g)) != i])[0] for i in range(len(g))]
        s.update(radius_deg=rad, element=e, loo_min=float(np.min(loo)), loo_max=float(np.max(loo)),
                 loo_dropped_site_at_min=g["site"].values[int(np.argmin(loo))],
                 rmse=float(np.sqrt(np.mean((g["class_ratio"] - g["ref_ratio"]) ** 2))),
                 bias=float(np.mean(g["class_ratio"] - g["ref_ratio"])),
                 median_site_se=float(np.median(g["class_se"])))
        stats.append(s)
    return t, pd.DataFrame(stats)


def read_grs(path):
    cols = ["PIXEL_INDEX", "MIN_LAT", "MAX_LAT", "MIN_LON", "MAX_LON", "AM", "NEUTRON_DEN",
            "W_MGO", "W_AL2O3", "W_SIO2", "W_CAO", "W_TIO2", "W_FEO", "W_K", "W_TH", "W_U"]
    g = pd.read_csv(path, sep=r"\s+", header=None)
    g = g.iloc[:, :61]
    names = cols + [f"E{i}" for i in range(61 - len(cols))]
    g.columns = names
    return g


def grs_compare(ms, rng):
    g = read_grs(os.path.join(ROOT, "data", "grs", "lpgrs_high1_elem_abundance_5deg.tab"))
    lat, lon = ms["clat"].values, ms["clon"].values
    rows = []
    for _, p in g.iterrows():
        lo1, lo2 = p.MIN_LON, p.MAX_LON
        inside = (lat >= p.MIN_LAT) & (lat < p.MAX_LAT) & (lon >= lo1) & (lon < lo2)
        row = {"pixel": int(p.PIXEL_INDEX), "clat": (p.MIN_LAT + p.MAX_LAT) / 2, "clon": (lo1 + lo2) / 2,
               "n_footprints": int(inside.sum())}
        si = p[OX["Si"][0]] * OX["Si"][1]
        for e in ELS:
            row[f"grs_{e}_Si"] = p[OX[e][0]] * OX[e][1] / si
            m, se, n = wmean(ms[f"r_{e}"].values[inside], ms[f"sd_r_{e}"].values[inside])
            row[f"class_{e}_Si"], row[f"class_{e}_se"] = m, se
        rows.append(row)
    t = pd.DataFrame(rows)
    t = t[t["n_footprints"] >= MIN_FP].reset_index(drop=True)
    blk = (np.floor((t["clat"] + 90) / BLOCK_DEG) * 100 + np.floor((t["clon"] + 180) / BLOCK_DEG)).values
    stats = []
    for e in ELS:
        ok = np.isfinite(t[f"class_{e}_Si"])
        if ok.sum() < 3:
            stats.append({"element": e, "n": int(ok.sum()), "pearson_r": np.nan, "pearson_ci95": [np.nan, np.nan]})
            continue
        s = corr_with_ci(t.loc[ok, f"grs_{e}_Si"].values, t.loc[ok, f"class_{e}_Si"].values, rng,
                         n_boot=2000, groups=blk[ok.values])
        s.update(element=e, min_footprints=MIN_FP, block_deg=BLOCK_DEG)
        stats.append(s)
    return t, pd.DataFrame(stats)


def plots(site_t, grs_t, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(9, 8))
    for j, e in enumerate(ELS):
        g = site_t[(site_t.radius_deg == PRIMARY_RADIUS) & (site_t.element == e)]
        ax = axes[0, j]
        ax.errorbar(g.ref_ratio, g.class_ratio, yerr=g.class_se, fmt="o", ms=4)
        for _, r in g.iterrows():
            ax.annotate(r.site, (r.ref_ratio, r.class_ratio), fontsize=7)
        ax.set_xlabel(f"Returned-soil {e}/Si")
        ax.set_ylabel(f"CLASS {e}/Si ({PRIMARY_RADIUS:.0f}° radius)")
        ax = axes[1, j]
        ax.scatter(grs_t[f"grs_{e}_Si"], grs_t[f"class_{e}_Si"], s=5, alpha=0.5)
        ax.set_xlabel(f"Lunar Prospector GRS {e}/Si (5° pixel)")
        ax.set_ylabel(f"CLASS {e}/Si (5° pixel)")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "validation_scatter.png"), dpi=120)
    fig.savefig(os.path.join(out, "validation_scatter.pdf"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mapping_set", default=os.path.join(ROOT, "outputs", "asr", "maps", "mapping_set.parquet"))
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(20261010)
    ms = pd.read_parquet(args.mapping_set)
    site_t, site_s = sites(ms, rng)
    grs_t, grs_s = grs_compare(ms, rng)
    site_t.to_csv(os.path.join(OUT, f"sites_table{args.tag}.csv"), index=False)
    site_s.to_csv(os.path.join(OUT, f"sites_stats{args.tag}.csv"), index=False)
    grs_t.to_csv(os.path.join(OUT, f"grs_table{args.tag}.csv"), index=False)
    grs_s.to_csv(os.path.join(OUT, f"grs_stats{args.tag}.csv"), index=False)
    claims = {}
    for e in ELS:
        s1 = site_s[(site_s.radius_deg == PRIMARY_RADIUS) & (site_s.element == e)]
        s2 = grs_s[grs_s.element == e]
        ok1 = len(s1) and s1.pearson_r.iloc[0] > 0 and s1.pearson_ci95.iloc[0][0] > 0
        ok2 = len(s2) and np.isfinite(s2.pearson_r.iloc[0]) and s2.pearson_r.iloc[0] > 0             and s2.pearson_ci95.iloc[0][0] > 0
        claims[e] = {"sites_pass": bool(ok1), "grs_pass": bool(ok2), "claimed": bool(ok1 and ok2)}
    with open(os.path.join(OUT, f"claims{args.tag}.json"), "w") as f:
        json.dump(claims, f, indent=1)
    plots(site_t, grs_t, OUT)
    pd.set_option("display.width", 250)
    print(site_s.to_string())
    print(grs_s.to_string())
    print(json.dumps(claims, indent=1))


if __name__ == "__main__":
    main()
