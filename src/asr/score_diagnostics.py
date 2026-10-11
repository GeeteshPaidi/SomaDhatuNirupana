"""
Why are the validation correlations what they are? (post-validation diagnostics,
no analysis choice is changed by this script)

1. Site-test ceiling: Lunar Prospector GRS (the reference orbital product) against
   the same nine returned-soil averages, using the GRS 5-deg pixel containing each site.
2. Internal reliability of the CLASS maps: split the mapping set into two halves by
   month (alternating), compute per-GRS-pixel means for each half, correlate them,
   and Spearman-Brown-correct to full-sample reliability.
3. GRS reliability from its own error matrix: 1 - mean(var_err)/var(observed ratio).
4. Attenuation ceiling for the CLASS-GRS correlation: sqrt(rel_CLASS * rel_GRS).
5. Earlier processing of the same data against GRS on the same pixels (median ratio of the
   855,190-row expanded result), as a baseline.
Output: outputs/asr/validation/score_diagnostics.json
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

sys.path.insert(0, os.path.dirname(__file__))
import validate_asr as va  # noqa: E402
import reference  # noqa: E402

ROOT = va.ROOT
OUT = va.OUT
ELS = ["Mg", "Al"]
M = {"Mg": 24.305 / 40.304, "Al": 2 * 26.982 / 101.961, "Si": 28.085 / 60.084}
IDX = {"Mg": 0, "Al": 1, "Si": 2}  # error-matrix indices (MgO, Al2O3, SiO2)


def err_index(i, j):
    """Column name for E[i,j] (upper triangle, 9x9) in read_grs naming (E0..E44)."""
    i, j = min(i, j), max(i, j)
    k = sum(9 - a for a in range(i)) + (j - i)
    return f"E{k}"


def grs_ratio_var(g, e):
    """Variance of (w_e*f_e)/(w_Si*f_Si) by first-order propagation of the error matrix."""
    a, b = g[f"W_{'MGO' if e == 'Mg' else 'AL2O3'}"], g["W_SIO2"]
    va_, vb = g[err_index(IDX[e], IDX[e])], g[err_index(2, 2)]
    cab = g[err_index(IDX[e], 2)]
    r = a / b
    var = r ** 2 * (va_ / a ** 2 + vb / b ** 2 - 2 * cab / (a * b))
    return (M[e] / M["Si"]) ** 2 * var


def pixel_of(g, lat, lon):
    m = (lat >= g.MIN_LAT) & (lat < g.MAX_LAT) & (lon >= g.MIN_LON) & (lon < g.MAX_LON)
    return g[m].iloc[0] if m.any() else None


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--mapping_set", default=os.path.join(ROOT, "outputs", "asr", "maps", "mapping_set.parquet"))
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    out = {}
    g = va.read_grs(os.path.join(ROOT, "data", "grs", "lpgrs_high1_elem_abundance_5deg.tab"))
    ref = reference.load_reference(os.path.join(ROOT, "data", "reference", "validation.xlsx"))
    # 1. GRS vs sites
    rows = []
    for s, info in ref.items():
        p = pixel_of(g, info["lat"], info["lon"])
        rows.append({"site": s, **{f"grs_{e}": p[f"W_{'MGO' if e=='Mg' else 'AL2O3'}"] * M[e] / (p.W_SIO2 * M["Si"]) for e in ELS},
                     **{f"ref_{e}": info["ratio_to_Si"][e] for e in ELS}})
    t = pd.DataFrame(rows)
    out["grs_vs_sites_r"] = {e: float(pearsonr(t[f"ref_{e}"], t[f"grs_{e}"])[0]) for e in ELS}
    t.to_csv(os.path.join(OUT, "grs_vs_sites.csv"), index=False)

    # 2. CLASS split-half reliability on GRS pixels
    ms = pd.read_parquet(args.mapping_set)
    months = sorted(ms.month.unique())
    half = ms.month.map({m: i % 2 for i, m in enumerate(months)}).values
    lat, lon = ms.clat.values, ms.clon.values
    per = {e: [[], []] for e in ELS}
    gvar = {e: [] for e in ELS}
    gval = {e: [] for e in ELS}
    for _, p in g.iterrows():
        inside = (lat >= p.MIN_LAT) & (lat < p.MAX_LAT) & (lon >= p.MIN_LON) & (lon < p.MAX_LON)
        for e in ELS:
            ok = True
            vals = []
            for h in (0, 1):
                sel = inside & (half == h)
                m, se, n = va.wmean(ms[f"r_{e}"].values[sel], ms[f"sd_r_{e}"].values[sel])
                if n < 10:
                    ok = False
                vals.append(m)
            if ok:
                per[e][0].append(vals[0])
                per[e][1].append(vals[1])
            gval[e].append(p[f"W_{'MGO' if e=='Mg' else 'AL2O3'}"] * M[e] / (p.W_SIO2 * M["Si"]))
            gvar[e].append(grs_ratio_var(p, e))
    rel_c, rel_g = {}, {}
    for e in ELS:
        r_half = pearsonr(per[e][0], per[e][1])[0]
        rel_c[e] = 2 * r_half / (1 + r_half)
        gv = np.array(gvar[e])
        gv = gv[np.isfinite(gv) & (gv > 0)]
        rel_g[e] = 1 - np.median(gv) / np.var(gval[e])
        out[f"class_split_half_r_{e}"] = float(r_half)
        out[f"class_reliability_{e}"] = float(rel_c[e])
        out[f"grs_reliability_{e}"] = float(rel_g[e])
        out[f"class_grs_ceiling_r_{e}"] = float(np.sqrt(max(rel_c[e], 0) * max(rel_g[e], 0)))
        out[f"n_pixels_split_{e}"] = len(per[e][0])

    # 5. old pipeline vs GRS
    old = pd.read_parquet(os.path.join(ROOT, "data", "processed", "catalogue_expanded_result.parquet"),
                          columns=["footprint", "success", "conc_Mg", "conc_Al", "conc_Si"])
    old = old[old.success == True]  # noqa: E712
    fp = np.array([np.mean(np.array(list(x)), axis=0) for x in old.footprint.values])
    olat, olon = fp[:, 0], ((fp[:, 1] + 180) % 360) - 180
    t2 = va.grs_compare.__wrapped__ if hasattr(va.grs_compare, "__wrapped__") else None  # noqa: F841
    res = {e: ([], []) for e in ELS}
    for _, p in g.iterrows():
        inside = (olat >= p.MIN_LAT) & (olat < p.MAX_LAT) & (olon >= p.MIN_LON) & (olon < p.MAX_LON)
        if inside.sum() < 20:
            continue
        for e in ELS:
            rr = (old[f"conc_{e}"].values[inside] / old["conc_Si"].values[inside])
            rr = rr[np.isfinite(rr)]
            res[e][0].append(p[f"W_{'MGO' if e=='Mg' else 'AL2O3'}"] * M[e] / (p.W_SIO2 * M["Si"]))
            res[e][1].append(np.median(rr))
    for e in ELS:
        out[f"old_pipeline_vs_grs_r_{e}"] = float(pearsonr(res[e][0], res[e][1])[0])
        out[f"old_pipeline_n_pixels_{e}"] = len(res[e][0])
    with open(os.path.join(OUT, f"score_diagnostics{args.tag}.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
