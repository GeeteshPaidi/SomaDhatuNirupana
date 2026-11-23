"""
Revision (R2-2): CLASS Mg/Si and Al/Si at post-Apollo in-situ and sample-return sites,
with the same 1/2/5 deg great-circle weighted means as validate_asr.py.
Reference oxide compositions are entered by hand from the cited papers.
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

from validate_asr import gc_dist_deg, wmean, RADII

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
F = {"Mg": 24.305 / 40.304, "Al": 2 * 26.982 / 101.961, "Si": 28.0855 / 60.084}  # element / oxide mass

# lat, lon (deg E), oxide wt% {SiO2, Al2O3, MgO} or None, source
SITES = {
    "Chang'E-3": (44.12, -19.51, {"SiO2": 40.9, "Al2O3": 11.1, "MgO": 10.7},
                  "Ng et al. 2023, Remote Sens. 15, 1643, Table 9 (APXS, soil LS1 mean)"),
    "Chang'E-5": (43.058, -51.916, {"SiO2": 42.2, "Al2O3": 10.8, "MgO": 6.48},
                  "Li et al. 2022, NSR 9, nwab188, Table 1 (soil CE5C0800YJFM002)"),
    "Chang'E-6": (-41.625, -153.978, {"SiO2": 45.6, "Al2O3": 14.3, "MgO": 7.08},
                  "Li et al. 2024, NSR 11, nwae328, Table 2 (soil CE6C0000YJFM00102)"),
    "Chandrayaan-3": (-69.3676, 32.3481, {"SiO2": 43.3, "Al2O3": 26.1, "MgO": 8.1},
                      "Vadawale et al. 2024, Nature 633, 327, Table 1 (APXS, mean of 23 observations)"),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mapping_set", default=os.path.join(ROOT, "outputs", "asr", "maps", "mapping_set_v5.parquet"))
    ap.add_argument("--out", default=os.path.join(ROOT, "revision", "new_sites.csv"))
    a = ap.parse_args()
    ms = pd.read_parquet(a.mapping_set)
    rows = []
    for site, (lat, lon, ox, src) in SITES.items():
        d = gc_dist_deg(ms["clat"].values, ms["clon"].values, lat, lon)
        ref = None
        if ox:
            si = ox["SiO2"] * F["Si"]
            ref = {"Mg": ox["MgO"] * F["Mg"] / si, "Al": ox["Al2O3"] * F["Al"] / si}
        for rad in RADII:
            sel = ms[d <= rad]
            for e in ("Mg", "Al"):
                m, se, n = wmean(sel[f"r_{e}"].values, sel[f"sd_r_{e}"].values)
                rows.append({"site": site, "lat": lat, "lon": lon, "radius_deg": rad, "element": e,
                             "n_footprints": n, "n_months": int(sel["month"].nunique()),
                             "class_ratio": m, "class_se": se,
                             "class_median": float(np.median(sel[f"r_{e}"])) if n else np.nan,
                             "ref_ratio": ref[e] if ref else np.nan, "ref_source": src})
    t = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    t.to_csv(a.out, index=False)
    print(t.drop(columns="ref_source").to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
