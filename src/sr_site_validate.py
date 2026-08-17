"""
Held-out check of the super-resolved product itself (reviewer question:
"how do you know A_SR is better, not just different, from the pre-SR
input?"), for the Al/Si headline result (Fig. 3 / Sec. VI-VII).

No retraining is involved -- the SR fusion (src/superres.py) never uses
Apollo/Luna data, so there is no leakage to guard against by excluding a
site. This script simply samples both the pre-SR (abundance_interp) and
post-SR (A_sr) Al/Si maps at the 9 Apollo/Luna site locations already
used in src/validate.py's landing-site validation, and compares each to
the same reference geochemistry (data/reference/validation.xlsx), using
the same +-5 deg matching window as validate.py's footprint radius.

Pearson r is the primary comparison (scale/offset invariant, so it is
unaffected by any difference between this map's native units and the
raw per-footprint ratio units used in Table V). RMSE/bias are also
reported for pre- vs. post-SR since both draw from the same map array
(A_sr is built from abundance_interp, so they share the same native
scale) -- but they are not expected to numerically match Table V, which
is computed independently from raw per-footprint fits, not the
map_making-gridded array used here.
"""
import argparse

import numpy as np
from scipy.stats import pearsonr

from validate import SITE_COORDS, load_reference


def sample_site(arr, lat0, lon0, radius_deg, grid_extent=((-90, 90), (-180, 180))):
    h, w = arr.shape
    (lat_min, lat_max), (lon_min, lon_max) = grid_extent
    lat_per_row = (lat_max - lat_min) / h
    lon_per_col = (lon_max - lon_min) / w
    row0 = (lat0 - lat_min) / lat_per_row
    col0 = (lon0 - lon_min) / lon_per_col
    dr = int(np.ceil(radius_deg / lat_per_row))
    dc = int(np.ceil(radius_deg / lon_per_col))
    r_lo, r_hi = max(0, int(row0 - dr)), min(h, int(row0 + dr) + 1)
    c_lo, c_hi = max(0, int(col0 - dc)), min(w, int(col0 + dc) + 1)
    window = arr[r_lo:r_hi, c_lo:c_hi]
    vals = window[~np.isnan(window)]
    if vals.size == 0:
        return None
    return float(np.median(vals))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sr_result", default="../outputs/figures/sr_ipynb_lat60/sr_result_Al.npz")
    ap.add_argument("--reference", default="../data/reference/validation.xlsx")
    ap.add_argument("--radius_deg", type=float, default=5.0)
    args = ap.parse_args()

    d = np.load(args.sr_result)
    pre_sr = d["abundance_interp"]
    post_sr = d["A_sr"]
    alpha = float(d["alpha"])
    print(f"Loaded SR result: alpha={alpha}, grid shape {pre_sr.shape}")

    ref = load_reference(args.reference)
    print(f"Loaded reference geochemistry for {len(ref)} sites")

    rows = []
    for site, info in ref.items():
        lat0, lon0 = info["lat"], info["lon"]
        ref_ratio = info["ratio_to_Si"]["Al"]
        pre_val = sample_site(pre_sr, lat0, lon0, args.radius_deg)
        post_val = sample_site(post_sr, lat0, lon0, args.radius_deg)
        if pre_val is None or post_val is None:
            print(f"  site {site}: no valid pixels within {args.radius_deg} deg -- skipped")
            continue
        rows.append({"site": site, "reference": ref_ratio, "pre_sr": pre_val, "post_sr": post_val})
        print(f"  site {site}: reference={ref_ratio:.4f}  pre-SR={pre_val:.4f}  post-SR={post_val:.4f}")

    ref_v = np.array([r["reference"] for r in rows])
    pre_v = np.array([r["pre_sr"] for r in rows])
    post_v = np.array([r["post_sr"] for r in rows])
    n = len(rows)

    def stats(fit_v):
        r, _ = pearsonr(ref_v, fit_v)
        rmse = np.sqrt(np.mean((ref_v - fit_v) ** 2))
        bias = np.mean(fit_v - ref_v)
        return r, rmse, bias

    r_pre, rmse_pre, bias_pre = stats(pre_v)
    r_post, rmse_post, bias_post = stats(post_v)

    print(f"\nAl/Si vs. Apollo/Luna reference, n={n} sites (+-{args.radius_deg} deg window, median pixel value):")
    print(f"  pre-SR  (input) : r={r_pre:+.3f}  RMSE={rmse_pre:.3f}  bias={bias_pre:+.3f}")
    print(f"  post-SR (A_sr)  : r={r_post:+.3f}  RMSE={rmse_post:.3f}  bias={bias_post:+.3f}")
    print(f"  delta r (post-pre): {r_post - r_pre:+.3f}")


if __name__ == "__main__":
    main()
