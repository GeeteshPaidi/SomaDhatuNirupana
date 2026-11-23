"""
Reviewer R1.2: compare Mg/Si and Al/Si between two independent epochs processed with the
same v5 pipeline: epoch A (Nov 2024 - May 2026, the mapping set) and epoch B
(Feb 2022 - Dec 2023). The epochs differ in solar activity and detector line width.

Inputs:  outputs/asr/maps/maps_2deg_{tagA}.npz, maps_2deg_{tagB}.npz,
         outputs/asr/maps/mapping_set_{tagA,tagB}.parquet, months_dir/*/calib.json
Outputs: outputs/asr/epochs/epoch_compare.json, epoch_compare.png/.pdf
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
MAPS = os.path.join(ROOT, "outputs", "asr", "maps")
OUT = os.path.join(ROOT, "outputs", "asr", "epochs")
ELS = ["Mg", "Al"]


def coarsen(m, e, f):
    """Inverse-variance combine 2-degree cells into f x f blocks."""
    mean, se, n = m[f"{e}_mean"], m[f"{e}_se"], m[f"{e}_n"]
    w = np.where(np.isfinite(mean) & np.isfinite(se) & (se > 0), 1 / se ** 2, 0.0)
    r0 = np.where(w > 0, mean, 0.0)
    nl, no = mean.shape[0] // f, mean.shape[1] // f

    def blk(a):
        return a[:nl * f, :no * f].reshape(nl, f, no, f).sum((1, 3))
    W, WR, N = blk(w), blk(w * r0), blk(np.nan_to_num(n))
    with np.errstate(invalid="ignore", divide="ignore"):
        return WR / W, 1 / np.sqrt(W), N


def epoch_context(months_dir, months, mapping_set):
    fw, tm = [], []
    for cj in sorted(glob.glob(os.path.join(months_dir, "*", "calib.json"))):
        tag = os.path.basename(os.path.dirname(cj))
        if len(tag) != 7 or tag not in months:
            continue
        c = json.load(open(cj))
        bins = c.get("temp_cal") or []
        if bins:  # TEMP-binned calibration: one FWHM per bin, weighted by spectra per bin
            for b in bins:
                fw += [b.get("fwhm_1.49keV_eV", np.nan)] * max(1, int(b.get("n", 1)) // 1000)
        else:
            fw.append(c.get("fwhm_1.49keV_eV", np.nan))
    ms = pd.read_parquet(mapping_set, columns=["goes_T_mk", "goes_xrsb", "month"])
    return {"n_months": int(ms["month"].nunique()), "n_footprints": int(len(ms)),
            "goes_T_mk_median": float(ms.goes_T_mk.median()),
            "goes_T_mk_p5_p95": [float(ms.goes_T_mk.quantile(.05)), float(ms.goes_T_mk.quantile(.95))],
            "goes_xrsb_median": float(ms.goes_xrsb.median()),
            "fwhm_Al_eV_median": float(np.nanmedian(fw)) if fw else None,
            "fwhm_Al_eV_range": [float(np.nanmin(fw)), float(np.nanmax(fw))] if fw else None}


def compare(A, B, e, f, min_n):
    if f == 1:
        ma, sa, na = A[f"{e}_mean"], A[f"{e}_se"], A[f"{e}_n"]
        mb, sb, nb = B[f"{e}_mean"], B[f"{e}_se"], B[f"{e}_n"]
    else:
        ma, sa, na = coarsen(A, e, f)
        mb, sb, nb = coarsen(B, e, f)
    ok = np.isfinite(ma) & np.isfinite(mb) & (na >= min_n) & (nb >= min_n)
    a, b = ma[ok], mb[ok]
    z = (b - a) / np.sqrt(sa[ok] ** 2 + sb[ok] ** 2)
    d = b - a
    rng = np.random.default_rng(0)
    boot = []
    for _ in range(1000):
        i = rng.integers(0, len(a), len(a))
        boot.append(np.median(b[i] / a[i]))
    return {"block_deg": 2 * f, "min_footprints_each": min_n, "n_cells": int(ok.sum()),
            "pearson_r": float(stats.pearsonr(a, b)[0]), "spearman_rho": float(stats.spearmanr(a, b)[0]),
            "median_A": float(np.median(a)), "median_B": float(np.median(b)),
            "median_ratio_B_over_A": float(np.median(b / a)),
            "median_ratio_ci95": [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
            "median_diff_B_minus_A": float(np.median(d)),
            "diff_p16_p84": [float(np.percentile(d, 16)), float(np.percentile(d, 84))],
            "robust_std_z": float(1.4826 * np.median(np.abs(z - np.median(z)))),
            "_a": a, "_b": b}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tagA", default="_v5_eA")
    ap.add_argument("--tagB", default="_v5_eB")
    ap.add_argument("--monthsA", default="2024_11:2026_05")
    ap.add_argument("--monthsB", default="2022_02:2023_12")
    ap.add_argument("--months_dir", default=os.path.join(ROOT, "outputs", "asr", "months_v5"))
    ap.add_argument("--min_n", type=int, default=20)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    from maps import parse_months
    A = dict(np.load(os.path.join(MAPS, f"maps_2deg{args.tagA}.npz")))
    B = dict(np.load(os.path.join(MAPS, f"maps_2deg{args.tagB}.npz")))
    allm = [os.path.basename(d) for d in glob.glob(os.path.join(args.months_dir, "*"))]
    res = {"epochs": {}}
    for k, tag, spec in [("A", args.tagA, args.monthsA), ("B", args.tagB, args.monthsB)]:
        keep = parse_months(spec)
        res["epochs"][k] = {"months": spec, **epoch_context(
            args.months_dir, [m for m in allm if len(m) == 7 and keep(m)],
            os.path.join(MAPS, f"mapping_set{tag}.parquet"))}
    res["compare"] = {}
    keep_pts = {}
    for e in ELS:
        rows = []
        for f, mn in [(1, args.min_n), (3, args.min_n), (5, args.min_n)]:
            r = compare(A, B, e, f, mn)
            if f == 5:
                keep_pts[e] = (r["_a"], r["_b"])
            rows.append({k: v for k, v in r.items() if not k.startswith("_")})
        res["compare"][e] = rows
    with open(os.path.join(OUT, "epoch_compare.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.3))
    for ax, e in zip(axes, ELS):
        a, b = keep_pts[e]
        lo, hi = np.percentile(np.r_[a, b], [0.5, 99.5])
        ax.scatter(a, b, s=6, alpha=0.5, color="#2a6fbb", lw=0)
        ax.plot([lo, hi], [lo, hi], color="0.3", lw=1)
        c = res["compare"][e][2]
        ax.set_title(f"{e}/Si, 10° blocks (n = {c['n_cells']}), r = {c['pearson_r']:.2f}", fontsize=10)
        ax.set_xlabel(f"{e}/Si, Nov 2024 – May 2026")
        ax.set_ylabel(f"{e}/Si, Feb 2022 – Dec 2023")
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect("equal")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"epoch_compare.{ext}"), dpi=150)


if __name__ == "__main__":
    import sys
    sys.path.insert(0, os.path.dirname(__file__))
    main()
