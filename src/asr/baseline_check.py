"""
Task 0.2: confirm the OLD pipeline (src/fitting) still runs and reproduces
the stored results, and record the old element/Si ratios as a baseline.

Steps
  1. Re-preprocess N raw FITS from one month with prepare_data_parallel
     (month nightside background, as run_pipeline.prepare_month does) and
     compare with the stored month catalogue (data/processed/months).
  2. Re-fit the same spectra with run_fit.fit_single_row and compare the
     concentrations with the stored expanded result.
  3. Global baseline: percentiles of X/Si over all successful fits in the
     stored expanded result (855,190 rows), and per-site validation with
     the old validate.py logic (5 deg radius) to check r(Al/Si)=0.738.

Writes outputs/asr/baseline/*.csv and baseline_summary.json.
"""
import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "src", "fitting"))
sys.path.insert(0, os.path.join(ROOT, "src"))
import prepare_data_parallel as prep  # noqa: E402
from run_fit import fit_single_row  # noqa: E402
import validate as oldval  # noqa: E402

prep.tqdm = lambda x, **k: x
ELS = ["Fe", "Al", "Mg", "Ca", "Ti", "O"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--month", default="2025/10")
    ap.add_argument("--n", type=int, default=300)
    ap.add_argument("--out_dir", default=os.path.join(ROOT, "outputs", "asr", "baseline"))
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    summary = {}

    # ---- 1. preprocessing reproduction ----
    mdir = os.path.join(ROOT, "data", "raw", *args.month.split("/"))
    tag = args.month.replace("/", "_")
    paths = sorted(glob.glob(os.path.join(mdir, "**", "*.fits"), recursive=True))
    bg, n_bg = prep.build_background(paths, max_bg_files=3000)
    prep._init_worker(bg)
    cat = pd.read_parquet(os.path.join(ROOT, "data", "processed", "months", f"cat_{tag}.parquet"))
    cat_idx = {s: i for i, s in enumerate(cat["source_file"])}
    rng = np.random.default_rng(0)
    # sample among dayside files present in the stored catalogue
    day = [p for p in paths if os.path.basename(p) in cat_idx]
    pick = [day[i] for i in sorted(rng.choice(len(day), size=min(args.n, len(day)), replace=False))]
    diffs, rows = [], []
    for p in pick:
        r = prep._process_one(p)
        stored = np.asarray(cat["counts"].iloc[cat_idx[os.path.basename(p)]], dtype=np.float64)
        diffs.append(np.max(np.abs(r["counts"].astype(np.float64) - stored)) / max(stored.max(), 1e-12))
        rows.append(r)
    summary["month"] = args.month
    summary["n_files_month"] = len(paths)
    summary["n_dayside_in_catalogue"] = len(cat)
    summary["n_nightside_bg"] = int(n_bg)
    summary["prep_max_rel_diff"] = float(np.max(diffs))

    # ---- 2. fit reproduction ----
    fits_new = [fit_single_row((r["counts"], r["footprint"])) for r in rows]
    new = pd.DataFrame(fits_new)
    new["source_file"] = [r["source_file"] for r in rows]
    old = pd.read_parquet(os.path.join(ROOT, "data", "processed", "catalogue_expanded_result.parquet"),
                          columns=["source_file", "success", "scale"] + [f"conc_{e}" for e in ELS + ["Si"]])
    m = new.merge(old, on="source_file", suffixes=("_new", "_old"))
    summary["n_refit"] = len(new)
    summary["n_refit_success"] = int(new["success"].sum())
    summary["n_refit_in_stored_result"] = len(m)
    for e in ELS + ["Si"]:
        a, b = m[f"conc_{e}_new"].values, m[f"conc_{e}_old"].values
        summary[f"refit_max_abs_diff_conc_{e}"] = float(np.nanmax(np.abs(a - b))) if len(m) else None
    new.to_csv(os.path.join(args.out_dir, f"refit_{tag}.csv"), index=False)
    ratios_new = pd.DataFrame({f"{e}/Si": new[f"conc_{e}"] / new["conc_Si"] for e in ELS})
    summary["refit_ratio_median"] = ratios_new.median().to_dict()

    # ---- 3. global baseline from stored expanded result ----
    old = old[old["success"] == True]  # noqa: E712
    rat = pd.DataFrame({f"{e}/Si": old[f"conc_{e}"] / old["conc_Si"] for e in ELS}).replace([np.inf, -np.inf], np.nan)
    q = rat.quantile([0.05, 0.25, 0.5, 0.75, 0.95]).T
    q["frac_conc_Si_lt_1e-3"] = float((old["conc_Si"] < 1e-3).mean())
    q.to_csv(os.path.join(args.out_dir, "baseline_ratio_percentiles.csv"))
    summary["n_stored_success"] = len(old)
    summary["baseline_ratio_percentiles"] = q.round(5).to_dict(orient="index")
    for e in ELS + ["Si"]:
        summary[f"frac_at_upper_bound_{e}"] = float((old[f"conc_{e}"] >= 199.9).mean())
        summary[f"frac_at_lower_bound_{e}"] = float((old[f"conc_{e}"] <= 1e-6).mean())

    # ---- 4. old validation on the expanded result (5 deg, median) ----
    fp = pd.read_parquet(os.path.join(ROOT, "data", "processed", "catalogue_expanded_result.parquet"),
                         columns=["footprint", "success"] + [f"conc_{e}" for e in ELS + ["Si"]])
    fp = fp[fp["success"] == True]  # noqa: E712
    ref = oldval.load_reference(os.path.join(ROOT, "data", "reference", "validation.xlsx"))
    vrows = []
    for site, info in ref.items():
        mt = oldval.match_footprints(fp, info["lat"], info["lon"], 5.0)
        for e in ["Al", "Fe", "Mg", "Ca", "Ti"]:
            fr = (mt[f"conc_{e}"] / mt["conc_Si"]).replace([np.inf, -np.inf], np.nan).dropna()
            vrows.append({"site": site, "element": e, "n": len(mt), "fit": fr.median(),
                          "ref": info["ratio_to_Si"][e]})
    vt = pd.DataFrame(vrows)
    vt.to_csv(os.path.join(args.out_dir, "baseline_validation_expanded_5deg.csv"), index=False)
    summary["baseline_validation_r"] = {e: float(pearsonr(g["ref"], g["fit"])[0])
                                        for e, g in vt.groupby("element")}

    with open(os.path.join(args.out_dir, "baseline_summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=float)
    print(json.dumps(summary, indent=2, default=float))


if __name__ == "__main__":
    main()
