"""
Reviewer 2, Major 5: observed 8 s spectra, fitted v5 models and deviance residuals at
the two locations discussed by Narendranath et al. (2024).

For each location, the quality-0 mapping-set spectra whose footprint centre lies within
--radius degrees are ranked by counts; the --n brightest are refitted with the v5
configuration (TEMP-binned detector calibration, ARF, gamma = 0, 8-17 keV proxy band)
and plotted with diag_plots.plot_fit.

Outputs: outputs/asr/maps/spectra_at_sites.{png,pdf}, spectra_at_sites.csv
"""
import argparse
import glob
import json
import os
import sys

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(__file__))
import diag_plots as dp  # noqa: E402
import run_month as rm  # noqa: E402

ROOT = rm.ROOT
SITES = [("57.51°N, 78.14°E", 57.51, 78.14), ("30.52°S, 132.46°W", -30.52, -132.46)]


def ang_dist(la0, lo0, la, lo):
    la0, lo0, la, lo = map(np.radians, (la0, lo0, la, lo))
    c = np.sin(la0) * np.sin(la) + np.cos(la0) * np.cos(la) * np.cos(lo - lo0)
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mapping_set", default=os.path.join(ROOT, "outputs", "asr", "maps", "mapping_set_v5.parquet"))
    ap.add_argument("--months_root", default=os.path.join(ROOT, "outputs", "asr", "months_v5"))
    ap.add_argument("--radius", type=float, default=0.5)
    ap.add_argument("--n", type=int, default=2)
    ap.add_argument("--out", default=os.path.join(ROOT, "outputs", "asr", "maps", "spectra_at_sites"))
    a = ap.parse_args()
    ms = pd.read_parquet(a.mapping_set)
    picks, summary = [], []
    for name, la, lo in SITES:
        d = ang_dist(la, lo, ms.clat.values, ms.clon.values)
        near = ms[d < a.radius]
        w = 1 / near.sd_r_Mg ** 2, 1 / near.sd_r_Al ** 2
        summary.append({"site": name, "n_quality0_within_radius": len(near), "radius_deg": a.radius,
                        "Mg_Si_wmean": float(np.sum(w[0] * near.r_Mg) / np.sum(w[0])),
                        "Mg_Si_se": float(1 / np.sqrt(np.sum(w[0]))),
                        "Al_Si_wmean": float(np.sum(w[1] * near.r_Al) / np.sum(w[1])),
                        "Al_Si_se": float(1 / np.sqrt(np.sum(w[1])))})
        for _, row in near.sort_values("n_counts_fit", ascending=False).head(a.n).iterrows():
            picks.append((name, row))
    pd.DataFrame(summary).to_csv(a.out + ".csv", index=False)
    print(pd.DataFrame(summary).round(3).to_string(index=False))

    fig, axes = plt.subplots(2 * len(picks), 1, figsize=(7.5, 2.6 * len(picks)),
                             gridspec_kw={"height_ratios": [3, 1] * len(picks)})
    for i, (name, row) in enumerate(picks):
        m = row["month"]
        mdir = os.path.join(a.months_root, m)
        calib = json.load(open(os.path.join(mdir, "calib.json")))
        bkg = np.load(os.path.join(mdir, "background.npy"))
        inst = np.load(os.path.join(mdir, "instrumental.npy"))
        rm._init(None, calib, bkg, inst, True, calib.get("response", "arf"), tuple(calib.get("hi_band", (8, 17))))
        path = glob.glob(os.path.join(ROOT, "data", "raw", *m.split("_"), "**", row["source_file"]), recursive=True)[0]
        r = rm.read_one(path)
        F, _ = rm._fitter_for(r["temp"])
        o = F.fit(r["counts"], r["exposure"], r["n_scd"], r["solarang"], T_mk=row["goes_T_mk"])
        t = "(%s) %s, %s" % ("abcdefgh"[i], name, row["source_file"][11:26])
        dp.plot_fit(F, r, o, axes[2 * i], axes[2 * i + 1], t)
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="upper center", ncol=len(lab), fontsize=7, frameon=False)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    for ext in ("png", "pdf"):
        fig.savefig(f"{a.out}.{ext}", dpi=150)
    print(a.out)


if __name__ == "__main__":
    main()
