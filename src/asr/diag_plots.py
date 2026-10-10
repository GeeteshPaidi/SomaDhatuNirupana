"""
Diagnostic plots of real-spectrum fits (TRACKER 1.8): re-fits selected
spectra of one month and plots data, total model, background, line and
scatter components, and deviance residuals.

Selection: high / median / low total counts in the fit window, or by
explicit --files.
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
import fitter as ft  # noqa: E402
import run_month as rm  # noqa: E402
import xrf_model as xm  # noqa: E402

ROOT = rm.ROOT


def plot_fit(F, r, o, ax_top, ax_bot, title):
    s1 = np.sin(np.radians(90 - r["solarang"]))
    th = np.array([o[k] for k in ft.PARAMS])
    bkg = F.background_counts(r["exposure"], r["n_scd"])
    d = r["counts"][F.sel]
    m = F.model(th, s1, bkg)
    lines_only = F.source(np.r_[th[:6], 0.0], s1)
    sc = F.source(th, s1) - lines_only
    E = F.E
    ax_top.step(E, d, where="mid", color="0.5", lw=0.7, label="data")
    ax_top.plot(E, m, "k", lw=1.2, label="model")
    ax_top.plot(E, th[7] * bkg[0], "C0", lw=0.9, label="particle bkg")
    ax_top.plot(E, th[8] * bkg[1], "C1", lw=0.9, label="instr. (solar)")
    ax_top.plot(E, lines_only, "C3", lw=0.9, label="lines")
    ax_top.plot(E, sc, "C2", lw=0.9, label="scatter")
    ax_top.set_yscale("symlog", linthresh=1)
    ax_top.set_xlim(1.1, 7.6)
    ax_top.set_title(title, fontsize=8)
    ax_top.legend(fontsize=6, ncol=5)
    dev = np.sign(d - m) * np.sqrt(np.maximum(ft._cstat_terms(d, np.maximum(m, 1e-12)), 0))
    ax_bot.plot(E, dev, ".", ms=2, color="k")
    ax_bot.axhline(0, color="r", lw=0.5)
    ax_bot.set_xlim(1.1, 7.6)
    ax_bot.set_ylabel("dev. resid", fontsize=7)
    ax_bot.set_xlabel("Energy (keV)", fontsize=7)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--month", default="2024/12")
    ap.add_argument("--kT", type=float, default=None)
    ap.add_argument("--subdir", default=None, help="month output subdir (default <YYYY_MM>)")
    ap.add_argument("--months_root", default=None, help="default outputs/asr/months; v3: outputs/asr/months_v3")
    ap.add_argument("--gamma_zero", action="store_true")
    ap.add_argument("--n_per", type=int, default=2)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    tag = args.month.replace("/", "_")
    mdir = os.path.join(args.months_root or os.path.join(ROOT, "outputs", "asr", "months"), args.subdir or tag)
    calib = json.load(open(os.path.join(mdir, "calib.json")))
    bkg = np.load(os.path.join(mdir, "background.npy"))
    inst = np.load(os.path.join(mdir, "instrumental.npy"))
    df = pd.read_parquet(os.path.join(mdir, "fits.parquet"))
    df = df[df["quality"] & ft.Q_BADGEOM == 0].sort_values("n_counts_fit").reset_index(drop=True)
    n = len(df)
    picks = {"high": df.iloc[[n - 1 - i * 5 for i in range(args.n_per)]],
             "median": df.iloc[[n // 2 + i * 5 for i in range(args.n_per)]],
             "low": df.iloc[[int(0.05 * n) + i * 5 for i in range(args.n_per)]]}
    F = rm.make_fitter(calib, bkg, inst, args.kT, args.gamma_zero)
    rows = sum(len(v) for v in picks.values())
    fig, axes = plt.subplots(2 * rows, 1, figsize=(9, 3.2 * rows),
                             gridspec_kw={"height_ratios": [3, 1] * rows})
    i = 0
    for lvl, sub in picks.items():
        for _, row in sub.iterrows():
            path = glob.glob(os.path.join(ROOT, "data", "raw", *args.month.split("/"), "**", row["source_file"]),
                             recursive=True)[0]
            r = rm.read_one(path)
            o = F.fit(r["counts"], r["exposure"], r["n_scd"], r["solarang"],
                      T_mk=None if args.kT else row.get("goes_T_mk"))
            ttl = (f"[{lvl}] {row['source_file']}  counts={o['n_counts_fit']:.0f}  SOLARANG={r['solarang']:.1f}  "
                   f"nSCD={r['n_scd']}  T={o.get('T_mk', float('nan')):.1f}MK  z={o['gof_z']:.1f}  q={o['quality']}\n" +
                   "  ".join(f"{e}/Si={o['r_'+e]:.3f}±{o['sd_r_'+e]:.3f}" for e in ["Mg", "Al", "Ca"]) +
                   f"  B={o['B']:.2f} beta={o['beta']:.2f} gamma={o['gamma']:.2f}")
            plot_fit(F, r, o, axes[2 * i], axes[2 * i + 1], ttl)
            i += 1
    fig.tight_layout()
    out = args.out or os.path.join(mdir, "diag_fits.png")
    fig.savefig(out, dpi=80)
    print(out)


if __name__ == "__main__":
    main()
