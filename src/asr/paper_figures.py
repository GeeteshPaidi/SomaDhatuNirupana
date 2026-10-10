"""
Figures for the ASR manuscript that are not produced by maps.py /
validate_asr.py / diag_plots.py. Output: paper_asr/figures/*.pdf (+ .png).

  fig_solar_spectra   CHIANTI isothermal spectra at 3, 5, 8, 15 MK with K edges
  fig_background      December 2024 particle (umbra) vs solar-driven instrumental templates,
                      and nightside Al Ka rate vs SOLARANG (sunlit / umbra)
  fig_goes            GOES temperature distribution of CLASS exposures by year,
                      and of the mapping set
"""
import json
import os
import sys

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(__file__))
import solar  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
FIG = os.path.join(ROOT, "paper_asr", "figures")
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8})
EDGES = {"Mg": 1.305, "Al": 1.560, "Si": 1.839, "Ca": 4.038, "Ti": 4.966, "Fe": 7.112}


def save(fig, name):
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, name + ".pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(FIG, name + ".png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_solar_spectra():
    st = solar.SolarTable()
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    for T, c in zip([3, 5, 8, 15], ["C0", "C1", "C2", "C3"]):
        ax.plot(st.e_mid, st.spectrum(T), color=c, lw=0.7, label=f"{T} MK")
    for el, e in EDGES.items():
        ax.axvline(e, color="k", ls=":", lw=0.6)
        ax.text(e, 1.5e8, el, fontsize=7, ha="left")
    ax.set_yscale("log")
    ax.set_ylim(1e-2, 1e9)
    ax.set_xlim(1.0, 9.0)
    ax.set_xlabel("Energy (keV)")
    ax.set_ylabel(r"Photon flux (ph s$^{-1}$ cm$^{-2}$ keV$^{-1}$, EM=10$^{48}$ cm$^{-3}$)")
    ax.legend(fontsize=7, title="Isothermal T", title_fontsize=7)
    save(fig, "fig_solar_spectra")


def fig_background(month="2024_12"):
    mdir = os.path.join(ROOT, "outputs", "asr", "months", month)
    cal = json.load(open(os.path.join(mdir, "calib.json")))
    P = np.load(os.path.join(mdir, "background.npy"))
    T = np.load(os.path.join(mdir, "instrumental.npy"))
    E = cal["gain"] * np.arange(len(P)) + cal["offset"]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
    ax = axes[0]
    ax.plot(E, P, lw=0.8, label="particle (umbra)")
    ax.plot(E, np.clip(T, 1e-5, None), lw=0.8, label="solar-driven instrumental")
    ax.set_yscale("log")
    ax.set_xlim(0.8, 10)
    ax.set_ylim(1e-3, 30)
    ax.axvline(1.4867, color="k", ls=":", lw=0.6)
    ax.text(1.52, 12, "Al K$\\alpha$", fontsize=7)
    ax.set_xlabel("Energy (keV)")
    ax.set_ylabel(r"counts s$^{-1}$ channel$^{-1}$ (16 SCD)")
    ax.legend(fontsize=7)
    f = os.path.join(ROOT, "outputs", "asr", "bkg_check3_2024_12.parquet")
    if os.path.exists(f):
        d = pd.read_parquet(f)
        n = d[d.sa >= 90]
        ax = axes[1]
        for lab, m, c in [("spacecraft sunlit", n.sc_sunlit, "C1"), ("spacecraft in umbra", ~n.sc_sunlit, "C0")]:
            ax.scatter(n.sa[m], n.al[m].clip(lower=0.01), s=1, alpha=0.3, color=c, label=lab, rasterized=True)
        ax.axvline(n.shadow_ang.median(), color="k", ls="--", lw=0.7)
        ax.set_yscale("log")
        ax.set_xlabel("SOLARANG (deg)")
        ax.set_ylabel(r"Al K$\alpha$ excess rate (counts s$^{-1}$)")
        ax.legend(fontsize=7, markerscale=6)
    fig.tight_layout()
    save(fig, "fig_background")


def fig_goes():
    ms_path = os.path.join(ROOT, "outputs", "asr", "maps", "mapping_set.parquet")
    tem = pd.read_parquet(solar.TEM)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
    ax = axes[0]
    tem = tem[tem.xrsb >= 1e-7]
    for y, g in tem.groupby(tem.index.year):
        ax.hist(g.temperature, bins=np.linspace(2, 20, 73), histtype="step", density=True, label=str(y))
    ax.set_xlabel("GOES isothermal temperature (MK), flux $\\geq$ B1")
    ax.set_ylabel("density")
    ax.legend(fontsize=7)
    if os.path.exists(ms_path):
        ms = pd.read_parquet(ms_path)
        ax = axes[1]
        ax.hist(ms.goes_T_mk, bins=np.linspace(2, 20, 73), color="0.4")
        ax.set_xlabel("GOES temperature of mapped exposures (MK)")
        ax.set_ylabel("exposures")
    fig.tight_layout()
    save(fig, "fig_goes")


if __name__ == "__main__":
    fig_solar_spectra()
    fig_background()
    fig_goes()
    print("figures ->", FIG)
