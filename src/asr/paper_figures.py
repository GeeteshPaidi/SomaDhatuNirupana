"""
Figures for the ASR manuscript that are not produced by maps.py /
validate_asr.py / diag_plots.py. Output: paper_asr/figures/*.pdf (+ .png).

  fig_solar_spectra   CHIANTI isothermal spectra at 3, 5, 8, 15 MK with K edges
  fig_background      December 2024 particle (umbra) vs solar-driven instrumental templates,
                      and nightside Al Ka rate vs SOLARANG (sunlit / umbra)
  fig_goes            GOES temperature distribution of CLASS exposures by year,
                      and of the mapping set
"""
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
import solar  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
FIG = os.path.join(ROOT, "paper_asr", "figures")
# Primary product (revision 1): v5 = ARF response, TEMP-binned detector calibration, 8-17 keV proxy band.
PT = "_v5"
PM = "months_v5"
LOLA = os.path.join(ROOT, "data", "lola", "ldem_4.img")
# Major maria and basins (approximate centres, deg N / deg E) for map orientation.
FEATURES = [("Imb", 32.8, -15.6), ("Proc", 18.4, -57.4), ("Ser", 28.0, 17.5), ("Tra", 8.5, 31.4),
            ("Cri", 17.0, 59.1), ("Fec", -7.8, 51.3), ("Nub", -21.3, -16.6), ("Fri", 56.0, 1.4),
            ("Ori", -19.4, -92.8), ("SPA", -53.0, -169.0), ("Mos", 27.3, 147.9)]


def lola_1deg():
    """LOLA LDEM_4 (4 px/deg, 0-360 E, row 0 = 90 N) block-averaged to 1 deg on the
    -180..180 E, south-up-from-row-0 grid used by the maps; heights in km."""
    a = np.fromfile(LOLA, dtype="<i2").reshape(720, 1440).astype(float) * 0.5 / 1000
    a = a.reshape(180, 4, 360, 4).mean((1, 3))[::-1]
    return np.roll(a, 180, axis=1)


def overlay_topo(ax, fontsize=5.5):
    if not os.path.exists(LOLA):
        return
    from scipy.ndimage import gaussian_filter
    h = gaussian_filter(lola_1deg(), sigma=(3, 3), mode=("nearest", "wrap"))  # ~3 deg smoothing: basin scale
    lon = np.arange(-179.5, 180, 1.0)
    lat = np.arange(-89.5, 90, 1.0)
    ax.contour(lon, lat, h, levels=[-2.5, 0, 2.5], colors="w", linewidths=0.45, alpha=0.9)
    for name, la, lo in FEATURES:
        ax.text(lo, la, name, fontsize=fontsize, color="k", ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.1", fc="w", ec="none", alpha=0.6))
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
    # The right panel needs outputs/asr/bkg_check3_2024_12.parquet, built from raw nightside
    # Level-1 spectra that are no longer on disk. Without it, keep the frozen figure
    # (products/figures_frozen/fig_background.png, a 600 dpi raster of the figure region of the
    # first-submission PDF; rasterised so that no hidden page text is embedded) instead of
    # writing a half-empty one.
    if not os.path.exists(os.path.join(ROOT, "outputs", "asr", "bkg_check3_2024_12.parquet")):
        import shutil
        shutil.copy(os.path.join(ROOT, "products", "figures_frozen", "fig_background.png"),
                    os.path.join(FIG, "fig_background.png"))
        print("fig_background: diagnostic parquet missing -> frozen figure used")
        return
    mdir = os.path.join(ROOT, "outputs", "asr", "months", month)
    cal = json.load(open(os.path.join(mdir, "calib.json")))
    P = np.load(os.path.join(mdir, "background.npy"))
    T = np.load(os.path.join(mdir, "instrumental.npy"))
    E = cal["gain"] * np.arange(len(P)) + cal["offset"]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
    ax = axes[0]
    ax.plot(E, P, "C0", ls="-", lw=0.8, label="particle (umbra)")
    ax.plot(E, np.clip(T, 1e-5, None), "C1", ls="--", lw=0.8, label="solar-driven instrumental")
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
        for lab, m, c, mk, s in [("spacecraft sunlit", n.sc_sunlit, "C1", "o", 1), ("spacecraft in umbra", ~n.sc_sunlit, "k", "x", 3)]:
            ax.scatter(n.sa[m], n.al[m].clip(lower=0.01), s=s, alpha=0.3, color=c, marker=mk, lw=0.4,
                       label=lab, rasterized=True)
        ax.axvline(n.shadow_ang.median(), color="k", ls="--", lw=0.7)
        ax.set_yscale("log")
        ax.set_xlabel("SOLARANG (deg)")
        ax.set_ylabel(r"Al K$\alpha$ excess rate (counts s$^{-1}$)")
        ax.legend(fontsize=7, markerscale=6)
    fig.tight_layout()
    save(fig, "fig_background")


def fig_goes():
    ms_path = os.path.join(ROOT, "outputs", "asr", "maps", f"mapping_set{PT}.parquet")
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


if __name__ == "__main__" and len(sys.argv) == 1:
    fig_solar_spectra()
    fig_background()
    fig_goes()
    print("figures ->", FIG)


def fig_maps(res=2, tag=PT):
    d = np.load(os.path.join(ROOT, "outputs", "asr", "maps", f"maps_{res}deg{tag}.npz"))
    ext = [-180, 180, -90, 90]
    fig, axes = plt.subplots(3, 2, figsize=(7.2, 7.0))
    for j, (e, lo, hi) in enumerate([("Mg", 0.1, 0.5), ("Al", 0.2, 1.3)]):
        for i, (key, cmap, lab, vlims) in enumerate([
                (f"{e}_mean", "viridis", f"{e}/Si (weight ratio)", (lo, hi)),
                (f"{e}_n", "magma", "log$_{10}$ footprints", None),
                (f"{e}_se", "cividis", f"SE of {e}/Si", None)]):
            ax = axes[i, j]
            a = d[key]
            if key.endswith("_n"):
                a = np.log10(np.maximum(a, 0.1))
                vlims = (0, np.nanpercentile(a, 99))
            if key.endswith("_se"):
                vlims = (0, np.nanpercentile(a, 95))
            im = ax.imshow(a, origin="lower", extent=ext, cmap=cmap, vmin=vlims[0], vmax=vlims[1], aspect="auto")
            cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
            cb.set_label(lab, fontsize=7)
            cb.ax.tick_params(labelsize=6)
            ax.tick_params(labelsize=6)
            ax.set_xticks([-180, -90, 0, 90, 180])
            ax.set_yticks([-90, -45, 0, 45, 90])
            if i == 0:
                overlay_topo(ax)
            if i == 2:
                ax.set_xlabel("Longitude (\u00b0E)", fontsize=7)
            if j == 0:
                ax.set_ylabel("Latitude (\u00b0)", fontsize=7)
    fig.tight_layout()
    save(fig, "fig_maps")


def fig_validation():
    import shutil
    src = os.path.join(ROOT, "outputs", "asr", "validation", f"validation_scatter{PT}.pdf")
    os.makedirs(FIG, exist_ok=True)
    shutil.copy(src, os.path.join(FIG, "fig_validation.pdf"))


def fig_fits():
    import subprocess
    subprocess.run([sys.executable, os.path.join(ROOT, "src", "asr", "diag_plots.py"), "--n_per", "1",
                    "--months_root", os.path.join(ROOT, "outputs", "asr", PM), "--gamma_zero",
                    "--out", os.path.join(FIG, "fig_fits.pdf"), "--paper"], check=True)


def fig_wmedian(tag=PT):
    """Sensitivity map: inverse-variance-weighted median per 2-deg cell, assigning each
    footprint to the cell of its centre (the primary maps share each footprint over
    a 5x5 sample grid). Same colour scales as fig_maps."""
    import alt_median_map as amm
    ms = pd.read_parquet(os.path.join(ROOT, "outputs", "asr", "maps", f"mapping_set{tag}.parquet"))
    ri = np.clip(((ms.clat.values + 90) / 2).astype(int), 0, 89)
    ci = np.clip(((ms.clon.values + 180) / 2).astype(int), 0, 179)
    cell = ri * 180 + ci
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.5))
    for ax, (e, lo, hi) in zip(axes, [("Mg", 0.1, 0.5), ("Al", 0.2, 1.3)]):
        w = 1 / np.maximum(ms[f"sd_r_{e}"].values, 1e-6) ** 2
        v = amm.group_stat(cell, ms[f"r_{e}"].values, w, "median")
        a = np.full(90 * 180, np.nan)
        a[v.index.values.astype(int)] = v.values
        im = ax.imshow(a.reshape(90, 180), origin="lower", extent=[-180, 180, -90, 90], cmap="viridis",
                       vmin=lo, vmax=hi, aspect="auto")
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
        cb.set_label(f"{e}/Si (weighted median)", fontsize=7)
        cb.ax.tick_params(labelsize=6)
        ax.set_xticks([-180, -90, 0, 90, 180])
        ax.set_yticks([-90, -45, 0, 45, 90])
        ax.tick_params(labelsize=6)
        ax.set_xlabel("Longitude (°E)", fontsize=7)
    axes[0].set_ylabel("Latitude (°)", fontsize=7)
    fig.tight_layout()
    save(fig, f"fig_maps_wmedian_posthoc{tag}")


def fig_xsm(months=("2024_12", "2025_06", "2025_12")):
    X = os.path.join(ROOT, "outputs", "asr", "xsm")
    j = pd.concat([pd.read_parquet(os.path.join(X, f"class_refits_{m}.parquet")) for m in months])
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.4))
    ax = axes[0]
    ax.scatter(j.goes_T_mk, j.T1_xsm, s=2, alpha=0.3, color="k", rasterized=True)
    lim = [2, max(j.goes_T_mk.max(), j.T1_xsm.max()) * 1.05]
    ax.plot(lim, lim, "r--", lw=0.8)
    ax.set_xlim(lim), ax.set_ylim(lim)
    ax.set_xlabel("GOES T (MK)"), ax.set_ylabel("XSM isothermal T (MK)")
    for ax, e in zip(axes[1:], ["Mg", "Al"]):
        bins = np.linspace(-60, 60, 61)
        for n, c, ls, lab in [("xsm1t", "C0", "-", "XSM 1T"), ("xsm2t", "C3", "--", "XSM 2T")]:
            ok = (j.quality_goes == 0) & (j[f"quality_{n}"] == 0)
            ax.hist(100 * (j.loc[ok, f"r_{e}_{n}"] / j.loc[ok, f"r_{e}_goes"] - 1), bins=bins,
                    histtype="step", color=c, ls=ls, lw=1.1, label=lab)
        ax.axvline(0, color="k", lw=0.6)
        ax.set_xlabel(f"{e}/Si change vs GOES (%)")
        ax.legend(fontsize=6)
    axes[1].set_ylabel("spectra")
    fig.tight_layout()
    save(fig, "fig_xsm")


def fig_fwhm_temp():
    """R1.5: Al Ka FWHM vs detector TEMP, one point per TEMP bin of each month's calibration."""
    rows = []
    for cj in sorted(glob.glob(os.path.join(ROOT, "outputs", "asr", PM, "*", "calib.json"))):
        tag = os.path.basename(os.path.dirname(cj))
        if len(tag) != 7:
            continue
        for b in json.load(open(cj)).get("temp_cal") or []:
            rows.append({"year": int(tag[:4]), "t": b["t_med"], "fwhm": b["fwhm_1.49keV_eV"]})
    d = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(3.6, 2.8))
    for y, g in d.groupby("year"):
        ax.plot(g.t, g.fwhm, "o", ms=3, alpha=0.8, label=str(y))
    ax.set_xlabel("Detector temperature TEMP (°C)", fontsize=8)
    ax.set_ylabel("Al Kα FWHM (eV)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=6, ncol=2)
    fig.tight_layout()
    save(fig, "fig_fwhm_temp")


def fig_epochs():
    import shutil
    shutil.copy(os.path.join(ROOT, "outputs", "asr", "epochs", "epoch_compare.pdf"), os.path.join(FIG, "fig_epochs.pdf"))


def fig_spectra_sites():
    import shutil
    shutil.copy(os.path.join(ROOT, "outputs", "asr", "maps", "spectra_at_sites.pdf"),
                os.path.join(FIG, "fig_spectra_sites.pdf"))


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "results":
    fig_maps()
    fig_wmedian()
    fig_epochs()
    fig_spectra_sites()
    fig_validation()
    fig_fits()
    fig_xsm()
    print("result figures ->", FIG)
