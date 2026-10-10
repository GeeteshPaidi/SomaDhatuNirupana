"""
Phase 3: element/Si maps with quality layers (TRACKER 3.1-3.3), plus the
Phase 2.1 pipeline-stage counts.

Mapping set (pre-declared, TRACKER decisions log): quality == 0 and GOES
XRS-B >= 1e-7 W/m^2 (B1). Each footprint quadrilateral is rasterised by
sampling a 5x5 bilinear grid of points inside its 4 corners; the footprint
contributes to every grid cell those points fall in (weight split equally
among its sample points). Per-cell value = inverse-variance-weighted mean
of the footprint ratios; quality layers = number of footprints, standard
error of the weighted mean, robust scatter. No smoothing or infill is used
for any statistic; an optional smoothed display layer is produced
separately and its infilled cells are masked explicitly.

Outputs: outputs/asr/maps/maps_<res>deg.npz, stage_counts.json,
footprint_size.json, figures under outputs/asr/maps/.
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.join(ROOT, "outputs", "asr", "maps")
R_MOON = 1737.4
ELS = ["Mg", "Al"]
Q_SOLAR_UNRELIABLE = 256
GOES_MIN_XRSB = 1e-7
PIN_EPS = 1e-3
COLS = ["source_file", "month", "lat", "lon", "solarang", "sat_alt", "n_scd", "exposure", "goes_T_mk",
        "goes_xrsb", "quality", "n_counts_fit", "gof_z", "success"] + \
       [f"r_{e}" for e in ELS] + [f"sd_r_{e}" for e in ELS] + ["r_Ca", "sd_r_Ca", "gamma", "A"]


def load_fits(months_dir=None, min_counts=300, min_xrsb=GOES_MIN_XRSB):
    months_dir = months_dir or os.path.join(ROOT, "outputs", "asr", "months")
    frames = []
    for d in sorted(glob.glob(os.path.join(months_dir, "*"))):
        f = os.path.join(d, "fits.parquet")
        tag = os.path.basename(d)
        if not os.path.exists(f) or not (len(tag) == 7 and tag[4] == "_"):
            continue  # skip sensitivity sub-runs (e.g. 2024_12_Tx1.2)
        df = pd.read_parquet(f)
        frames.append(df[[c for c in COLS if c in df.columns]])
    df = pd.concat(frames, ignore_index=True)
    q = df["quality"].fillna(-1).astype(int).values
    unreliable = ~(df["goes_xrsb"].values >= min_xrsb)
    q = np.where(unreliable & (q >= 0), q | Q_SOLAR_UNRELIABLE, q)
    # Boundary bug fix (found after the first validation, logged in TRACKER): a ratio
    # pinned within 1e-3 of zero has no valid Fisher uncertainty (sd collapses), but the
    # fitter's at-bound test used 1e-6. Treat r < 1e-3 (or > 4.999) for Mg/Si or Al/Si
    # as at bound (quality bit 8).
    pinned = np.zeros(len(df), bool)
    for e in ELS:
        r = df[f"r_{e}"].values
        pinned |= np.isfinite(r) & ((r < PIN_EPS) | (r > 4.999))
    q = np.where(pinned & (q >= 0), q | 8, q)
    df["quality_final"] = q
    df["use"] = (q == 0) & (df["n_counts_fit"].values >= min_counts)
    lat = np.vstack(df["lat"].values)
    lon = np.vstack(df["lon"].values)
    df["clat"] = lat.mean(1)
    # circular mean of corner longitudes
    df["clon"] = np.degrees(np.arctan2(np.sin(np.radians(lon)).mean(1), np.cos(np.radians(lon)).mean(1)))
    return df


def stage_counts(df, n_raw_files):
    return {
        "raw_L1_files": int(n_raw_files),
        "dayside_valid_geometry": int(len(df)),
        "goes_T_reliable_B1": int(np.sum(df["goes_xrsb"].values >= GOES_MIN_XRSB)),
        "fit_converged": int(np.sum(df["success"].fillna(False).values.astype(bool))),
        "passes_QC_quality0_and_B1": int(np.sum(df["quality_final"].values == 0)),
        "mapping_set": int(df["use"].sum()),
    }


def footprint_samples(lat, lon, n=5):
    """Bilinear sample points inside each quadrilateral (corners V0..V3 in order)."""
    u = (np.arange(n) + 0.5) / n
    U, V = np.meshgrid(u, u)
    U, V = U.ravel(), V.ravel()
    # unwrap longitudes around the first corner
    lon = lon.copy()
    lon = lon[:, :1] + ((lon - lon[:, :1] + 180) % 360 - 180)
    w = np.stack([(1 - U) * (1 - V), (1 - U) * V, U * V, U * (1 - V)], 1)  # V0,V1,V2,V3
    plat = lat @ w.T
    plon = ((lon @ w.T) + 180) % 360 - 180
    return plat, plon  # (n_fp, n*n)


def footprint_size_km(lat, lon):
    def gc(a1, o1, a2, o2):
        a1, o1, a2, o2 = map(np.radians, (a1, o1, a2, o2))
        return R_MOON * 2 * np.arcsin(np.sqrt(np.sin((a2 - a1) / 2) ** 2 +
                                              np.cos(a1) * np.cos(a2) * np.sin((o2 - o1) / 2) ** 2))
    s01 = gc(lat[:, 0], lon[:, 0], lat[:, 1], lon[:, 1])
    s12 = gc(lat[:, 1], lon[:, 1], lat[:, 2], lon[:, 2])
    return s01, s12


def grid_maps(df, res=2.0, n_sub=5):
    nlat, nlon = int(180 / res), int(360 / res)
    lat = np.vstack(df["lat"].values)
    lon = np.vstack(df["lon"].values)
    plat, plon = footprint_samples(lat, lon, n_sub)
    ri = np.clip(((plat + 90) / res).astype(int), 0, nlat - 1)
    ci = np.clip(((plon + 180) / res).astype(int), 0, nlon - 1)
    cell = ri * nlon + ci  # (n_fp, n_sub^2)
    frac = 1.0 / cell.shape[1]
    out = {"res": res, "nlat": nlat, "nlon": nlon}
    # coverage by all dayside footprints vs mapping set
    for name, mask in [("all", np.ones(len(df), bool)), ("used", df["use"].values)]:
        hit = np.zeros(nlat * nlon, bool)
        hit[np.unique(cell[mask])] = True
        out[f"hit_{name}"] = hit.reshape(nlat, nlon)
    m = df["use"].values
    cm = cell[m]
    for e in ELS:
        r = df[f"r_{e}"].values[m]
        sd = df[f"sd_r_{e}"].values[m]
        ok = np.isfinite(r) & np.isfinite(sd) & (sd > 0)
        w = np.where(ok, 1.0 / np.maximum(sd, 1e-6) ** 2, 0.0)
        r0 = np.where(ok, r, 0.0)
        W = np.bincount(cm.ravel(), weights=np.repeat(w * frac, cm.shape[1]), minlength=nlat * nlon)
        WR = np.bincount(cm.ravel(), weights=np.repeat(w * r0 * frac, cm.shape[1]), minlength=nlat * nlon)
        WR2 = np.bincount(cm.ravel(), weights=np.repeat(w * r0 ** 2 * frac, cm.shape[1]), minlength=nlat * nlon)
        N = np.bincount(cm.ravel(), weights=np.repeat(ok * frac, cm.shape[1]), minlength=nlat * nlon)
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = WR / W
            se = 1.0 / np.sqrt(W)
            wstd = np.sqrt(np.maximum(WR2 / W - mean ** 2, 0))
        mean[W == 0] = np.nan
        se[W == 0] = np.nan
        wstd[W == 0] = np.nan
        out[f"{e}_mean"] = mean.reshape(nlat, nlon)
        out[f"{e}_se"] = se.reshape(nlat, nlon)
        out[f"{e}_wstd"] = wstd.reshape(nlat, nlon)
        out[f"{e}_n"] = N.reshape(nlat, nlon)
    return out


def coverage_stats(maps):
    res = maps["res"]
    lat_c = -90 + res * (np.arange(maps["nlat"]) + 0.5)
    area = np.repeat(np.cos(np.radians(lat_c))[:, None], maps["nlon"], 1)
    s = {}
    for k in ["all", "used"]:
        h = maps[f"hit_{k}"]
        s[f"coverage_{k}_cells_pct"] = float(100 * h.mean())
        s[f"coverage_{k}_area_pct"] = float(100 * (area * h).sum() / area.sum())
    for e in ELS:
        n = maps[f"{e}_n"]
        s[f"{e}_cells_with_ge5_footprints_area_pct"] = float(100 * (area * (n >= 5)).sum() / area.sum())
    return s


def plot_maps(maps, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ext = [-180, 180, -90, 90]
    for e, (lo, hi) in [("Mg", (0.1, 0.5)), ("Al", (0.2, 1.6))]:
        fig, axes = plt.subplots(3, 1, figsize=(9, 12))
        im = axes[0].imshow(maps[f"{e}_mean"], origin="lower", extent=ext, cmap="viridis", vmin=lo, vmax=hi)
        fig.colorbar(im, ax=axes[0], label=f"{e}/Si (weight ratio)")
        axes[0].set_title(f"{e}/Si, inverse-variance mean per {maps['res']:.0f}° cell (observed cells only)")
        im = axes[1].imshow(np.log10(np.maximum(maps[f"{e}_n"], 0.1)), origin="lower", extent=ext, cmap="magma")
        fig.colorbar(im, ax=axes[1], label="log10 footprints per cell")
        axes[1].set_title("Quality layer: observation count")
        im = axes[2].imshow(maps[f"{e}_se"], origin="lower", extent=ext, cmap="cividis",
                            vmax=np.nanpercentile(maps[f"{e}_se"], 95))
        fig.colorbar(im, ax=axes[2], label=f"standard error of {e}/Si")
        axes[2].set_title("Quality layer: statistical uncertainty of the cell mean")
        for a in axes:
            a.set_xlabel("Longitude (°E)")
            a.set_ylabel("Latitude (°)")
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, f"map_{e}_Si_{maps['res']:.0f}deg.png"), dpi=110)
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", type=float, default=2.0)
    ap.add_argument("--min_counts", type=int, default=300)
    ap.add_argument("--min_xrsb", type=float, default=GOES_MIN_XRSB)
    ap.add_argument("--tag", default="")
    ap.add_argument("--months_dir", default=None, help="default outputs/asr/months (v2); v3: outputs/asr/months_v3")
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    df = load_fits(months_dir=args.months_dir, min_counts=args.min_counts, min_xrsb=args.min_xrsb)
    n_raw = 0
    for cj in glob.glob(os.path.join(args.months_dir or os.path.join(ROOT, "outputs", "asr", "months"), "*", "calib.json")):
        tag = os.path.basename(os.path.dirname(cj))
        if len(tag) == 7:
            n_raw += json.load(open(cj))["n_files"]
    stages = stage_counts(df, n_raw)
    lat, lon = np.vstack(df["lat"].values), np.vstack(df["lon"].values)
    a, b = footprint_size_km(lat, lon)
    fps = {"side01_km_median": float(np.median(a)), "side12_km_median": float(np.median(b)),
           "side01_km_p5_p95": [float(np.percentile(a, 5)), float(np.percentile(a, 95))],
           "side12_km_p5_p95": [float(np.percentile(b, 5)), float(np.percentile(b, 95))],
           "km_per_deg_equator": float(np.pi * R_MOON / 180)}
    maps = grid_maps(df, args.res)
    cov = coverage_stats(maps)
    tag = f"{args.res:.0f}deg{args.tag}"
    np.savez_compressed(os.path.join(OUT, f"maps_{tag}.npz"), **{k: v for k, v in maps.items()})
    summary = {"stages": stages, "footprint_size": fps, "coverage": cov,
               "mapping_set_ratio_percentiles": {e: [float(x) for x in np.nanpercentile(
                   df.loc[df.use, f"r_{e}"], [5, 25, 50, 75, 95])] for e in ELS + ["Ca"]},
               "params": vars(args)}
    with open(os.path.join(OUT, f"summary_{tag}.json"), "w") as f:
        json.dump(summary, f, indent=1)
    df.loc[df.use, ["source_file", "month", "clat", "clon", "solarang", "goes_T_mk", "goes_xrsb", "n_counts_fit"] +
           [f"r_{e}" for e in ELS] + [f"sd_r_{e}" for e in ELS]].to_parquet(
        os.path.join(OUT, f"mapping_set{args.tag}.parquet"), index=False)
    if not args.tag:
        plot_maps(maps, OUT)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
