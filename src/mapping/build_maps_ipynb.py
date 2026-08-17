"""
Faithful reproduction of src/notebooks/map.ipynb's mapping recipe -- the
method that produced the clean, high-resolution Al/Si map the user liked --
applied to the expanded fit result.

Differences from build_maps.py (which used map_making_v1.GaussianArray at
1 deg on raw concentration): here we map the ELEMENT/Si RATIO at 1080x540
(~0.33 deg), rasterize by center-point binning + gaussian_filter smoothing
+ 3x3 neighbour infill, and tail-cut with log1p+IQR plus an optional hard
percentile cap -- exactly as the notebook does.
"""
import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.ndimage import gaussian_filter, convolve

BACKGROUND_COLOR = "#0d0d0d"

# Fixed display ranges matched to outputs/figures/legacy_ratio_maps/conc_<El>_Si.png
# colorbars. Confirmed by direct comparison: Fe/Ca/Ti/O's ratio-value
# distribution on the new, much larger fit sits in the same regime as the
# legacy maps, so using the identical fixed vmin/vmax reproduces their look
# almost exactly. Al/Mg's distribution sits on a different absolute scale in
# the new fit (larger sample, different composition of months) -- forcing
# the legacy numbers there washes the map out, so those two fall back to
# percentile-based stretch (see --vmin_pct/--vmax_pct/--stretch_maxlat).
LEGACY_SCALE = {
    "Fe": (0.013, 0.048),
    "Ca": (0.004, 0.024),
    "Ti": (0.001, 0.013),
    "O": (55.0, 230.0),
}


def log_iqr_filter(df, col, k=1.5):
    logv = np.log1p(df[col].values)
    q1, q3 = np.percentile(logv, 25), np.percentile(logv, 75)
    iqr = q3 - q1
    lo, hi = q1 - k * iqr, q3 + k * iqr
    return df[(logv >= lo) & (logv <= hi)].copy()


def rasterize_footprints(clat, clon, values, width, height, sigma=1.1):
    """Vectorised center-point binning (same result as the notebook's
    per-row loop, but via np.add.at) + gaussian smoothing + renormalise."""
    clon = ((clon + 180.0) % 360.0) - 180.0
    ci = np.clip(((clon + 180.0) / 360.0 * width).astype(int), 0, width - 1)
    ri = np.clip(((clat + 90.0) / 180.0 * height).astype(int), 0, height - 1)
    value_grid = np.zeros((height, width), dtype=np.float64)
    weight_grid = np.zeros((height, width), dtype=np.float64)
    np.add.at(value_grid, (ri, ci), values)
    np.add.at(weight_grid, (ri, ci), 1.0)
    cov = weight_grid > 0
    avg = np.where(cov, value_grid / np.where(weight_grid > 0, weight_grid, 1), np.nan)
    coverage_pct = cov.sum() / (width * height) * 100.0
    vals = gaussian_filter(np.nan_to_num(avg, nan=0.0), sigma=sigma)
    wts = gaussian_filter(cov.astype(float), sigma=sigma)
    final = np.where(cov, vals / np.where(wts > 1e-6, wts, 1.0), np.nan)
    return final, coverage_pct


def fill_empty_pixels(grid):
    filled = grid.copy()
    empty = np.isnan(grid)
    kernel = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]])
    safe = np.nan_to_num(grid, nan=0.0)
    valid = (~empty).astype(float)
    ncount = convolve(valid, kernel, mode="constant", cval=0.0)
    nsum = convolve(safe, kernel, mode="constant", cval=0.0)
    loc = empty & (ncount > 0)
    filled[loc] = nsum[loc] / ncount[loc]
    return filled


def render(grid, coverage_pct, title, out_path, dpi=200, vmin_pct=2, vmax_pct=98,
           stretch_maxlat=None, vmin=None, vmax=None):
    if vmin is not None and vmax is not None:
        # Explicit fixed scale (e.g. matched to a reference/legacy map).
        pass
    else:
        # Compute the display stretch from a chosen latitude band so the
        # bright polar grazing-geometry tail does not compress the
        # equatorial bulk into the dark end (affects ONLY the PNG colour
        # scale, not the data).
        ref = grid
        if stretch_maxlat is not None:
            h = grid.shape[0]
            lat = (np.arange(h) / h) * 180 - 90
            band = np.abs(lat) <= stretch_maxlat
            ref = grid[band, :]
        vmin, vmax = np.nanpercentile(ref, vmin_pct), np.nanpercentile(ref, vmax_pct)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap("viridis")
    normed = norm(np.nan_to_num(grid, nan=vmin))
    rgba = cmap(normed)
    rgba[..., 3] = np.where(np.isnan(grid), 0.0, 1.0)

    fw = 18
    fh = fw / 2 + 1.5
    fig, ax = plt.subplots(figsize=(fw, fh), dpi=dpi)
    fig.patch.set_facecolor(BACKGROUND_COLOR)
    ax.set_facecolor(BACKGROUND_COLOR)
    ax.imshow(rgba, extent=[-180, 180, -90, 90], origin="lower", aspect="equal",
              interpolation="bilinear")
    ax.set_xlim(-180, 180); ax.set_ylim(-90, 90)
    ax.set_xticks(range(-180, 181, 30)); ax.set_yticks(range(-90, 91, 30))
    ax.tick_params(colors="white", labelsize=7)
    ax.grid(color="white", linewidth=0.3, alpha=0.25, linestyle="--")
    for sp in ax.spines.values():
        sp.set_edgecolor("white"); sp.set_linewidth(0.5)
    ax.set_xlabel("Longitude (deg)", color="white", fontsize=9)
    ax.set_ylabel("Latitude (deg)", color="white", fontsize=9)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.025, pad=0.08, aspect=50)
    cbar.set_label(title, color="white", fontsize=10)
    cbar.ax.tick_params(colors="white", labelsize=8)
    cbar.outline.set_edgecolor("white")
    ax.set_title(f"Lunar {title} Map", color="white", fontsize=13, pad=10, fontweight="bold")
    ax.text(0.99, 0.02, f"Map Coverage: {coverage_pct:.2f}%", transform=ax.transAxes,
            color="white", fontsize=8, ha="right", va="bottom",
            bbox=dict(boxstyle="round,pad=0.3", fc=BACKGROUND_COLOR, ec="white", alpha=0.6, lw=0.5))
    plt.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit_result", default="../../data/processed/catalogue_expanded_result.parquet")
    ap.add_argument("--element", default="Al")
    ap.add_argument("--out_dir", default="../../outputs/figures/maps_ipynb")
    ap.add_argument("--grid", type=int, default=1080, help="grid width (height = grid/2)")
    ap.add_argument("--ratio", action="store_true", default=True)
    ap.add_argument("--no_ratio", dest="ratio", action="store_false")
    ap.add_argument("--flux_top_pct", type=float, default=0.0,
                    help="keep top X%% by flux_<element>; 0 disables (use all fit data, the default)")
    ap.add_argument("--cap_pct", type=float, default=99.0,
                    help="hard upper percentile cap on the mapped value after log-IQR (notebook's manual x); 0 disables")
    ap.add_argument("--save_npy", action="store_true", help="also save the raster as .npy for SR")
    ap.add_argument("--max_abs_lat", type=float, default=None,
                    help="drop footprints with |lat| above this (polar grazing-geometry cut for SR analysis)")
    ap.add_argument("--vmin_pct", type=float, default=2.0, help="display: lower percentile for colour stretch")
    ap.add_argument("--vmax_pct", type=float, default=98.0, help="display: upper percentile for colour stretch")
    ap.add_argument("--stretch_maxlat", type=float, default=70.0,
                    help="display: compute the colour stretch from |lat|<this band (0/None = whole map)")
    ap.add_argument("--vmin", type=float, default=None, help="display: explicit fixed lower bound (overrides vmin_pct)")
    ap.add_argument("--vmax", type=float, default=None, help="display: explicit fixed upper bound (overrides vmax_pct)")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    el = args.element
    df = pd.read_parquet(args.fit_result)
    df = df[df["success"] == True].copy()

    # value column: Element/Si ratio (default) or raw concentration
    if args.ratio and el != "Si":
        with np.errstate(divide="ignore", invalid="ignore"):
            df["val"] = df[f"conc_{el}"] / df["conc_Si"]
        col, label = "val", f"{el}/Si"
    else:
        df["val"] = df[f"conc_{el}"]
        col, label = "val", f"{el}"
    df = df[np.isfinite(df["val"]) & (df["val"] > 0)]
    n0 = len(df)

    if args.flux_top_pct and f"flux_{el}" in df.columns:
        thr = np.nanpercentile(df[f"flux_{el}"].values, 100 - args.flux_top_pct)
        df = df[df[f"flux_{el}"].values >= thr]
    n1 = len(df)

    df = log_iqr_filter(df, col)
    n2 = len(df)

    if args.cap_pct:
        cap = np.percentile(df[col].values, args.cap_pct)
        df = df[df[col].values <= cap]
    n3 = len(df)

    # vectorised footprint centers
    fps = df["footprint"].values
    clat = np.array([np.mean([c[0] for c in fp]) for fp in fps])
    clon = np.array([np.mean([c[1] for c in fp]) for fp in fps])
    values = df[col].values.astype(np.float64)

    if args.max_abs_lat is not None:
        keep = np.abs(clat) <= args.max_abs_lat
        clat, clon, values = clat[keep], clon[keep], values[keep]
    n4 = len(values)

    print(f"{el}: {n0} -> flux {n1} -> log-IQR {n2} -> cap {n3} -> lat {n4} spectra "
          f"({label}, grid {args.grid}x{args.grid//2})")

    grid, cov = rasterize_footprints(clat, clon, values, args.grid, args.grid // 2)
    grid = fill_empty_pixels(grid)
    out_png = os.path.join(args.out_dir, f"{el}_Si_map.png" if args.ratio else f"{el}_map.png")
    smaxlat = args.stretch_maxlat if args.stretch_maxlat and args.stretch_maxlat > 0 else None
    render(grid, cov, label, out_png, vmin_pct=args.vmin_pct, vmax_pct=args.vmax_pct,
           stretch_maxlat=smaxlat, vmin=args.vmin, vmax=args.vmax)
    print(f"coverage {cov:.2f}%  saved {out_png}")
    if args.save_npy:
        np.save(os.path.join(args.out_dir, f"{el}_abundance.npy"), grid)


if __name__ == "__main__":
    main()
