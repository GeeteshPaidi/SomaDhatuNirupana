"""
Builds elemental abundance maps from fitted concentrations (Sec 5 of the paper).

Data array: uses map_making_v1.GaussianArray - footprint-quadrilateral
rasterization + adaptive 2D Gaussian weighting with the exact
target_diagonal=17.625, base_value=2.1739 constants Sec 5 quotes as
d_target/sigma_base. This is what the paper's equations actually describe,
so the figure and the text stay consistent.

Rendering: percentile-stretched (2-98%) RGBA overlay on a dark background,
matching the look from src/notebooks/map.ipynb (which produced the maps
the user liked), but applied to the array from the method above rather than
the notebook's simpler center-point-binning rasterizer.
"""
import argparse
import os

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from map_making_v1 import GaussianArray

ELEMENTS = ["Fe", "Al", "Mg", "Si", "Ca", "Ti", "O"]
BLOCK_LAT = [-90, 90, 90, -90]
BLOCK_LON = [-180, -180, 180, 180]


def log_iqr_filter(values, k=1.5):
    """log1p + IQR outlier removal (matches src/notebooks/map.ipynb).

    Cuts the long high-value tail so the color scale reflects the bulk of the
    distribution rather than a handful of extreme fits. Returns a boolean mask
    of rows to KEEP.
    """
    v = np.asarray(values, dtype=np.float64)
    finite = np.isfinite(v) & (v > 0)
    logv = np.full_like(v, np.nan)
    logv[finite] = np.log1p(v[finite])
    q1 = np.nanpercentile(logv, 25)
    q3 = np.nanpercentile(logv, 75)
    iqr = q3 - q1
    lo, hi = q1 - k * iqr, q3 + k * iqr
    return finite & (logv >= lo) & (logv <= hi)


def build_abundance_grid(df, column, grid_size=(180, 360)):
    grid = GaussianArray(grid_size=grid_size)
    n_used = 0
    for _, row in df.iterrows():
        fp = row["footprint"]
        img_lat = [p[0] for p in fp]
        img_lon = [((p[1] + 180) % 360) - 180 for p in fp]
        if max(img_lon) - min(img_lon) > 180:
            continue
        val = row[column]
        if pd.isna(val) or val <= 0:
            continue
        grid.add_gaussian_box(img_lat, img_lon, BLOCK_LAT, BLOCK_LON, val)
        n_used += 1
    coverage = grid.check_coverage()
    return grid, n_used, coverage


def render_dark_map(arr_2d, title, out_path, cmap_name="viridis"):
    height, width = arr_2d.shape
    valid = arr_2d[~np.isnan(arr_2d)]
    if len(valid) == 0:
        raise ValueError(f"No valid pixels to render for {title}")
    vmin, vmax = np.nanpercentile(arr_2d, 2), np.nanpercentile(arr_2d, 98)
    if vmax <= vmin:
        vmax = vmin + 1e-9
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap(cmap_name)

    normed = norm(np.nan_to_num(arr_2d, nan=vmin))
    rgba = cmap(normed)
    rgba[..., 3] = np.where(np.isnan(arr_2d), 0.0, 1.0)

    fig_w = 12
    fig_h = fig_w / 2 + 1.3
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=300)
    fig.patch.set_facecolor("#0d0d0d")
    ax.set_facecolor("#0d0d0d")

    extent = [-180, 180, -90, 90]
    ax.imshow(rgba, extent=extent, origin="lower", aspect="equal", interpolation="bilinear")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xticks(range(-180, 181, 30))
    ax.set_yticks(range(-90, 91, 30))
    ax.tick_params(colors="white", labelsize=7)
    ax.grid(color="white", linewidth=0.3, alpha=0.25, linestyle="--")
    for spine in ax.spines.values():
        spine.set_edgecolor("white")
        spine.set_linewidth(0.5)
    ax.set_xlabel("Longitude (deg)", color="white", fontsize=9)
    ax.set_ylabel("Latitude (deg)", color="white", fontsize=9)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.025, pad=0.1, aspect=50)
    cbar.set_label(title, color="white", fontsize=10)
    cbar.ax.tick_params(colors="white", labelsize=8)
    cbar.outline.set_edgecolor("white")
    ax.set_title(title, color="white", fontsize=13, pad=10, fontweight="bold")

    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit_result", default="../../data/processed/catalogue_full_result.parquet")
    ap.add_argument("--out_dir", default="../../outputs/figures/maps")
    ap.add_argument("--grid_rows", type=int, default=180)
    ap.add_argument("--grid_cols", type=int, default=360)
    ap.add_argument("--cost_percentile_cutoff", type=float, default=None,
                     help="Drop fits with cost above this percentile (quality control for poor TRF convergence)")
    ap.add_argument("--max_abs_lat", type=float, default=None,
                     help="Drop footprints with |lat| above this (avoids polar projection distortion / grazing geometry)")
    ap.add_argument("--ratio", action="store_true",
                     help="Map Element/Si ratio instead of raw relative concentration (as in map.ipynb)")
    ap.add_argument("--log_iqr", action="store_true",
                     help="Apply log1p+IQR outlier removal per element before rasterizing (tail-cut, as in map.ipynb)")
    ap.add_argument("--flux_top_pct", type=float, default=None,
                     help="Per element, keep only rows in the top X%% by flux_<element> "
                          "(needs flux_<el> columns in the fit result, e.g. 50 for top half). "
                          "This is the user's top-Al-flux quality selection generalised per element.")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    df = pd.read_parquet(args.fit_result)
    df = df[df["success"] == True].copy()
    print(f"Loaded {len(df)} successful fits")

    if args.cost_percentile_cutoff is not None:
        cutoff = np.percentile(df["cost"], args.cost_percentile_cutoff)
        n_before = len(df)
        df = df[df["cost"] < cutoff]
        print(f"Cost QC filter (<{args.cost_percentile_cutoff}th pct = {cutoff:.2f}): {n_before} -> {len(df)}")

    if args.max_abs_lat is not None:
        lat = df["footprint"].apply(lambda pts: np.mean([p[0] for p in pts]))
        n_before = len(df)
        df = df[lat.abs() < args.max_abs_lat]
        print(f"Latitude QC filter (|lat|<{args.max_abs_lat}): {n_before} -> {len(df)}")

    coverage_report = []
    for el in ELEMENTS:
        col = f"conc_{el}"
        df_el = df

        # Per-element top-flux quality selection (the user's top-50%-flux idea).
        flux_col = f"flux_{el}"
        if args.flux_top_pct is not None and flux_col in df.columns:
            thr = np.nanpercentile(df[flux_col], 100 - args.flux_top_pct)
            n_before = len(df_el)
            df_el = df_el[df_el[flux_col] >= thr]
            print(f"{el}: flux top-{args.flux_top_pct:.0f}% (flux_{el}>={thr:.4g}): {n_before} -> {len(df_el)}")

        # Ratio to Si, if requested (Element/Si, as in map.ipynb).
        value_col = col
        if args.ratio and el != "Si":
            df_el = df_el.copy()
            with np.errstate(divide="ignore", invalid="ignore"):
                df_el["ratio_val"] = df_el[col] / df_el["conc_Si"]
            df_el = df_el[np.isfinite(df_el["ratio_val"])]
            value_col = "ratio_val"

        # Log-IQR tail-cut on the value distribution (as in map.ipynb).
        if args.log_iqr:
            keep = log_iqr_filter(df_el[value_col].values)
            n_before = len(df_el)
            df_el = df_el[keep]
            print(f"{el}: log-IQR tail-cut: {n_before} -> {len(df_el)}")

        grid, n_used, coverage = build_abundance_grid(df_el, value_col, grid_size=(args.grid_rows, args.grid_cols))
        arr = grid.arr[:, :, 0].copy()
        counts = grid.arr[:, :, 1]
        arr[counts == 0] = np.nan

        np.save(os.path.join(args.out_dir, f"{el}_abundance.npy"), arr)
        render_dark_map(arr, f"Lunar {el} Relative Abundance", os.path.join(args.out_dir, f"{el}_map.png"))
        print(f"{el}: used {n_used} footprints, coverage {coverage*100:.2f}%")
        coverage_report.append({"element": el, "n_footprints": n_used, "coverage_pct": coverage * 100})

    pd.DataFrame(coverage_report).to_csv(os.path.join(args.out_dir, "coverage_report.csv"), index=False)
    print("Saved coverage_report.csv")


if __name__ == "__main__":
    main()
