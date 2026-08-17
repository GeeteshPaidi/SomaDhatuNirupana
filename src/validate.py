"""
Validation of fitted elemental abundances against Apollo/Luna landing-site
geochemistry (Sec 7.4 of the paper, Fig 4).

IMPORTANT CAVEAT (flagged, not silently worked around): the fitted
concentrations from src/fitting are unitless relative scale factors, not
wt% (see README "Known gaps"). A direct wt% comparison would therefore be
meaningless. Instead this script compares in RATIO SPACE (Element/Si),
which cancels the unknown per-element sensitivity constant to first order,
the same reasoning documented in docs/usage.md. Results should be reported
to the user as ratio-space validation, not absolute wt% validation, until
a real calibration step exists.

validation.xlsx layout (0-indexed openpyxl rows, from manual inspection):
  row 34 (1-indexed) header: 'Element', then site labels
  rows 35-40: Ti, Al, Fe, Mg, Ca, Si element wt% per site
Site label columns (1-indexed col B..T): 11,12,14,15a,15b,15c,15,16a,16b,16c,
  16,17a,17b,17c,17d,17,L16,L20,L24
"""
import argparse
import os

import numpy as np
import openpyxl
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

# Public Apollo/Luna landing site coordinates (deg, lon in [-180,180])
SITE_COORDS = {
    "11": (0.674, 23.473),
    "12": (-3.012, -23.42),
    "14": (-3.645, -17.471),
    "15": (26.132, 3.633),
    "16": (8.973, 15.499),
    "17": (20.19, 30.77),
    "L16": (-0.68, 56.3),
    "L20": (3.53, 56.55),
    "L24": (12.714, 62.2),
}

ELEMENT_ROWS = {"Ti": 33, "Al": 34, "Fe": 35, "Mg": 36, "Ca": 37, "Si": 38}  # openpyxl 1-indexed row
SITE_COLS = {  # 0-indexed column in the sheet
    "11": 1, "12": 2, "14": 3, "15a": 4, "15b": 5, "15c": 6, "15": 7,
    "16a": 8, "16b": 9, "16c": 10, "16": 11, "17a": 12, "17b": 13, "17c": 14,
    "17d": 15, "17": 16, "L16": 17, "L20": 18, "L24": 19,
}


def load_reference(xlsx_path):
    wb = openpyxl.load_workbook(xlsx_path, data_only=True)
    ws = wb["Sheet1"]
    ref = {}
    for site, coords in SITE_COORDS.items():
        col = SITE_COLS[site]
        wt = {el: ws.cell(row=r, column=col + 1).value for el, r in ELEMENT_ROWS.items()}
        if wt.get("Si") in (None, 0):
            continue
        ratios = {el: wt[el] / wt["Si"] for el in wt if el != "Si"}
        ref[site] = {"lat": coords[0], "lon": coords[1], "wt_pct": wt, "ratio_to_Si": ratios}
    return ref


def match_footprints(fit_df, site_lat, site_lon, radius_deg):
    fp = fit_df["footprint"]
    centers_lat = fp.apply(lambda pts: np.mean([p[0] for p in pts]))
    centers_lon = fp.apply(lambda pts: np.mean([p[1] for p in pts]))
    dlat = centers_lat - site_lat
    dlon = (centers_lon - site_lon + 180) % 360 - 180
    dist = np.sqrt(dlat**2 + (dlon * np.cos(np.radians(site_lat)))**2)
    return fit_df[dist.values <= radius_deg]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit_result", default="../data/processed/catalogue_full_result.parquet")
    ap.add_argument("--reference", default="../data/reference/validation.xlsx")
    ap.add_argument("--out_dir", default="../outputs/figures")
    ap.add_argument("--radius_deg", type=float, default=5.0)
    args = ap.parse_args()

    fit_df = pd.read_parquet(args.fit_result)
    fit_df = fit_df[fit_df["success"] == True].copy()
    print(f"Loaded {len(fit_df)} successful fits")

    ref = load_reference(args.reference)
    print(f"Loaded reference geochemistry for {len(ref)} sites: {list(ref.keys())}")

    elements = ["Al", "Fe", "Mg", "Ca", "Ti"]
    rows = []
    scatter_points = {el: {"ref": [], "fit": []} for el in elements}

    for site, info in ref.items():
        matched = match_footprints(fit_df, info["lat"], info["lon"], args.radius_deg)
        n = len(matched)
        if n == 0:
            print(f"  site {site}: 0 footprints within {args.radius_deg} deg — skipped")
            continue
        for el in elements:
            fit_ratio = (matched[f"conc_{el}"] / matched["conc_Si"]).replace([np.inf, -np.inf], np.nan).dropna()
            if len(fit_ratio) == 0:
                continue
            # Median, not mean: a handful of near-zero conc_Si fits blow the
            # ratio up to ~1e18 and would otherwise dominate the average.
            fit_mean = fit_ratio.median()
            ref_ratio = info["ratio_to_Si"][el]
            rows.append({
                "site": site, "element": el, "n_footprints": n,
                "fitted_ratio_to_Si_mean": fit_mean,
                "reference_ratio_to_Si": ref_ratio,
            })
            scatter_points[el]["ref"].append(ref_ratio)
            scatter_points[el]["fit"].append(fit_mean)
        print(f"  site {site}: {n} footprints matched")

    table = pd.DataFrame(rows)
    os.makedirs(args.out_dir, exist_ok=True)
    table.to_csv(os.path.join(args.out_dir, "validation_table.csv"), index=False)
    print(f"\nSaved validation_table.csv ({len(table)} rows)")

    # Stats (Table V) are computed for all 5 elements as before. The FIGURE,
    # however, only plots Al: it is the sole element with a real diagonal
    # trend (r=+0.738); Fe/Mg/Ca/Ti scatter near-flat and add no visual
    # information beyond what Table V's r/RMSE/bias columns already report,
    # so plotting them would just repeat the table in picture form.
    stats_rows = []
    for el in elements:
        ref_v = np.array(scatter_points[el]["ref"])
        fit_v = np.array(scatter_points[el]["fit"])
        if len(ref_v) < 2:
            continue
        r, _ = pearsonr(ref_v, fit_v)
        rmse = np.sqrt(np.mean((ref_v - fit_v) ** 2))
        bias = np.mean(fit_v - ref_v)
        stats_rows.append({"element": el, "n_sites": len(ref_v), "pearson_r": r, "rmse": rmse, "bias": bias})

    plot_el = "Al"
    ref_v = np.array(scatter_points[plot_el]["ref"])
    fit_v = np.array(scatter_points[plot_el]["fit"])
    r, _ = pearsonr(ref_v, fit_v)
    rmse = np.sqrt(np.mean((ref_v - fit_v) ** 2))

    plt.rcParams.update({"font.size": 9})
    fig, ax = plt.subplots(figsize=(3.5, 2.7), dpi=600)
    ax.scatter(ref_v, fit_v, color="tab:blue", s=18)
    lims = [min(ref_v.min(), fit_v.min()), max(ref_v.max(), fit_v.max())]
    ax.plot(lims, lims, "k--", linewidth=1.0, label="1:1")
    ax.set_xlabel(f"Reference {plot_el}/Si", fontsize=8.5)
    ax.set_ylabel(f"Fitted {plot_el}/Si", fontsize=8.5)
    ax.tick_params(labelsize=8)
    ax.set_title(f"{plot_el}/Si validation (r={r:.3f}, RMSE={rmse:.3f}, n={len(ref_v)})", fontsize=8)
    ax.legend(fontsize=7, frameon=True, handlelength=1.5)

    fig.tight_layout()
    out_dir_final = os.path.join(args.out_dir, "final")
    os.makedirs(out_dir_final, exist_ok=True)
    fig_path_pdf = os.path.join(out_dir_final, "fig4_validation_scatter.pdf")
    fig_path_png = os.path.join(out_dir_final, "fig4_validation_scatter.png")
    plt.savefig(fig_path_pdf, bbox_inches="tight")
    plt.savefig(fig_path_png, dpi=600, bbox_inches="tight")
    print(f"Saved {fig_path_pdf} and {fig_path_png}")

    stats_df = pd.DataFrame(stats_rows)
    stats_df.to_csv(os.path.join(args.out_dir, "validation_stats.csv"), index=False)
    print("\nValidation stats (ratio-to-Si space):")
    print(stats_df.to_string(index=False))


if __name__ == "__main__":
    main()
