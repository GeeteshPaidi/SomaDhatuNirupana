"""
Coverage figure (Sec "Spatial Coverage"): binary occupancy of the
1,056,619 dayside footprints on a global 1-deg x 1-deg grid.

Reproduces the coverage computation in src/mapping/map_making_v1.py's
GaussianArray.check_coverage() (each footprint's quadrilateral is rasterised
into the grid; a cell is "covered" if any footprint touched it), applied to
ALL dayside footprints from data/processed/months/*.parquet (not just the
flux-selected fitted subset used for the abundance maps), matching what the
paper reports: "Projecting the 1,056,619 dayside footprints onto a global
grid yields 99.0% coverage at 1x1 resolution."
"""
import glob
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "mapping"))
from map_making_v1 import GaussianArray

BLOCK_LAT = [-90, 90, 90, -90]
BLOCK_LON = [-180, -180, 180, 180]


def main():
    grid = GaussianArray(grid_size=(180, 360))
    n_total = 0
    for f in sorted(glob.glob("../../data/processed/months/*.parquet")):
        df = pd.read_parquet(f, columns=["footprint"])
        for fp in df["footprint"]:
            img_lat = [p[0] for p in fp]
            img_lon = [((p[1] + 180) % 360) - 180 for p in fp]
            if max(img_lon) - min(img_lon) > 180:
                continue
            grid.add_gaussian_box(img_lat, img_lon, BLOCK_LAT, BLOCK_LON, 1.0)
            n_total += 1
        print(f"  {os.path.basename(f)}: {len(df)} footprints (running total {n_total})")

    coverage = grid.check_coverage()
    print(f"\nTotal footprints: {n_total}, coverage: {coverage*100:.2f}%")

    # The 99.0% statistic is a binary occupied/unoccupied count over the
    # same grid.arr[:, :, 1] count array used below -- at that saturation
    # level a plain boolean occ map is ~99% one colour and shows essentially
    # no visible structure. Plotting the actual per-cell footprint COUNT
    # (same underlying accumulator, just not thresholded to a boolean) shows
    # the real orbital-track density variation, and the empty cells that
    # drive the 1.0% gap are still exactly the cells with count == 0 -- the
    # single darkest colour on the scale, unambiguous with a colorbar.
    counts = grid.arr[:, :, 1]
    n_empty = int(np.sum(counts == 0))
    print(f"Empty cells: {n_empty} / {counts.size} ({100 * n_empty / counts.size:.2f}%)")
    log_counts = np.log10(counts + 1)

    plt.rcParams.update({"font.size": 9})
    fig, ax = plt.subplots(figsize=(3.5, 2.3), dpi=600)

    im = ax.imshow(log_counts, extent=[-180, 180, -90, 90], origin="lower",
                    cmap="viridis", aspect="equal", interpolation="nearest")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xticks(range(-150, 151, 75))
    ax.set_yticks(range(-90, 91, 45))
    ax.tick_params(labelsize=7)
    ax.set_xlabel("Longitude (deg)", fontsize=8)
    ax.set_ylabel("Latitude (deg)", fontsize=8)
    # No in-axes title: the LaTeX \caption carries the description (matches
    # the style used for Fig. 1), and a full title line overflows a
    # 3.5in single-column canvas at a readable font size.

    cbar = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.06, pad=0.32, aspect=30)
    cbar.set_label(r"$\log_{10}$(footprints per $1^\circ\times1^\circ$ cell $+\,1$)", fontsize=7.5)
    cbar.ax.tick_params(labelsize=7)

    fig.tight_layout()
    os.makedirs("../../outputs/figures/supplementary", exist_ok=True)
    fig.savefig("../../outputs/figures/supplementary/supp_coverage.pdf", bbox_inches="tight")
    fig.savefig("../../outputs/figures/supplementary/supp_coverage.png", dpi=600, bbox_inches="tight")
    print("Saved supp_coverage.pdf/.png to outputs/figures/supplementary/")


if __name__ == "__main__":
    main()
