"""
Regenerate Table 5 and Figures 1, 4 and 7 of the paper from a clean checkout, using only the
archived products in products/v5 (no CLASS Level-1 data or per-spectrum fit tables needed).

    python src/asr/reproduce_release.py

Steps:
1. copy products/v5/* into outputs/asr/ (files already present are left alone);
2. derive the epoch-A mapping set (Nov 2024 - May 2026) from the primary mapping set;
3. Table 5: new_sites.py on the primary mapping set, then the table writer of make_numbers.py;
4. Figure 1: the frozen first-submission figure (products/figures_frozen/fig_background.png);
5. Figure 4: maps from maps_2deg_v5.npz with LOLA contours (data/lola/ldem_4.img);
6. Figure 7: epoch_compare.py on maps_2deg_v5_eA/eB.npz.
"""
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.join(ROOT, "outputs", "asr")
sys.path.insert(0, HERE)


def run(*args):
    print("+", " ".join(args), flush=True)
    subprocess.run([sys.executable, *args], check=True, cwd=ROOT)


def main():
    src = os.path.join(ROOT, "products", "v5")
    for d, _, fs in os.walk(src):
        for f in fs:
            dst = os.path.join(OUT, os.path.relpath(os.path.join(d, f), src))
            if not os.path.exists(dst):
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                shutil.copy2(os.path.join(d, f), dst)

    import pandas as pd
    from maps import parse_months
    ea = os.path.join(OUT, "maps", "mapping_set_v5_eA.parquet")
    if not os.path.exists(ea):
        ms = pd.read_parquet(os.path.join(OUT, "maps", "mapping_set_v5.parquet"))
        ms[ms.month.map(parse_months("2024_11:2026_05"))].reset_index(drop=True).to_parquet(ea)

    # Table 5
    run(os.path.join(HERE, "new_sites.py"), "--mapping_set", os.path.join(OUT, "maps", "mapping_set_v5.parquet"),
        "--out", os.path.join(ROOT, "revision", "new_sites_v5.csv"))
    import make_numbers
    make_numbers.new_sites_numbers()

    # Figure 7 input
    run(os.path.join(HERE, "epoch_compare.py"))

    # Figures 1, 4, 7
    import paper_figures as pf
    pf.fig_background()
    pf.fig_maps()
    pf.fig_epochs()
    print("Table 5 -> paper_asr/tables/new_sites.tex; Figs 1, 4, 7 -> paper_asr/figures/")


if __name__ == "__main__":
    main()
