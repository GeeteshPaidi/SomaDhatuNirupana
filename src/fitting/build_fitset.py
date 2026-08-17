"""
Build the fit-set catalogue = union of each element's top-flux spectra.

Rationale (per the user): fitting the full 100% is the safe baseline, but a
file that lands in the bottom half of the flux distribution for EVERY element
carries no strong line for any map, so it can be dropped with no loss to any
element's map. Keeping the union of the per-element top-P% flux sets does
exactly that: every element still keeps all of its own best-signal spectra,
while the globally-faint spectra are removed to cut fitting cost.

Because line fluxes are strongly correlated across elements (a bright exposure
is bright in every line), the union is far smaller than 6*50%; in practice it
is close to the single top-P% set. The exact reduction is reported.
"""
import argparse

import numpy as np
import pandas as pd

# Elements whose maps we care about (K-alpha lines that actually carry signal).
UNION_ELEMENTS = ["Mg", "Al", "Si", "Ca", "Ti", "Fe"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalogue", default="../../data/processed/catalogue_expanded.parquet")
    ap.add_argument("--out", default="../../data/processed/catalogue_fitset.parquet")
    ap.add_argument("--top_pct", type=float, default=50.0,
                    help="Per element, keep rows in the top X%% by flux_<element>")
    args = ap.parse_args()

    df = pd.read_parquet(args.catalogue)
    n = len(df)
    print(f"Full catalogue: {n} dayside spectra")

    keep = np.zeros(n, dtype=bool)
    per_el = {}
    for el in UNION_ELEMENTS:
        col = f"flux_{el}"
        if col not in df.columns:
            print(f"  WARNING: {col} missing, skipping")
            continue
        thr = np.nanpercentile(df[col].values, 100 - args.top_pct)
        mask = df[col].values >= thr
        per_el[el] = int(mask.sum())
        keep |= mask
        print(f"  {el}: top-{args.top_pct:.0f}% flux threshold {thr:.4g} -> {mask.sum()} spectra")

    kept = int(keep.sum())
    print(f"\nUnion keep-set: {kept} / {n} = {kept/n*100:.1f}% "
          f"(dropped {n-kept} spectra that are bottom-{100-args.top_pct:.0f}% in every element)")

    df_keep = df[keep].reset_index(drop=True)
    df_keep.to_parquet(args.out, index=False)
    print(f"Saved fit-set -> {args.out}")


if __name__ == "__main__":
    main()
