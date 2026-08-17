"""
Attach the per-element flux columns (and source_file) from the fit-set
catalogue onto the fit result, so downstream mapping can do per-element
top-flux selection.

run_fit.run_all_fits reads the catalogue and writes results with
multiprocessing.Pool.imap (ORDER PRESERVED), one result row per input row.
So result row i corresponds to catalogue row i and we can join by position.
A length check guards against any mismatch.
"""
import argparse

import pandas as pd

FLUX_COLS = ["flux_O", "flux_Mg", "flux_Al", "flux_Si", "flux_Ca", "flux_Ti", "flux_Fe", "flux_total"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalogue", default="../../data/processed/catalogue_fitset.parquet")
    ap.add_argument("--result", default="../../data/processed/catalogue_fitset_result.parquet")
    ap.add_argument("--out", default="../../data/processed/catalogue_fitset_result_enriched.parquet")
    args = ap.parse_args()

    cat = pd.read_parquet(args.catalogue)
    res = pd.read_parquet(args.result)
    if len(cat) != len(res):
        raise SystemExit(f"Length mismatch: catalogue {len(cat)} vs result {len(res)} "
                         f"- cannot positionally join. Re-check the fit input.")

    for c in FLUX_COLS + ["source_file"]:
        if c in cat.columns:
            res[c] = cat[c].values

    res.to_parquet(args.out, index=False)
    print(f"Enriched result ({len(res)} rows, +{len(FLUX_COLS)} flux cols) -> {args.out}")


if __name__ == "__main__":
    main()
