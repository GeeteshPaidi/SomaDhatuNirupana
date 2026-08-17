"""
Timing wrapper around the existing run_fit.run_all_fits (unmodified) so we
can report computational cost (total wall time, ms/spectrum, workers used)
alongside the fit results, without touching the fitting logic itself.
"""
import argparse
import os
import time

import pandas as pd
from run_fit import run_all_fits

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="../../data/processed/catalogue_full.parquet")
    ap.add_argument("--output_dir", default="../../data/processed")
    ap.add_argument("--max_workers", type=int, default=14)
    args = ap.parse_args()

    n_rows = len(pd.read_parquet(args.input, columns=["exposure"]))
    print(f"About to fit {n_rows} spectra with {args.max_workers} workers")

    t0 = time.perf_counter()
    run_all_fits(
        input_path=args.input,
        output_format="parquet",
        max_workers=args.max_workers,
        output_dir=args.output_dir,
    )
    elapsed = time.perf_counter() - t0

    timing_path = os.path.join(args.output_dir, "fit_timing.txt")
    with open(timing_path, "w") as f:
        f.write(f"n_spectra: {n_rows}\n")
        f.write(f"max_workers: {args.max_workers}\n")
        f.write(f"elapsed_sec: {elapsed:.3f}\n")
        f.write(f"sec_per_spectrum_wallclock: {elapsed / max(n_rows,1):.4f}\n")
        f.write(f"effective_spectra_per_sec: {n_rows / max(elapsed,1e-9):.2f}\n")
    print(f"Total fitting wall time: {elapsed:.1f}s for {n_rows} spectra "
          f"({elapsed/max(n_rows,1)*1000:.1f} ms/spectrum wall-clock, {args.max_workers} workers)")
    print(f"Timing saved to {timing_path}")
