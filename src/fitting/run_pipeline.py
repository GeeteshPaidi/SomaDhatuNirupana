"""
Memory-bounded, restartable end-to-end pipeline for the expanded dataset
(~1.5M raw FITS across ~24 months). Processes ONE MONTH AT A TIME so the
machine never holds more than one month of spectra in RAM.

Stages
------
LOOP 1 (prepare):  for each month  data/raw/YYYY/MM
    - build nightside background from that month's sample
    - parallel denoise all dayside spectra + compute per-element K-alpha flux
    - write a float32 month catalogue to data/processed/months/cat_YYYY_MM.parquet
    - accumulate only the small flux columns in memory (a few floats/spectrum)

THRESHOLDS:  global per-element top-P% flux thresholds from the accumulated
    flux; a spectrum is kept if it is in the top P% for AT LEAST ONE element
    (union). Files that are bottom-(100-P)% in every element are dropped -
    they carry no strong line for any map, so no map loses anything.

LOOP 2 (fit):  for each month
    - load its catalogue, keep only union-top-P% rows
    - parallel FP/TRF fit (reuses run_fit.fit_single_row)
    - carry the flux columns onto the result, append to the master result

Output: data/processed/catalogue_expanded_result.parquet  (footprint, success,
cost, conc_*, scale, std_dev, flux_*, source_file) - small, ready for mapping.
Also writes a coverage + cost report and reuses/refreshes the Fig.1 example.
"""
import argparse
import glob
import os
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd

import prepare_data_parallel as prep
from run_fit import fit_single_row

# Silence the per-month nightside-scan progress bar (24 months x 3000 files
# of tqdm carriage-return spam would bury the real log).
prep.tqdm = lambda x, **k: x

UNION_ELEMENTS = ["Mg", "Al", "Si", "Ca", "Ti", "Fe"]
FLUX_COLS = ["flux_O", "flux_Mg", "flux_Al", "flux_Si", "flux_Ca", "flux_Ti", "flux_Fe", "flux_total"]


def list_months(raw_dir):
    months = []
    for y in sorted(glob.glob(os.path.join(raw_dir, "*"))):
        if not os.path.isdir(y):
            continue
        for m in sorted(glob.glob(os.path.join(y, "*"))):
            if os.path.isdir(m):
                months.append(m)
    return months


def prepare_month(month_dir, workers):
    paths = sorted(glob.glob(os.path.join(month_dir, "**", "*.fits"), recursive=True))
    if not paths:
        return None, None
    background, _ = prep.build_background(paths, max_bg_files=3000)
    rows = []
    with Pool(workers, initializer=prep._init_worker, initargs=(background,)) as pool:
        for r in pool.imap_unordered(prep._process_one, paths, chunksize=64):
            if r is not None:
                rows.append(r)
    return pd.DataFrame(rows), background


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw_dir", default="../../data/raw")
    ap.add_argument("--months_dir", default="../../data/processed/months")
    ap.add_argument("--out", default="../../data/processed/catalogue_expanded_result.parquet")
    ap.add_argument("--example_out", default="../../outputs/figures/preprocessing_example.npz")
    ap.add_argument("--workers", type=int, default=14)
    ap.add_argument("--top_pct", type=float, default=50.0)
    ap.add_argument("--skip_prepare", action="store_true",
                    help="reuse existing month catalogues in months_dir")
    ap.add_argument("--prepare_only", action="store_true",
                    help="stop after LOOP1 + thresholds + coverage (checkpoint before the long fit)")
    args = ap.parse_args()

    os.makedirs(args.months_dir, exist_ok=True)
    months = list_months(args.raw_dir)
    print(f"Found {len(months)} months under {args.raw_dir}", flush=True)

    # ---------- LOOP 1: prepare per month ----------
    flux_frames = []
    coverage_hit = np.zeros((180, 360), dtype=bool)
    best_example = None  # (flux_total, source_file, month_dir, background)
    t_prep0 = time.perf_counter()
    for i, mdir in enumerate(months, 1):
        tag = "_".join(mdir.replace("\\", "/").split("/")[-2:])
        cat_path = os.path.join(args.months_dir, f"cat_{tag}.parquet")
        if args.skip_prepare and os.path.exists(cat_path):
            df = pd.read_parquet(cat_path, columns=["footprint"] + FLUX_COLS + ["source_file"])
            print(f"[{i}/{len(months)}] {tag}: reuse {len(df)} spectra", flush=True)
        else:
            df, background = prepare_month(mdir, args.workers)
            if df is None or len(df) == 0:
                print(f"[{i}/{len(months)}] {tag}: no dayside spectra", flush=True)
                continue
            df.to_parquet(cat_path, index=False)
            print(f"[{i}/{len(months)}] {tag}: {len(df)} spectra -> {cat_path}", flush=True)
            # track best example spectrum for Fig.1
            top = df.loc[df["flux_total"].idxmax()]
            if best_example is None or top["flux_total"] > best_example[0]:
                best_example = (float(top["flux_total"]), top["source_file"], mdir, background)

        flux_frames.append(df[["source_file"] + FLUX_COLS])
        # coverage from footprint centers
        for fp in df["footprint"].values:
            clat = np.mean([p[0] for p in fp])
            clon = ((np.mean([p[1] for p in fp]) + 180) % 360) - 180
            r = int(np.clip((clat + 90) / 180 * 180, 0, 179))
            c = int(np.clip((clon + 180) / 360 * 360, 0, 359))
            coverage_hit[r, c] = True

    flux_all = pd.concat(flux_frames, ignore_index=True)
    n_total = len(flux_all)
    cov_pct = coverage_hit.sum() / (180 * 360) * 100
    print(f"\nLOOP1 done in {time.perf_counter()-t_prep0:.0f}s: {n_total} dayside spectra, "
          f"coverage {cov_pct:.2f}% of 1-deg grid", flush=True)

    # ---------- thresholds: union of per-element top-P% ----------
    thresholds = {}
    keep_masks = {}
    keep_union = np.zeros(n_total, dtype=bool)
    for el in UNION_ELEMENTS:
        col = f"flux_{el}"
        thr = np.nanpercentile(flux_all[col].values, 100 - args.top_pct)
        thresholds[el] = thr
        m = flux_all[col].values >= thr
        keep_masks[el] = m
        keep_union |= m
        print(f"  {el}: top-{args.top_pct:.0f}% flux>={thr:.4g} -> {m.sum()} spectra", flush=True)
    n_keep = int(keep_union.sum())
    print(f"Union keep-set: {n_keep}/{n_total} = {n_keep/n_total*100:.1f}% "
          f"(dropped {n_total-n_keep} globally-faint spectra)", flush=True)
    np.save(os.path.join(args.months_dir, "coverage_hit.npy"), coverage_hit)

    # ---------- Fig.1 example (highest-flux clean spectrum) ----------
    if best_example is not None:
        _, src, mdir, background = best_example
        cand = glob.glob(os.path.join(mdir, "**", src), recursive=True)
        if cand:
            meta = prep.read_fits_spectrum(cand[0])
            rate = meta["counts"] / meta["exposure"]
            bgsub = np.clip(rate - background, 0.0, None)
            den = np.clip(prep.wavelet_savgol_denoise(bgsub), 0.0, None)
            os.makedirs(os.path.dirname(args.example_out), exist_ok=True)
            np.savez(args.example_out, raw=rate, bgsub=bgsub, denoised=den,
                     background=background, energy_factor=prep.ENERGY_FACTOR,
                     n_channels=prep.N_CHANNELS, source_file=src)
            print(f"Saved Fig.1 example ({src}) -> {args.example_out}", flush=True)

    if args.prepare_only:
        print("prepare_only: stopping before fit (LOOP2). "
              f"Re-run with --skip_prepare to fit the {n_keep} kept spectra.", flush=True)
        return

    # ---------- LOOP 2: fit per month (kept rows only) ----------
    result_frames = []
    n_fit_total = 0
    t_fit0 = time.perf_counter()
    for i, mdir in enumerate(months, 1):
        tag = "_".join(mdir.replace("\\", "/").split("/")[-2:])
        cat_path = os.path.join(args.months_dir, f"cat_{tag}.parquet")
        if not os.path.exists(cat_path):
            continue
        df = pd.read_parquet(cat_path)
        # per-element union keep within this month
        mask = np.zeros(len(df), dtype=bool)
        for el in UNION_ELEMENTS:
            mask |= df[f"flux_{el}"].values >= thresholds[el]
        dfk = df[mask].reset_index(drop=True)
        if len(dfk) == 0:
            continue
        args_list = [(dfk["counts"].iloc[j], dfk["footprint"].iloc[j]) for j in range(len(dfk))]
        with Pool(args.workers) as pool:
            results = pool.map(fit_single_row, args_list, chunksize=32)
        res = pd.DataFrame(results)
        for c in FLUX_COLS + ["source_file"]:
            res[c] = dfk[c].values
        result_frames.append(res)
        n_fit_total += len(res)
        print(f"[fit {i}/{len(months)}] {tag}: fit {len(dfk)} (cum {n_fit_total}), "
              f"{time.perf_counter()-t_fit0:.0f}s", flush=True)

    result = pd.concat(result_frames, ignore_index=True)
    fit_elapsed = time.perf_counter() - t_fit0
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    result.to_parquet(args.out, index=False)
    ok = int(result["success"].sum())
    print(f"\nFITTED {len(result)} spectra ({ok} ok) in {fit_elapsed:.0f}s "
          f"({fit_elapsed/max(len(result),1)*1000:.1f} ms/spectrum wall, {args.workers} workers)", flush=True)
    print(f"Saved result -> {args.out}", flush=True)

    # timing + coverage report
    with open(os.path.splitext(args.out)[0] + "_timing.txt", "w") as f:
        f.write(f"n_total_dayside: {n_total}\n")
        f.write(f"n_fit (union top-{args.top_pct:.0f}%): {len(result)}\n")
        f.write(f"coverage_pct_1deg: {cov_pct:.3f}\n")
        f.write(f"fit_elapsed_sec: {fit_elapsed:.1f}\n")
        f.write(f"fit_ms_per_spectrum_wall: {fit_elapsed/max(len(result),1)*1000:.3f}\n")
        f.write(f"workers: {args.workers}\n")


if __name__ == "__main__":
    main()
