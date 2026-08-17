"""
Parallel raw CLASS FITS -> preprocessed spectrum catalogue (Sec 2).

Identical preprocessing to prepare_data.py (nightside background subtraction,
non-negativity clip, SWT sym8/level-2 soft-threshold denoise + Savitzky-Golay
window-7/order-3), but:

  * multiprocessing over dayside spectra (the wavelet+SavGol step is CPU-bound
    and embarrassingly parallel) -> ~14x faster on this machine, so ~900k
    spectra prepare in ~20 min instead of ~5 h;
  * adds per-element K-alpha-window integrated flux columns
    (flux_O/Mg/Al/Si/Ca/Ti/Fe and flux_total) so a "keep only the strongest
    signal spectra" quality filter (the user's top-50%-Al-flux catalogue,
    e.g. catalogue_Al_50_*.parquet) can be applied before fitting.

Background is built from a nightside sample (SOLARANG >= 90) exactly like the
serial version's build_background (bounded scan), then held constant for the
parallel pass.
"""
import argparse
import glob
import os
import time

import numpy as np
import pandas as pd
import pywt
from astropy.io import fits
from scipy.signal import savgol_filter
from tqdm import tqdm
from multiprocessing import Pool

N_CHANNELS = 2048
ENERGY_FACTOR = 0.0135  # keV per channel (as used in Fig.1 / paper)

# K-alpha line energies (keV) -> channel = E / ENERGY_FACTOR
KALPHA_KEV = {
    "O": 0.525, "Mg": 1.254, "Al": 1.487, "Si": 1.740,
    "Ca": 3.692, "Ti": 4.511, "Fe": 6.404,
}
FLUX_HALF_WIN = 6  # +/- channels around each K-alpha centroid

_BACKGROUND = None  # set in each worker via initializer


def read_fits_spectrum(path):
    with fits.open(path) as hdul:
        hdr = hdul[1].header
        data = hdul[1].data
        counts = np.zeros(N_CHANNELS, dtype=np.float64)
        ch = np.clip(data["CHANNEL"].astype(int), 0, N_CHANNELS - 1)
        counts[ch] = data["COUNTS"]
        footprint = [(hdr[f"V{k}_LAT"], hdr[f"V{k}_LON"]) for k in range(4)]
        return {
            "counts": counts,
            "footprint": footprint,
            "exposure": float(hdr["EXPOSURE"]),
            "solarang": float(hdr["SOLARANG"]),
            "sat_lat": float(hdr["SAT_LAT"]),
            "sat_lon": float(hdr["SAT_LON"]),
            "mid_utc": str(hdr["MID_UTC"]),
        }


def wavelet_savgol_denoise(rate, wavelet="sym8", level=2, sg_window=7, sg_order=3):
    n = len(rate)
    coeffs = pywt.swt(rate, wavelet, level=level)
    out = []
    for cA, cD in coeffs:
        sigma = np.median(np.abs(cD - np.median(cD))) / 0.6745 if len(cD) else 0.0
        thresh = sigma * np.sqrt(2 * np.log(max(n, 2)))
        cD_soft = np.sign(cD) * np.maximum(np.abs(cD) - thresh, 0.0)
        out.append((cA, cD_soft))
    rec = pywt.iswt(out, wavelet)[:n]
    if n >= sg_window:
        rec = savgol_filter(rec, window_length=sg_window, polyorder=sg_order)
    return rec


def _flux_windows():
    win = {}
    for el, e in KALPHA_KEV.items():
        c = int(round(e / ENERGY_FACTOR))
        lo, hi = max(0, c - FLUX_HALF_WIN), min(N_CHANNELS, c + FLUX_HALF_WIN + 1)
        win[el] = (lo, hi)
    return win


FLUX_WIN = _flux_windows()


def build_background(fits_paths, solarang_threshold=90.0, max_bg_files=3000):
    bg_sum = np.zeros(N_CHANNELS, dtype=np.float64)
    n_bg = 0
    for p in tqdm(fits_paths[:max_bg_files], desc="Scanning for nightside background"):
        try:
            meta = read_fits_spectrum(p)
        except Exception:
            continue
        if meta["solarang"] >= solarang_threshold:
            bg_sum += meta["counts"] / meta["exposure"]
            n_bg += 1
    if n_bg == 0:
        return np.zeros(N_CHANNELS), 0
    return bg_sum / n_bg, n_bg


def _init_worker(background):
    global _BACKGROUND
    _BACKGROUND = background


def _process_one(path):
    try:
        meta = read_fits_spectrum(path)
    except Exception:
        return None
    if meta["solarang"] >= 90.0:
        return None  # nightside
    rate = meta["counts"] / meta["exposure"]
    bgsub = np.clip(rate - _BACKGROUND, 0.0, None)
    denoised = np.clip(wavelet_savgol_denoise(bgsub), 0.0, None)

    row = {
        # float32 halves memory/disk vs float64; fit precision is unaffected
        # (counts are already denoised, low dynamic range).
        "counts": denoised.astype(np.float32),
        "footprint": meta["footprint"],
        "exposure": meta["exposure"],
        "solarang": meta["solarang"],
        "sat_lat": meta["sat_lat"],
        "sat_lon": meta["sat_lon"],
        "mid_utc": meta["mid_utc"],
        "source_file": os.path.basename(path),
    }
    tot = 0.0
    for el, (lo, hi) in FLUX_WIN.items():
        f = float(denoised[lo:hi].sum())
        row[f"flux_{el}"] = f
        tot += f
    row["flux_total"] = float(denoised.sum())
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw_dir", default="data/raw")
    ap.add_argument("--out", default="data/processed/catalogue.parquet")
    ap.add_argument("--example_out", default="outputs/figures/preprocessing_example.npz")
    ap.add_argument("--workers", type=int, default=14)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    fits_paths = sorted(glob.glob(os.path.join(args.raw_dir, "**", "*.fits"), recursive=True))
    if args.limit:
        fits_paths = fits_paths[: args.limit]
    print(f"Found {len(fits_paths)} raw FITS files under {args.raw_dir}", flush=True)

    background_rate, n_bg = build_background(fits_paths)
    print(f"Background from {n_bg} nightside exposures (rate sum {background_rate.sum():.4f})", flush=True)

    t0 = time.perf_counter()
    rows = []
    with Pool(processes=args.workers, initializer=_init_worker, initargs=(background_rate,)) as pool:
        for r in tqdm(pool.imap_unordered(_process_one, fits_paths, chunksize=64),
                      total=len(fits_paths), desc="Preprocessing (parallel)"):
            if r is not None:
                rows.append(r)
    elapsed = time.perf_counter() - t0

    df = pd.DataFrame(rows)
    print(f"Preprocessed {len(df)} dayside spectra in {elapsed:.1f}s "
          f"({elapsed/max(len(df),1)*1000:.2f} ms/spectrum, {args.workers} workers)", flush=True)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df.to_parquet(args.out, index=False)
    print(f"Saved catalogue -> {args.out}", flush=True)

    # Example spectrum for Fig.1: highest-flux clean dayside spectrum (re-read raw/bgsub).
    if len(df):
        top = df.sort_values("flux_total", ascending=False).iloc[0]
        # locate its full path
        cand = glob.glob(os.path.join(args.raw_dir, "**", top["source_file"]), recursive=True)
        if cand:
            meta = read_fits_spectrum(cand[0])
            rate = meta["counts"] / meta["exposure"]
            bgsub = np.clip(rate - background_rate, 0.0, None)
            denoised = np.clip(wavelet_savgol_denoise(bgsub), 0.0, None)
            os.makedirs(os.path.dirname(args.example_out), exist_ok=True)
            np.savez(args.example_out, raw=rate, bgsub=bgsub, denoised=denoised,
                     background=background_rate, energy_factor=ENERGY_FACTOR,
                     n_channels=N_CHANNELS, source_file=top["source_file"])
            print(f"Saved Fig.1 example ({top['source_file']}) -> {args.example_out}", flush=True)

    timing_path = os.path.splitext(args.out)[0] + "_timing.txt"
    with open(timing_path, "w") as f:
        f.write(f"n_spectra: {len(df)}\n")
        f.write(f"elapsed_sec: {elapsed:.3f}\n")
        f.write(f"workers: {args.workers}\n")
        f.write(f"sec_per_spectrum_wallclock: {elapsed/max(len(df),1):.6f}\n")
    print(f"Saved timing -> {timing_path}", flush=True)


if __name__ == "__main__":
    main()
