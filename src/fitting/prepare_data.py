"""
Raw CLASS FITS -> preprocessed spectrum catalogue (Sec 2 of the paper).

Each raw FITS file (data/raw/2026/<month>/<day>/*.fits) contains a single
8-16s exposure: a 2048-channel counts spectrum plus footprint geometry
(V0-V3 lat/lon corners) and SOLARANG (sun-surface-normal angle). Nightside
exposures (SOLARANG >= 90 deg) receive negligible solar-driven XRF signal
and are used here as the empirical background reference, exactly as
described in Sec 2.2 of the paper. There is no separate background FITS
file bundled with this data release, so the background spectrum is built
from the nightside population of the same batch rather than a dedicated
background-allevents.fits.

Pipeline per dayside spectrum:
  1. Convert counts -> count rate (counts / EXPOSURE).
  2. Subtract the mean nightside count-rate spectrum (background subtraction).
  3. Clip negative residuals to zero (non-negativity constraint).
  4. Denoise: Stationary Wavelet Transform (sym8, level 2, soft threshold)
     followed by Savitzky-Golay smoothing (window 7, order 3).

Output: a parquet catalogue with columns [counts, footprint, mid_utc,
sat_lat, sat_lon, solarang, exposure] ready for src/fitting/run_fit.py.
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

N_CHANNELS = 2048


def read_fits_spectrum(path):
    with fits.open(path) as hdul:
        hdr = hdul[1].header
        data = hdul[1].data
        counts = np.zeros(N_CHANNELS, dtype=np.float64)
        ch = data["CHANNEL"].astype(int)
        ch = np.clip(ch, 0, N_CHANNELS - 1)
        counts[ch] = data["COUNTS"]

        footprint = [
            (hdr["V0_LAT"], hdr["V0_LON"]),
            (hdr["V1_LAT"], hdr["V1_LON"]),
            (hdr["V2_LAT"], hdr["V2_LON"]),
            (hdr["V3_LAT"], hdr["V3_LON"]),
        ]
        meta = {
            "counts": counts,
            "footprint": footprint,
            "exposure": float(hdr["EXPOSURE"]),
            "solarang": float(hdr["SOLARANG"]),
            "sat_lat": float(hdr["SAT_LAT"]),
            "sat_lon": float(hdr["SAT_LON"]),
            "mid_utc": str(hdr["MID_UTC"]),
        }
        return meta


def wavelet_savgol_denoise(rate, wavelet="sym8", level=2, sg_window=7, sg_order=3):
    n = len(rate)
    coeffs = pywt.swt(rate, wavelet, level=level)
    # coeffs: list of (cA, cD) per level, coarsest first in pywt.swt output order
    denoised_coeffs = []
    for cA, cD in coeffs:
        sigma = np.median(np.abs(cD - np.median(cD))) / 0.6745 if len(cD) else 0.0
        thresh = sigma * np.sqrt(2 * np.log(max(n, 2)))
        cD_soft = np.sign(cD) * np.maximum(np.abs(cD) - thresh, 0.0)
        denoised_coeffs.append((cA, cD_soft))
    rec = pywt.iswt(denoised_coeffs, wavelet)
    rec = rec[:n]
    if n >= sg_window:
        rec = savgol_filter(rec, window_length=sg_window, polyorder=sg_order)
    return rec


def build_background(fits_paths, solarang_threshold=90.0, max_bg_files=2000):
    """Average count-rate spectrum over nightside exposures."""
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


def process_batch(fits_paths, background_rate, solarang_threshold=90.0):
    rows = []
    example_raw = example_bgsub = example_denoised = None
    t0 = time.perf_counter()
    n_fit_attempted = 0
    for p in tqdm(fits_paths, desc="Preprocessing dayside spectra"):
        try:
            meta = read_fits_spectrum(p)
        except Exception:
            continue
        if meta["solarang"] >= solarang_threshold:
            continue  # nightside, not a science spectrum for fitting

        rate = meta["counts"] / meta["exposure"]
        bgsub = np.clip(rate - background_rate, 0.0, None)
        denoised = wavelet_savgol_denoise(bgsub)
        denoised = np.clip(denoised, 0.0, None)

        if example_raw is None:
            example_raw, example_bgsub, example_denoised = rate.copy(), bgsub.copy(), denoised.copy()

        rows.append({
            "counts": denoised,
            "footprint": meta["footprint"],
            "exposure": meta["exposure"],
            "solarang": meta["solarang"],
            "sat_lat": meta["sat_lat"],
            "sat_lon": meta["sat_lon"],
            "mid_utc": meta["mid_utc"],
            "source_file": os.path.basename(p),
        })
        n_fit_attempted += 1

    elapsed = time.perf_counter() - t0
    df = pd.DataFrame(rows)
    timing = {
        "n_spectra": n_fit_attempted,
        "elapsed_sec": elapsed,
        "sec_per_spectrum": elapsed / max(n_fit_attempted, 1),
    }
    example = {"raw": example_raw, "bgsub": example_bgsub, "denoised": example_denoised}
    return df, timing, example


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw_dir", default="data/raw/2026")
    ap.add_argument("--out", default="data/processed/catalogue.parquet")
    ap.add_argument("--example_out", default="outputs/figures/preprocessing_example.npz")
    ap.add_argument("--limit", type=int, default=None, help="cap number of FITS files (debug)")
    args = ap.parse_args()

    fits_paths = sorted(glob.glob(os.path.join(args.raw_dir, "**", "*.fits"), recursive=True))
    if args.limit:
        fits_paths = fits_paths[: args.limit]
    print(f"Found {len(fits_paths)} raw FITS files under {args.raw_dir}")

    background_rate, n_bg = build_background(fits_paths)
    print(f"Background built from {n_bg} nightside exposures, total rate = {background_rate.sum():.4f} counts/s")

    df, timing, example = process_batch(fits_paths, background_rate)
    print(f"Preprocessed {timing['n_spectra']} dayside spectra in {timing['elapsed_sec']:.2f}s "
          f"({timing['sec_per_spectrum']*1000:.2f} ms/spectrum)")

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df.to_parquet(args.out, index=False)
    print(f"Saved catalogue to {args.out}")

    os.makedirs(os.path.dirname(args.example_out), exist_ok=True)
    np.savez(args.example_out, **example, background=background_rate,
             energy_factor=0.0135, n_channels=N_CHANNELS)
    print(f"Saved Fig.1 example arrays to {args.example_out}")

    timing_path = os.path.splitext(args.out)[0] + "_timing.txt"
    with open(timing_path, "w") as f:
        for k, v in timing.items():
            f.write(f"{k}: {v}\n")
    print(f"Saved timing info to {timing_path}")


if __name__ == "__main__":
    main()
