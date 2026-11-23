"""
R1-5: dependence of the resolution term sigma0 (and gain/offset) on detector
temperature within one month.

Dayside 16-SCD spectra of the month are split into detector-temperature (TEMP
header) quantile bins; each bin is stacked, the nightside particle background
subtracted, and calibrated with the same line fit as run_month.calibrate.

Usage: python src/asr/sigma_temp.py --month 2024/12 [--nbins 5]
"""
import argparse
import glob
import json
import os
import sys
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_month as rm  # noqa: E402

ROOT = rm.ROOT


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--month", default="2024/12")
    ap.add_argument("--nbins", type=int, default=5)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    tag = a.month.replace("/", "_")
    paths = sorted(glob.glob(os.path.join(ROOT, "data", "raw", *a.month.split("/"), "**", "*.fits"), recursive=True))
    with Pool(a.workers) as p:
        rows = [r for r in p.imap(rm.read_one, paths, chunksize=256) if r is not None]
    bkg, _, _ = rm.background_templates(rows)
    day = [r for r in rows if r["sat_alt"] > 0 and r["exposure"] > 0 and r["n_scd"] == 16 and r["solarang"] < 90]
    T = np.array([r["temp"] for r in day])
    edges = np.quantile(T, np.linspace(0, 1, a.nbins + 1))
    out = {"month": tag, "n_day16": len(day), "bins": []}
    allc = rm.calibrate(np.sum([r["counts"] for r in day], axis=0) - bkg * sum(r["exposure"] for r in day))
    out["all"] = allc
    print(f"{tag}: {len(day)} dayside 16-SCD spectra, TEMP {T.min():.1f} .. {T.max():.1f} C")
    print(f"  all      : FWHM(1.49) {allc['fwhm_1.49keV_eV']:.1f} eV  gain {allc['gain']:.5f}  offset {allc['offset']*1e3:+.1f} eV")
    for k in range(a.nbins):
        sel = (T >= edges[k]) & ((T < edges[k + 1]) if k < a.nbins - 1 else (T <= edges[k + 1]))
        rr = [r for r, s in zip(day, sel) if s]
        net = np.sum([r["counts"] for r in rr], axis=0) - bkg * sum(r["exposure"] for r in rr)
        c = rm.calibrate(net)
        c.update(t_lo=float(edges[k]), t_hi=float(edges[k + 1]), t_med=float(np.median(T[sel])), n=len(rr))
        out["bins"].append(c)
        print(f"  {edges[k]:6.1f}..{edges[k+1]:6.1f} C (n={len(rr):5d}): FWHM(1.49) {c['fwhm_1.49keV_eV']:.1f} eV  "
              f"gain {c['gain']:.5f}  offset {c['offset']*1e3:+.1f} eV  ok={c['success']}")
    os.makedirs(os.path.join(ROOT, "outputs", "asr", "sigma_temp"), exist_ok=True)
    with open(os.path.join(ROOT, "outputs", "asr", "sigma_temp", f"{tag}.json"), "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
