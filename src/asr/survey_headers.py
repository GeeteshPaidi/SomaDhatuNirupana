"""
Task 0.3: survey CLASS L1 FITS headers across all months.

Samples up to N files per month (evenly spaced in sorted order) and tabulates
the keywords relevant to the forward model: GAIN, DETCHANS, EXPOSURE, TEMP,
SOLARANG (incidence), EMISNANG (emergence), PHASEANG, SCD_USED, SW_VERSN,
RESPFILE/ANCRFILE. Also counts total FITS files per month (full file listing).

Writes outputs/asr/header_survey.csv (one row per sampled file) and
outputs/asr/header_survey_summary.txt.
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd
from astropy.io import fits

KEYS = ["GAIN", "DETCHANS", "EXPOSURE", "TEMP", "SOLARANG", "EMISNANG", "PHASEANG",
        "SCD_USED", "SW_VERSN", "RESPFILE", "ANCRFILE", "SAT_ALT", "MID_UTC"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw_dir", default="data/raw")
    ap.add_argument("--per_month", type=int, default=200)
    ap.add_argument("--out_dir", default="outputs/asr")
    args = ap.parse_args()

    rows, counts = [], {}
    for mdir in sorted(glob.glob(os.path.join(args.raw_dir, "*", "*"))):
        tag = "_".join(mdir.replace("\\", "/").split("/")[-2:])
        paths = sorted(glob.glob(os.path.join(mdir, "**", "*.fits"), recursive=True))
        counts[tag] = len(paths)
        if not paths:
            continue
        idx = np.unique(np.linspace(0, len(paths) - 1, min(args.per_month, len(paths))).astype(int))
        for i in idx:
            with fits.open(paths[i]) as h:
                hdr, d = h[1].header, h[1].data
                r = {k: hdr.get(k) for k in KEYS}
                r["month"] = tag
                r["nrows"] = len(d)
                r["total_counts"] = float(d["COUNTS"].sum())
                r["n_scd"] = len([s for s in str(hdr.get("SCD_USED", "")).split(",") if s.strip()])
                rows.append(r)
    df = pd.DataFrame(rows)
    os.makedirs(args.out_dir, exist_ok=True)
    df.to_csv(os.path.join(args.out_dir, "header_survey.csv"), index=False)

    lines = [f"Total FITS files: {sum(counts.values())}",
             "Files per month: " + ", ".join(f"{k}={v}" for k, v in counts.items()),
             f"Sampled: {len(df)}"]
    for k in ["GAIN", "DETCHANS", "EXPOSURE", "SW_VERSN", "n_scd", "nrows", "RESPFILE", "ANCRFILE"]:
        lines.append(f"{k}: {df[k].value_counts(dropna=False).to_dict()}")
    for k in ["TEMP", "SOLARANG", "EMISNANG", "PHASEANG", "SAT_ALT", "total_counts"]:
        v = pd.to_numeric(df[k], errors="coerce")
        lines.append(f"{k}: min={v.min():.4g} median={v.median():.4g} max={v.max():.4g}")
    lines.append("SOLARANG==PHASEANG fraction: "
                 f"{np.mean(np.isclose(df.SOLARANG.astype(float), df.PHASEANG.astype(float))):.3f}")
    lines.append("Per-month: " + df.groupby("month").agg(
        sw=("SW_VERSN", lambda s: sorted(set(s))), nscd=("n_scd", "median"),
        temp=("TEMP", "median")).to_string())
    txt = "\n".join(lines)
    with open(os.path.join(args.out_dir, "header_survey_summary.txt"), "w") as f:
        f.write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
